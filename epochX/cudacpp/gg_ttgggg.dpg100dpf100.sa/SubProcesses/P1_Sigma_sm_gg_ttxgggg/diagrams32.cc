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
  diagramgroup32( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 45 );
    retrieveWf( wfs, w_cx, nevt, 84 );
    retrieveWf( wfs, w_cx, nevt, 86 );
    retrieveWf( wfs, w_cx, nevt, 88 );
    retrieveWf( wfs, w_cx, nevt, 100 );
    retrieveWf( wfs, w_cx, nevt, 102 );
    retrieveWf( wfs, w_cx, nevt, 113 );
    retrieveWf( wfs, w_cx, nevt, 116 );
    retrieveWf( wfs, w_cx, nevt, 125 );
    retrieveWf( wfs, w_cx, nevt, 131 );
    retrieveWf( wfs, w_cx, nevt, 151 );
    retrieveWf( wfs, w_cx, nevt, 152 );
    retrieveWf( wfs, w_cx, nevt, 153 );
    retrieveWf( wfs, w_cx, nevt, 156 );
    retrieveWf( wfs, w_cx, nevt, 158 );
    retrieveWf( wfs, w_cx, nevt, 161 );
    retrieveWf( wfs, w_cx, nevt, 163 );
    retrieveWf( wfs, w_cx, nevt, 169 );
    retrieveWf( wfs, w_cx, nevt, 185 );
    retrieveWf( wfs, w_cx, nevt, 187 );
    retrieveWf( wfs, w_cx, nevt, 189 );
    retrieveWf( wfs, w_cx, nevt, 190 );
    retrieveWf( wfs, w_cx, nevt, 191 );
    retrieveWf( wfs, w_cx, nevt, 193 );
    retrieveWf( wfs, w_cx, nevt, 207 );
    retrieveWf( wfs, w_cx, nevt, 209 );
    retrieveWf( wfs, w_cx, nevt, 359 );
    retrieveWf( wfs, w_cx, nevt, 360 );
    retrieveWf( wfs, w_cx, nevt, 361 );
    retrieveWf( wfs, w_cx, nevt, 362 );
    retrieveWf( wfs, w_cx, nevt, 435 );
    retrieveWf( wfs, w_cx, nevt, 437 );
    retrieveWf( wfs, w_cx, nevt, 438 );
    retrieveWf( wfs, w_cx, nevt, 442 );
    retrieveWf( wfs, w_cx, nevt, 444 );
    retrieveWf( wfs, w_cx, nevt, 445 );
    retrieveWf( wfs, w_cx, nevt, 447 );
    retrieveWf( wfs, w_cx, nevt, 448 );
    retrieveWf( wfs, w_cx, nevt, 449 );
    retrieveWf( wfs, w_cx, nevt, 450 );
    retrieveWf( wfs, w_cx, nevt, 451 );
    retrieveWf( wfs, w_cx, nevt, 452 );
    retrieveWf( wfs, w_cx, nevt, 453 );
    retrieveWf( wfs, w_cx, nevt, 454 );
    retrieveWf( wfs, w_cx, nevt, 455 );
    retrieveWf( wfs, w_cx, nevt, 458 );
    retrieveWf( wfs, w_cx, nevt, 459 );
    retrieveWf( wfs, w_cx, nevt, 461 );
    retrieveWf( wfs, w_cx, nevt, 464 );
    retrieveWf( wfs, w_cx, nevt, 477 );
    retrieveWf( wfs, w_cx, nevt, 478 );
    retrieveWf( wfs, w_cx, nevt, 488 );
    retrieveWf( wfs, w_cx, nevt, 497 );
    retrieveWf( wfs, w_cx, nevt, 505 );
    retrieveWf( wfs, w_cx, nevt, 509 );
    retrieveWf( wfs, w_cx, nevt, 510 );
#endif
#endif

    // *** DIAGRAM 3101 OF 15495 ***
    // Wavefunction(s) for diagram number 3101
    // (none)
    // Amplitude(s) for diagram number 3101
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[509], w_fp[163], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[353] -= amp_sv[0];

    // *** DIAGRAM 3102 OF 15495 ***
    // Wavefunction(s) for diagram number 3102
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[156], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[512] );
    // Amplitude(s) for diagram number 3102
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[444], w_fp[512], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[287] -= amp_sv[0];

    // *** DIAGRAM 3103 OF 15495 ***
    // Wavefunction(s) for diagram number 3103
    // (none)
    // Amplitude(s) for diagram number 3103
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[445], w_fp[512], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[281] -= amp_sv[0];

    // *** DIAGRAM 3104 OF 15495 ***
    // Wavefunction(s) for diagram number 3104
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[442], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[513] );
    // Amplitude(s) for diagram number 3104
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[513], w_fp[161], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[333] -= amp_sv[0];

    // *** DIAGRAM 3105 OF 15495 ***
    // Wavefunction(s) for diagram number 3105
    // (none)
    // Amplitude(s) for diagram number 3105
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[445], w_fp[161], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[323] -= amp_sv[0];

    // *** DIAGRAM 3106 OF 15495 ***
    // Wavefunction(s) for diagram number 3106
    // (none)
    // Amplitude(s) for diagram number 3106
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[513], w_fp[163], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[357] -= amp_sv[0];

    // *** DIAGRAM 3107 OF 15495 ***
    // Wavefunction(s) for diagram number 3107
    // (none)
    // Amplitude(s) for diagram number 3107
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[444], w_fp[163], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[347] -= amp_sv[0];

    // *** DIAGRAM 3108 OF 15495 ***
    // Wavefunction(s) for diagram number 3108
    // (none)
    // Amplitude(s) for diagram number 3108
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[448], w_fp[512], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[285] -= amp_sv[0];

    // *** DIAGRAM 3109 OF 15495 ***
    // Wavefunction(s) for diagram number 3109
    // (none)
    // Amplitude(s) for diagram number 3109
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[449], w_fp[512], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[275] -= amp_sv[0];

    // *** DIAGRAM 3110 OF 15495 ***
    // Wavefunction(s) for diagram number 3110
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[447], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[514] );
    // Amplitude(s) for diagram number 3110
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[514], w_fp[158], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[309] -= amp_sv[0];

    // *** DIAGRAM 3111 OF 15495 ***
    // Wavefunction(s) for diagram number 3111
    // (none)
    // Amplitude(s) for diagram number 3111
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[449], w_fp[158], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[299] -= amp_sv[0];

    // *** DIAGRAM 3112 OF 15495 ***
    // Wavefunction(s) for diagram number 3112
    // (none)
    // Amplitude(s) for diagram number 3112
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[514], w_fp[163], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[351] -= amp_sv[0];

    // *** DIAGRAM 3113 OF 15495 ***
    // Wavefunction(s) for diagram number 3113
    // (none)
    // Amplitude(s) for diagram number 3113
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[448], w_fp[163], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[345] -= amp_sv[0];

    // *** DIAGRAM 3114 OF 15495 ***
    // Wavefunction(s) for diagram number 3114
    // (none)
    // Amplitude(s) for diagram number 3114
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[452], w_fp[512], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[279] -= amp_sv[0];

    // *** DIAGRAM 3115 OF 15495 ***
    // Wavefunction(s) for diagram number 3115
    // (none)
    // Amplitude(s) for diagram number 3115
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[453], w_fp[512], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[273] -= amp_sv[0];

    // *** DIAGRAM 3116 OF 15495 ***
    // Wavefunction(s) for diagram number 3116
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[451], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[515] );
    // Amplitude(s) for diagram number 3116
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[515], w_fp[158], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[303] -= amp_sv[0];

    // *** DIAGRAM 3117 OF 15495 ***
    // Wavefunction(s) for diagram number 3117
    // (none)
    // Amplitude(s) for diagram number 3117
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[453], w_fp[158], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[297] -= amp_sv[0];

    // *** DIAGRAM 3118 OF 15495 ***
    // Wavefunction(s) for diagram number 3118
    // (none)
    // Amplitude(s) for diagram number 3118
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[515], w_fp[161], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[327] -= amp_sv[0];

    // *** DIAGRAM 3119 OF 15495 ***
    // Wavefunction(s) for diagram number 3119
    // (none)
    // Amplitude(s) for diagram number 3119
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[452], w_fp[161], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[321] -= amp_sv[0];

    // *** DIAGRAM 3120 OF 15495 ***
    // Wavefunction(s) for diagram number 3120
    // (none)
    // Amplitude(s) for diagram number 3120
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[185], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[305] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3121 OF 15495 ***
    // Wavefunction(s) for diagram number 3121
    // (none)
    // Amplitude(s) for diagram number 3121
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[163], w_fp[113], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[353] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3122 OF 15495 ***
    // Wavefunction(s) for diagram number 3122
    // (none)
    // Amplitude(s) for diagram number 3122
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[156], w_fp[116], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[305] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];

    // *** DIAGRAM 3123 OF 15495 ***
    // Wavefunction(s) for diagram number 3123
    // (none)
    // Amplitude(s) for diagram number 3123
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[497], w_fp[359], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[345] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3124 OF 15495 ***
    // Wavefunction(s) for diagram number 3124
    // (none)
    // Amplitude(s) for diagram number 3124
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[497], w_fp[1], w_fp[116], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[305] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3125 OF 15495 ***
    // Wavefunction(s) for diagram number 3125
    // (none)
    // Amplitude(s) for diagram number 3125
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[113], w_fp[7], w_fp[497], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[345] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[113], w_fp[7], w_fp[497], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[285] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[305] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[345] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[113], w_fp[7], w_fp[497], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[273] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[305] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3126 OF 15495 ***
    // Wavefunction(s) for diagram number 3126
    // (none)
    // Amplitude(s) for diagram number 3126
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[464], w_fp[512], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[285] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3127 OF 15495 ***
    // Wavefunction(s) for diagram number 3127
    // (none)
    // Amplitude(s) for diagram number 3127
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[464], w_fp[163], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[345] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3128 OF 15495 ***
    // Wavefunction(s) for diagram number 3128
    // (none)
    // Amplitude(s) for diagram number 3128
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[512], w_fp[113], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[273] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3129 OF 15495 ***
    // Wavefunction(s) for diagram number 3129
    // (none)
    // Amplitude(s) for diagram number 3129
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[156], w_fp[359], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[273] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[327] += amp_sv[0];

    // *** DIAGRAM 3130 OF 15495 ***
    // Wavefunction(s) for diagram number 3130
    // (none)
    // Amplitude(s) for diagram number 3130
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[185], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[303] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3131 OF 15495 ***
    // Wavefunction(s) for diagram number 3131
    // (none)
    // Amplitude(s) for diagram number 3131
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[45], w_fp[512], w_fp[116], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[273] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];

    // *** DIAGRAM 3132 OF 15495 ***
    // Wavefunction(s) for diagram number 3132
    // (none)
    // Amplitude(s) for diagram number 3132
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[45], w_fp[163], w_fp[359], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[345] += amp_sv[0];
    jamp_sv[347] -= amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];

    // *** DIAGRAM 3133 OF 15495 ***
    // Wavefunction(s) for diagram number 3133
    // (none)
    // Amplitude(s) for diagram number 3133
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[187], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[311] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3134 OF 15495 ***
    // Wavefunction(s) for diagram number 3134
    // (none)
    // Amplitude(s) for diagram number 3134
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[161], w_fp[100], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[329] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3135 OF 15495 ***
    // Wavefunction(s) for diagram number 3135
    // (none)
    // Amplitude(s) for diagram number 3135
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[156], w_fp[125], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[311] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];

    // *** DIAGRAM 3136 OF 15495 ***
    // Wavefunction(s) for diagram number 3136
    // (none)
    // Amplitude(s) for diagram number 3136
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[497], w_fp[360], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[275] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[309] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[351] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3137 OF 15495 ***
    // Wavefunction(s) for diagram number 3137
    // (none)
    // Amplitude(s) for diagram number 3137
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[497], w_fp[1], w_fp[125], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[275] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[311] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3138 OF 15495 ***
    // Wavefunction(s) for diagram number 3138
    // (none)
    // Amplitude(s) for diagram number 3138
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[100], w_fp[6], w_fp[497], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[275] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[309] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[351] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[100], w_fp[6], w_fp[497], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[309] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[311] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[351] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[100], w_fp[6], w_fp[497], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[275] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[311] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3139 OF 15495 ***
    // Wavefunction(s) for diagram number 3139
    // (none)
    // Amplitude(s) for diagram number 3139
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[461], w_fp[512], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[279] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3140 OF 15495 ***
    // Wavefunction(s) for diagram number 3140
    // (none)
    // Amplitude(s) for diagram number 3140
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[461], w_fp[161], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3141 OF 15495 ***
    // Wavefunction(s) for diagram number 3141
    // (none)
    // Amplitude(s) for diagram number 3141
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[447], w_fp[512], w_fp[100], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[275] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3142 OF 15495 ***
    // Wavefunction(s) for diagram number 3142
    // (none)
    // Amplitude(s) for diagram number 3142
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[447], w_fp[156], w_fp[360], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[275] += amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[351] += amp_sv[0];

    // *** DIAGRAM 3143 OF 15495 ***
    // Wavefunction(s) for diagram number 3143
    // (none)
    // Amplitude(s) for diagram number 3143
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[447], w_fp[187], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[309] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[351] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3144 OF 15495 ***
    // Wavefunction(s) for diagram number 3144
    // (none)
    // Amplitude(s) for diagram number 3144
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[45], w_fp[512], w_fp[125], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[275] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[281] += amp_sv[0];
    jamp_sv[285] -= amp_sv[0];

    // *** DIAGRAM 3145 OF 15495 ***
    // Wavefunction(s) for diagram number 3145
    // (none)
    // Amplitude(s) for diagram number 3145
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[45], w_fp[161], w_fp[360], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[321] += amp_sv[0];
    jamp_sv[323] -= amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];

    // *** DIAGRAM 3146 OF 15495 ***
    // Wavefunction(s) for diagram number 3146
    // (none)
    // Amplitude(s) for diagram number 3146
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[158], w_fp[84], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[305] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[311] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3147 OF 15495 ***
    // Wavefunction(s) for diagram number 3147
    // (none)
    // Amplitude(s) for diagram number 3147
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[189], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[335] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3148 OF 15495 ***
    // Wavefunction(s) for diagram number 3148
    // (none)
    // Amplitude(s) for diagram number 3148
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[156], w_fp[131], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[305] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];

    // *** DIAGRAM 3149 OF 15495 ***
    // Wavefunction(s) for diagram number 3149
    // (none)
    // Amplitude(s) for diagram number 3149
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[497], w_fp[361], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[281] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[299] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[305] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[311] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3150 OF 15495 ***
    // Wavefunction(s) for diagram number 3150
    // (none)
    // Amplitude(s) for diagram number 3150
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[497], w_fp[1], w_fp[131], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[305] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[311] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3151 OF 15495 ***
    // Wavefunction(s) for diagram number 3151
    // (none)
    // Amplitude(s) for diagram number 3151
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[5], w_fp[84], w_fp[497], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[299] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[5], w_fp[84], w_fp[497], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[281] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[299] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[305] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[311] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[5], w_fp[84], w_fp[497], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[273] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[305] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[311] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3152 OF 15495 ***
    // Wavefunction(s) for diagram number 3152
    // (none)
    // Amplitude(s) for diagram number 3152
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[442], w_fp[512], w_fp[84], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3153 OF 15495 ***
    // Wavefunction(s) for diagram number 3153
    // (none)
    // Amplitude(s) for diagram number 3153
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[442], w_fp[156], w_fp[361], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[281] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[357] += amp_sv[0];

    // *** DIAGRAM 3154 OF 15495 ***
    // Wavefunction(s) for diagram number 3154
    // (none)
    // Amplitude(s) for diagram number 3154
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[442], w_fp[189], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[333] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3155 OF 15495 ***
    // Wavefunction(s) for diagram number 3155
    // (none)
    // Amplitude(s) for diagram number 3155
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[458], w_fp[512], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[273] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3156 OF 15495 ***
    // Wavefunction(s) for diagram number 3156
    // (none)
    // Amplitude(s) for diagram number 3156
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[458], w_fp[158], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[297] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[299] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3157 OF 15495 ***
    // Wavefunction(s) for diagram number 3157
    // (none)
    // Amplitude(s) for diagram number 3157
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[45], w_fp[512], w_fp[131], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[273] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];

    // *** DIAGRAM 3158 OF 15495 ***
    // Wavefunction(s) for diagram number 3158
    // (none)
    // Amplitude(s) for diagram number 3158
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[45], w_fp[158], w_fp[361], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[297] += amp_sv[0];
    jamp_sv[299] -= amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];

    // *** DIAGRAM 3159 OF 15495 ***
    // Wavefunction(s) for diagram number 3159
    // (none)
    // Amplitude(s) for diagram number 3159
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[156], w_fp[151], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[305] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[156], w_fp[152], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[353] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[156], w_fp[153], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[353] += amp_sv[0];
    jamp_sv[359] -= amp_sv[0];

    // *** DIAGRAM 3160 OF 15495 ***
    // Wavefunction(s) for diagram number 3160
    // (none)
    // Amplitude(s) for diagram number 3160
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[151], w_fp[497], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[305] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[311] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[152], w_fp[497], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[275] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[311] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[153], w_fp[497], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[273] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[305] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3161 OF 15495 ***
    // Wavefunction(s) for diagram number 3161
    // (none)
    // Amplitude(s) for diagram number 3161
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[455], w_fp[156], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[273] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[477], w_fp[156], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[478], w_fp[156], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];

    // *** DIAGRAM 3162 OF 15495 ***
    // Wavefunction(s) for diagram number 3162
    // (none)
    // Amplitude(s) for diagram number 3162
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[509], w_fp[190], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[431] -= amp_sv[0];

    // *** DIAGRAM 3163 OF 15495 ***
    // Wavefunction(s) for diagram number 3163
    // (none)
    // Amplitude(s) for diagram number 3163
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[510], w_fp[190], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[425] -= amp_sv[0];

    // *** DIAGRAM 3164 OF 15495 ***
    // Wavefunction(s) for diagram number 3164
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[505], w_fp[4], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[478] );
    // Amplitude(s) for diagram number 3164
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[478], w_fp[191], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[455] -= amp_sv[0];

    // *** DIAGRAM 3165 OF 15495 ***
    // Wavefunction(s) for diagram number 3165
    // (none)
    // Amplitude(s) for diagram number 3165
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[510], w_fp[191], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[449] -= amp_sv[0];

    // *** DIAGRAM 3166 OF 15495 ***
    // Wavefunction(s) for diagram number 3166
    // (none)
    // Amplitude(s) for diagram number 3166
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[478], w_fp[193], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[479] -= amp_sv[0];

    // *** DIAGRAM 3167 OF 15495 ***
    // Wavefunction(s) for diagram number 3167
    // (none)
    // Amplitude(s) for diagram number 3167
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[509], w_fp[193], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[473] -= amp_sv[0];

    // *** DIAGRAM 3168 OF 15495 ***
    // Wavefunction(s) for diagram number 3168
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[169], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[477] );
    // Amplitude(s) for diagram number 3168
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[437], w_fp[477], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[407] -= amp_sv[0];

    // *** DIAGRAM 3169 OF 15495 ***
    // Wavefunction(s) for diagram number 3169
    // (none)
    // Amplitude(s) for diagram number 3169
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[438], w_fp[477], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[401] -= amp_sv[0];

    // *** DIAGRAM 3170 OF 15495 ***
    // Wavefunction(s) for diagram number 3170
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[435], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[455] );
    // Amplitude(s) for diagram number 3170
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[455], w_fp[191], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[453] -= amp_sv[0];

    // *** DIAGRAM 3171 OF 15495 ***
    // Wavefunction(s) for diagram number 3171
    // (none)
    // Amplitude(s) for diagram number 3171
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[438], w_fp[191], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[443] -= amp_sv[0];

    // *** DIAGRAM 3172 OF 15495 ***
    // Wavefunction(s) for diagram number 3172
    // (none)
    // Amplitude(s) for diagram number 3172
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[455], w_fp[193], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[477] -= amp_sv[0];

    // *** DIAGRAM 3173 OF 15495 ***
    // Wavefunction(s) for diagram number 3173
    // (none)
    // Amplitude(s) for diagram number 3173
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[437], w_fp[193], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[467] -= amp_sv[0];

    // *** DIAGRAM 3174 OF 15495 ***
    // Wavefunction(s) for diagram number 3174
    // (none)
    // Amplitude(s) for diagram number 3174
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[450], w_fp[477], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[405] -= amp_sv[0];

    // *** DIAGRAM 3175 OF 15495 ***
    // Wavefunction(s) for diagram number 3175
    // (none)
    // Amplitude(s) for diagram number 3175
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[449], w_fp[477], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[395] -= amp_sv[0];

    // *** DIAGRAM 3176 OF 15495 ***
    // Wavefunction(s) for diagram number 3176
    // (none)
    // Amplitude(s) for diagram number 3176
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[514], w_fp[190], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[429] -= amp_sv[0];

    // *** DIAGRAM 3177 OF 15495 ***
    // Wavefunction(s) for diagram number 3177
    // (none)
    // Amplitude(s) for diagram number 3177
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[449], w_fp[190], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[419] -= amp_sv[0];

    // *** DIAGRAM 3178 OF 15495 ***
    // Wavefunction(s) for diagram number 3178
    // (none)
    // Amplitude(s) for diagram number 3178
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[514], w_fp[193], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[471] -= amp_sv[0];

    // *** DIAGRAM 3179 OF 15495 ***
    // Wavefunction(s) for diagram number 3179
    // (none)
    // Amplitude(s) for diagram number 3179
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[450], w_fp[193], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[465] -= amp_sv[0];

    // *** DIAGRAM 3180 OF 15495 ***
    // Wavefunction(s) for diagram number 3180
    // (none)
    // Amplitude(s) for diagram number 3180
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[454], w_fp[477], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[399] -= amp_sv[0];

    // *** DIAGRAM 3181 OF 15495 ***
    // Wavefunction(s) for diagram number 3181
    // (none)
    // Amplitude(s) for diagram number 3181
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[453], w_fp[477], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[393] -= amp_sv[0];

    // *** DIAGRAM 3182 OF 15495 ***
    // Wavefunction(s) for diagram number 3182
    // (none)
    // Amplitude(s) for diagram number 3182
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[515], w_fp[190], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[423] -= amp_sv[0];

    // *** DIAGRAM 3183 OF 15495 ***
    // Wavefunction(s) for diagram number 3183
    // (none)
    // Amplitude(s) for diagram number 3183
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[453], w_fp[190], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[417] -= amp_sv[0];

    // *** DIAGRAM 3184 OF 15495 ***
    // Wavefunction(s) for diagram number 3184
    // (none)
    // Amplitude(s) for diagram number 3184
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[515], w_fp[191], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[447] -= amp_sv[0];

    // *** DIAGRAM 3185 OF 15495 ***
    // Wavefunction(s) for diagram number 3185
    // (none)
    // Amplitude(s) for diagram number 3185
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[454], w_fp[191], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[441] -= amp_sv[0];

    // *** DIAGRAM 3186 OF 15495 ***
    // Wavefunction(s) for diagram number 3186
    // (none)
    // Amplitude(s) for diagram number 3186
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[207], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[425] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3187 OF 15495 ***
    // Wavefunction(s) for diagram number 3187
    // (none)
    // Amplitude(s) for diagram number 3187
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[193], w_fp[86], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[473] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3188 OF 15495 ***
    // Wavefunction(s) for diagram number 3188
    // (none)
    // Amplitude(s) for diagram number 3188
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[169], w_fp[88], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[425] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];

    // *** DIAGRAM 3189 OF 15495 ***
    // Wavefunction(s) for diagram number 3189
    // (none)
    // Amplitude(s) for diagram number 3189
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[362], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[467] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3190 OF 15495 ***
    // Wavefunction(s) for diagram number 3190
    // (none)
    // Amplitude(s) for diagram number 3190
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[1], w_fp[88], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3191 OF 15495 ***
    // Wavefunction(s) for diagram number 3191
    // (none)
    // Amplitude(s) for diagram number 3191
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[86], w_fp[7], w_fp[488], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[467] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[86], w_fp[7], w_fp[488], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[405] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[467] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[86], w_fp[7], w_fp[488], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[393] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3192 OF 15495 ***
    // Wavefunction(s) for diagram number 3192
    // (none)
    // Amplitude(s) for diagram number 3192
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[459], w_fp[477], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[405] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3193 OF 15495 ***
    // Wavefunction(s) for diagram number 3193
    // (none)
    // Amplitude(s) for diagram number 3193
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[459], w_fp[193], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[465] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[467] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3194 OF 15495 ***
    // Wavefunction(s) for diagram number 3194
    // (none)
    // Amplitude(s) for diagram number 3194
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[477], w_fp[86], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[393] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3195 OF 15495 ***
    // Wavefunction(s) for diagram number 3195
    // (none)
    // Amplitude(s) for diagram number 3195
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[169], w_fp[362], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[393] += amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[447] += amp_sv[0];

    // *** DIAGRAM 3196 OF 15495 ***
    // Wavefunction(s) for diagram number 3196
    // (none)
    // Amplitude(s) for diagram number 3196
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[207], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[423] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3197 OF 15495 ***
    // Wavefunction(s) for diagram number 3197
    // (none)
    // Amplitude(s) for diagram number 3197
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[45], w_fp[477], w_fp[88], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[393] += amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];

    // *** DIAGRAM 3198 OF 15495 ***
    // Wavefunction(s) for diagram number 3198
    // (none)
    // Amplitude(s) for diagram number 3198
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[45], w_fp[193], w_fp[362], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[465] += amp_sv[0];
    jamp_sv[467] -= amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];

    // *** DIAGRAM 3199 OF 15495 ***
    // Wavefunction(s) for diagram number 3199
    // (none)
    // Amplitude(s) for diagram number 3199
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[209], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[431] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3200 OF 15495 ***
    // Wavefunction(s) for diagram number 3200
    // (none)
    // Amplitude(s) for diagram number 3200
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[191], w_fp[102], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[449] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] -= cxtype( 0, 1 ) * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 455 );
    storeWf( wfs, w_cx, nevt, 477 );
    storeWf( wfs, w_cx, nevt, 478 );
    storeWf( wfs, w_cx, nevt, 512 );
    storeWf( wfs, w_cx, nevt, 513 );
    storeWf( wfs, w_cx, nevt, 514 );
    storeWf( wfs, w_cx, nevt, 515 );
#endif
#endif
  }

  //--------------------------------------------------------------------------
}
