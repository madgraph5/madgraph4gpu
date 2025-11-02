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
  diagramgroup88( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 84 );
    retrieveWf( wfs, w_cx, nevt, 86 );
    retrieveWf( wfs, w_cx, nevt, 88 );
    retrieveWf( wfs, w_cx, nevt, 99 );
    retrieveWf( wfs, w_cx, nevt, 100 );
    retrieveWf( wfs, w_cx, nevt, 102 );
    retrieveWf( wfs, w_cx, nevt, 104 );
    retrieveWf( wfs, w_cx, nevt, 112 );
    retrieveWf( wfs, w_cx, nevt, 113 );
    retrieveWf( wfs, w_cx, nevt, 116 );
    retrieveWf( wfs, w_cx, nevt, 123 );
    retrieveWf( wfs, w_cx, nevt, 125 );
    retrieveWf( wfs, w_cx, nevt, 129 );
    retrieveWf( wfs, w_cx, nevt, 130 );
    retrieveWf( wfs, w_cx, nevt, 131 );
    retrieveWf( wfs, w_cx, nevt, 151 );
    retrieveWf( wfs, w_cx, nevt, 152 );
    retrieveWf( wfs, w_cx, nevt, 153 );
    retrieveWf( wfs, w_cx, nevt, 168 );
    retrieveWf( wfs, w_cx, nevt, 170 );
    retrieveWf( wfs, w_cx, nevt, 171 );
    retrieveWf( wfs, w_cx, nevt, 196 );
    retrieveWf( wfs, w_cx, nevt, 198 );
    retrieveWf( wfs, w_cx, nevt, 199 );
    retrieveWf( wfs, w_cx, nevt, 200 );
    retrieveWf( wfs, w_cx, nevt, 216 );
    retrieveWf( wfs, w_cx, nevt, 218 );
    retrieveWf( wfs, w_cx, nevt, 243 );
    retrieveWf( wfs, w_cx, nevt, 259 );
    retrieveWf( wfs, w_cx, nevt, 436 );
    retrieveWf( wfs, w_cx, nevt, 439 );
    retrieveWf( wfs, w_cx, nevt, 441 );
    retrieveWf( wfs, w_cx, nevt, 443 );
    retrieveWf( wfs, w_cx, nevt, 449 );
    retrieveWf( wfs, w_cx, nevt, 453 );
    retrieveWf( wfs, w_cx, nevt, 460 );
    retrieveWf( wfs, w_cx, nevt, 463 );
    retrieveWf( wfs, w_cx, nevt, 465 );
    retrieveWf( wfs, w_cx, nevt, 466 );
    retrieveWf( wfs, w_cx, nevt, 467 );
    retrieveWf( wfs, w_cx, nevt, 468 );
    retrieveWf( wfs, w_cx, nevt, 525 );
    retrieveWf( wfs, w_cx, nevt, 536 );
    retrieveWf( wfs, w_cx, nevt, 541 );
    retrieveWf( wfs, w_cx, nevt, 549 );
    retrieveWf( wfs, w_cx, nevt, 585 );
    retrieveWf( wfs, w_cx, nevt, 588 );
#endif
#endif

    // *** DIAGRAM 8701 OF 15495 ***
    // Wavefunction(s) for diagram number 8701
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[441], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[505] );
    // Amplitude(s) for diagram number 8701
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[216], w_fp[505], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[221] -= amp_sv[0];

    // *** DIAGRAM 8702 OF 15495 ***
    // Wavefunction(s) for diagram number 8702
    // (none)
    // Amplitude(s) for diagram number 8702
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[588], w_fp[441], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[235] -= amp_sv[0];

    // *** DIAGRAM 8703 OF 15495 ***
    // Wavefunction(s) for diagram number 8703
    // (none)
    // Amplitude(s) for diagram number 8703
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[198], w_fp[505], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[219] -= amp_sv[0];

    // *** DIAGRAM 8704 OF 15495 ***
    // Wavefunction(s) for diagram number 8704
    // (none)
    // Amplitude(s) for diagram number 8704
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[541], w_fp[441], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[229] -= amp_sv[0];

    // *** DIAGRAM 8705 OF 15495 ***
    // Wavefunction(s) for diagram number 8705
    // (none)
    // Amplitude(s) for diagram number 8705
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[243], w_fp[449], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[141] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8706 OF 15495 ***
    // Wavefunction(s) for diagram number 8706
    // (none)
    // Amplitude(s) for diagram number 8706
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[199], w_fp[449], w_fp[113], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8707 OF 15495 ***
    // Wavefunction(s) for diagram number 8707
    // (none)
    // Amplitude(s) for diagram number 8707
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[449], w_fp[116], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];

    // *** DIAGRAM 8708 OF 15495 ***
    // Wavefunction(s) for diagram number 8708
    // (none)
    // Amplitude(s) for diagram number 8708
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[536], w_fp[465], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[184] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8709 OF 15495 ***
    // Wavefunction(s) for diagram number 8709
    // (none)
    // Amplitude(s) for diagram number 8709
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[536], w_fp[441], w_fp[113], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[232] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8710 OF 15495 ***
    // Wavefunction(s) for diagram number 8710
    // (none)
    // Amplitude(s) for diagram number 8710
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[536], w_fp[259], w_fp[116], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[184] += amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];

    // *** DIAGRAM 8711 OF 15495 ***
    // Wavefunction(s) for diagram number 8711
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[113], COUPs[0], 1.0, 0., 0., w_fp[559] );
    // Amplitude(s) for diagram number 8711
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[559], w_fp[549], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[181] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[221] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8712 OF 15495 ***
    // Wavefunction(s) for diagram number 8712
    // (none)
    // Amplitude(s) for diagram number 8712
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[441], w_fp[559], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[219] += amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];

    // *** DIAGRAM 8713 OF 15495 ***
    // Wavefunction(s) for diagram number 8713
    // (none)
    // Amplitude(s) for diagram number 8713
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[199], w_fp[259], w_fp[559], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[181] -= amp_sv[0];
    jamp_sv[205] += amp_sv[0];

    // *** DIAGRAM 8714 OF 15495 ***
    // Wavefunction(s) for diagram number 8714
    // (none)
    // Amplitude(s) for diagram number 8714
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[549], w_fp[116], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[184] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8715 OF 15495 ***
    // Wavefunction(s) for diagram number 8715
    // (none)
    // Amplitude(s) for diagram number 8715
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[199], w_fp[465], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[181] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8716 OF 15495 ***
    // Wavefunction(s) for diagram number 8716
    // (none)
    // Amplitude(s) for diagram number 8716
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[243], w_fp[441], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[221] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8717 OF 15495 ***
    // Wavefunction(s) for diagram number 8717
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[113], w_fp[7], COUPs[2], 1.0, 0., 0., w_fp[479] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[113], w_fp[7], COUPs[2], 1.0, 0., 0., w_fp[474] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[113], w_fp[7], COUPs[2], 1.0, 0., 0., w_fp[444] );
    // Amplitude(s) for diagram number 8717
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[259], w_fp[479], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[184] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[259], w_fp[474], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[141] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[181] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[184] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[221] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[259], w_fp[444], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[181] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[221] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8718 OF 15495 ***
    // Wavefunction(s) for diagram number 8718
    // (none)
    // Amplitude(s) for diagram number 8718
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[123], w_fp[449], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[135] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8719 OF 15495 ***
    // Wavefunction(s) for diagram number 8719
    // (none)
    // Amplitude(s) for diagram number 8719
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[198], w_fp[449], w_fp[100], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8720 OF 15495 ***
    // Wavefunction(s) for diagram number 8720
    // (none)
    // Amplitude(s) for diagram number 8720
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[449], w_fp[125], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[131] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[141] -= amp_sv[0];

    // *** DIAGRAM 8721 OF 15495 ***
    // Wavefunction(s) for diagram number 8721
    // (none)
    // Amplitude(s) for diagram number 8721
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[536], w_fp[466], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[190] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8722 OF 15495 ***
    // Wavefunction(s) for diagram number 8722
    // (none)
    // Amplitude(s) for diagram number 8722
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[536], w_fp[439], w_fp[100], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[208] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8723 OF 15495 ***
    // Wavefunction(s) for diagram number 8723
    // (none)
    // Amplitude(s) for diagram number 8723
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[536], w_fp[259], w_fp[125], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[190] += amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[232] -= amp_sv[0];

    // *** DIAGRAM 8724 OF 15495 ***
    // Wavefunction(s) for diagram number 8724
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[100], COUPs[0], 1.0, 0., 0., w_fp[522] );
    // Amplitude(s) for diagram number 8724
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[522], w_fp[549], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[187] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8725 OF 15495 ***
    // Wavefunction(s) for diagram number 8725
    // (none)
    // Amplitude(s) for diagram number 8725
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[439], w_fp[522], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[195] += amp_sv[0];
    jamp_sv[197] -= amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[214] += amp_sv[0];

    // *** DIAGRAM 8726 OF 15495 ***
    // Wavefunction(s) for diagram number 8726
    // (none)
    // Amplitude(s) for diagram number 8726
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[198], w_fp[259], w_fp[522], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[131] += amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[229] += amp_sv[0];

    // *** DIAGRAM 8727 OF 15495 ***
    // Wavefunction(s) for diagram number 8727
    // (none)
    // Amplitude(s) for diagram number 8727
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[549], w_fp[125], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8728 OF 15495 ***
    // Wavefunction(s) for diagram number 8728
    // (none)
    // Amplitude(s) for diagram number 8728
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[198], w_fp[466], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[187] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8729 OF 15495 ***
    // Wavefunction(s) for diagram number 8729
    // (none)
    // Amplitude(s) for diagram number 8729
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[123], w_fp[439], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8730 OF 15495 ***
    // Wavefunction(s) for diagram number 8730
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[100], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[523] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[100], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[580] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[100], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[528] );
    // Amplitude(s) for diagram number 8730
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[259], w_fp[523], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[259], w_fp[580], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[135] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[187] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[259], w_fp[528], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[187] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8731 OF 15495 ***
    // Wavefunction(s) for diagram number 8731
    // (none)
    // Amplitude(s) for diagram number 8731
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[216], w_fp[449], w_fp[84], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[137] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8732 OF 15495 ***
    // Wavefunction(s) for diagram number 8732
    // (none)
    // Amplitude(s) for diagram number 8732
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[129], w_fp[449], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8733 OF 15495 ***
    // Wavefunction(s) for diagram number 8733
    // (none)
    // Amplitude(s) for diagram number 8733
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[449], w_fp[131], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];

    // *** DIAGRAM 8734 OF 15495 ***
    // Wavefunction(s) for diagram number 8734
    // (none)
    // Amplitude(s) for diagram number 8734
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[536], w_fp[436], w_fp[84], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[184] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8735 OF 15495 ***
    // Wavefunction(s) for diagram number 8735
    // (none)
    // Amplitude(s) for diagram number 8735
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[536], w_fp[467], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8736 OF 15495 ***
    // Wavefunction(s) for diagram number 8736
    // (none)
    // Amplitude(s) for diagram number 8736
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[536], w_fp[259], w_fp[131], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[184] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];

    // *** DIAGRAM 8737 OF 15495 ***
    // Wavefunction(s) for diagram number 8737
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[84], COUPs[0], 1.0, 0., 0., w_fp[516] );
    // Amplitude(s) for diagram number 8737
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[516], w_fp[549], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[137] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[184] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8738 OF 15495 ***
    // Wavefunction(s) for diagram number 8738
    // (none)
    // Amplitude(s) for diagram number 8738
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[436], w_fp[516], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[171] += amp_sv[0];
    jamp_sv[173] -= amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];

    // *** DIAGRAM 8739 OF 15495 ***
    // Wavefunction(s) for diagram number 8739
    // (none)
    // Amplitude(s) for diagram number 8739
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[216], w_fp[259], w_fp[516], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[137] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[235] += amp_sv[0];

    // *** DIAGRAM 8740 OF 15495 ***
    // Wavefunction(s) for diagram number 8740
    // (none)
    // Amplitude(s) for diagram number 8740
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[549], w_fp[131], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[184] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8741 OF 15495 ***
    // Wavefunction(s) for diagram number 8741
    // (none)
    // Amplitude(s) for diagram number 8741
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[129], w_fp[436], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8742 OF 15495 ***
    // Wavefunction(s) for diagram number 8742
    // (none)
    // Amplitude(s) for diagram number 8742
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[216], w_fp[467], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8743 OF 15495 ***
    // Wavefunction(s) for diagram number 8743
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[5], w_fp[84], COUPs[2], 1.0, 0., 0., w_fp[549] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[5], w_fp[84], COUPs[2], 1.0, 0., 0., w_fp[451] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[5], w_fp[84], COUPs[2], 1.0, 0., 0., w_fp[571] );
    // Amplitude(s) for diagram number 8743
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[259], w_fp[549], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[184] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[259], w_fp[451], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[137] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[184] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[259], w_fp[571], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8744 OF 15495 ***
    // Wavefunction(s) for diagram number 8744
    // (none)
    // Amplitude(s) for diagram number 8744
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[449], w_fp[151], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[449], w_fp[152], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[449], w_fp[153], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];

    // *** DIAGRAM 8745 OF 15495 ***
    // Wavefunction(s) for diagram number 8745
    // (none)
    // Amplitude(s) for diagram number 8745
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[536], w_fp[259], w_fp[151], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[184] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[536], w_fp[259], w_fp[152], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[536], w_fp[259], w_fp[153], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];

    // *** DIAGRAM 8746 OF 15495 ***
    // Wavefunction(s) for diagram number 8746
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[151], COUPs[0], 1.0, 0., 0., w_fp[452] );
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[152], COUPs[0], 1.0, 0., 0., w_fp[488] );
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[153], COUPs[0], 1.0, 0., 0., w_fp[437] );
    // Amplitude(s) for diagram number 8746
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[259], w_fp[452], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[184] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[259], w_fp[488], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[259], w_fp[437], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[184] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8747 OF 15495 ***
    // Wavefunction(s) for diagram number 8747
    // (none)
    // Amplitude(s) for diagram number 8747
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[218], w_fp[453], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[136] -= amp_sv[0];

    // *** DIAGRAM 8748 OF 15495 ***
    // Wavefunction(s) for diagram number 8748
    // (none)
    // Amplitude(s) for diagram number 8748
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[218], w_fp[585], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[142] -= amp_sv[0];

    // *** DIAGRAM 8749 OF 15495 ***
    // Wavefunction(s) for diagram number 8749
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[449], w_fp[4], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[546] );
    // Amplitude(s) for diagram number 8749
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[170], w_fp[546], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[125] -= amp_sv[0];

    // *** DIAGRAM 8750 OF 15495 ***
    // Wavefunction(s) for diagram number 8750
    // (none)
    // Amplitude(s) for diagram number 8750
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[170], w_fp[585], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[139] -= amp_sv[0];

    // *** DIAGRAM 8751 OF 15495 ***
    // Wavefunction(s) for diagram number 8751
    // (none)
    // Amplitude(s) for diagram number 8751
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[171], w_fp[546], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] -= amp_sv[0];

    // *** DIAGRAM 8752 OF 15495 ***
    // Wavefunction(s) for diagram number 8752
    // (none)
    // Amplitude(s) for diagram number 8752
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[171], w_fp[453], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[133] -= amp_sv[0];

    // *** DIAGRAM 8753 OF 15495 ***
    // Wavefunction(s) for diagram number 8753
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[168], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[527] );
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[527], w_fp[6], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[538] );
    // Amplitude(s) for diagram number 8753
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[538], w_fp[443], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[166] -= amp_sv[0];

    // *** DIAGRAM 8754 OF 15495 ***
    // Wavefunction(s) for diagram number 8754
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[527], w_fp[7], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[543] );
    // Amplitude(s) for diagram number 8754
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[543], w_fp[443], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] -= amp_sv[0];

    // *** DIAGRAM 8755 OF 15495 ***
    // Wavefunction(s) for diagram number 8755
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[527], w_fp[4], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[509] );
    // Amplitude(s) for diagram number 8755
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[509], w_fp[439], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[212] -= amp_sv[0];

    // *** DIAGRAM 8756 OF 15495 ***
    // Wavefunction(s) for diagram number 8756
    // (none)
    // Amplitude(s) for diagram number 8756
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[543], w_fp[439], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[202] -= amp_sv[0];

    // *** DIAGRAM 8757 OF 15495 ***
    // Wavefunction(s) for diagram number 8757
    // (none)
    // Amplitude(s) for diagram number 8757
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[509], w_fp[441], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[236] -= amp_sv[0];

    // *** DIAGRAM 8758 OF 15495 ***
    // Wavefunction(s) for diagram number 8758
    // (none)
    // Amplitude(s) for diagram number 8758
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[538], w_fp[441], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[226] -= amp_sv[0];

    // *** DIAGRAM 8759 OF 15495 ***
    // Wavefunction(s) for diagram number 8759
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[443], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[520] );
    // Amplitude(s) for diagram number 8759
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[170], w_fp[520], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[149] -= amp_sv[0];

    // *** DIAGRAM 8760 OF 15495 ***
    // Wavefunction(s) for diagram number 8760
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[170], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[535] );
    // Amplitude(s) for diagram number 8760
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[535], w_fp[443], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[163] -= amp_sv[0];

    // *** DIAGRAM 8761 OF 15495 ***
    // Wavefunction(s) for diagram number 8761
    // (none)
    // Amplitude(s) for diagram number 8761
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[171], w_fp[520], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[147] -= amp_sv[0];

    // *** DIAGRAM 8762 OF 15495 ***
    // Wavefunction(s) for diagram number 8762
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[171], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[450] );
    // Amplitude(s) for diagram number 8762
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[450], w_fp[443], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[157] -= amp_sv[0];

    // *** DIAGRAM 8763 OF 15495 ***
    // Wavefunction(s) for diagram number 8763
    // (none)
    // Amplitude(s) for diagram number 8763
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[218], w_fp[525], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[196] -= amp_sv[0];

    // *** DIAGRAM 8764 OF 15495 ***
    // Wavefunction(s) for diagram number 8764
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[218], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[481] );
    // Amplitude(s) for diagram number 8764
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[481], w_fp[439], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[210] -= amp_sv[0];

    // *** DIAGRAM 8765 OF 15495 ***
    // Wavefunction(s) for diagram number 8765
    // (none)
    // Amplitude(s) for diagram number 8765
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[171], w_fp[525], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[193] -= amp_sv[0];

    // *** DIAGRAM 8766 OF 15495 ***
    // Wavefunction(s) for diagram number 8766
    // (none)
    // Amplitude(s) for diagram number 8766
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[450], w_fp[439], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[199] -= amp_sv[0];

    // *** DIAGRAM 8767 OF 15495 ***
    // Wavefunction(s) for diagram number 8767
    // (none)
    // Amplitude(s) for diagram number 8767
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[218], w_fp[505], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[220] -= amp_sv[0];

    // *** DIAGRAM 8768 OF 15495 ***
    // Wavefunction(s) for diagram number 8768
    // (none)
    // Amplitude(s) for diagram number 8768
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[481], w_fp[441], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[234] -= amp_sv[0];

    // *** DIAGRAM 8769 OF 15495 ***
    // Wavefunction(s) for diagram number 8769
    // (none)
    // Amplitude(s) for diagram number 8769
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[170], w_fp[505], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[217] -= amp_sv[0];

    // *** DIAGRAM 8770 OF 15495 ***
    // Wavefunction(s) for diagram number 8770
    // (none)
    // Amplitude(s) for diagram number 8770
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[535], w_fp[441], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[223] -= amp_sv[0];

    // *** DIAGRAM 8771 OF 15495 ***
    // Wavefunction(s) for diagram number 8771
    // (none)
    // Amplitude(s) for diagram number 8771
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[200], w_fp[449], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[139] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8772 OF 15495 ***
    // Wavefunction(s) for diagram number 8772
    // (none)
    // Amplitude(s) for diagram number 8772
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[171], w_fp[449], w_fp[86], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8773 OF 15495 ***
    // Wavefunction(s) for diagram number 8773
    // (none)
    // Amplitude(s) for diagram number 8773
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[168], w_fp[449], w_fp[88], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];

    // *** DIAGRAM 8774 OF 15495 ***
    // Wavefunction(s) for diagram number 8774
    // (none)
    // Amplitude(s) for diagram number 8774
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[527], w_fp[460], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8775 OF 15495 ***
    // Wavefunction(s) for diagram number 8775
    // (none)
    // Amplitude(s) for diagram number 8775
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[527], w_fp[441], w_fp[86], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8776 OF 15495 ***
    // Wavefunction(s) for diagram number 8776
    // (none)
    // Amplitude(s) for diagram number 8776
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[527], w_fp[259], w_fp[88], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] += amp_sv[0];
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];

    // *** DIAGRAM 8777 OF 15495 ***
    // Wavefunction(s) for diagram number 8777
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[86], COUPs[0], 1.0, 0., 0., w_fp[576] );
    // Amplitude(s) for diagram number 8777
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[576], w_fp[468], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[217] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8778 OF 15495 ***
    // Wavefunction(s) for diagram number 8778
    // (none)
    // Amplitude(s) for diagram number 8778
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[168], w_fp[441], w_fp[576], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[217] += amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];

    // *** DIAGRAM 8779 OF 15495 ***
    // Wavefunction(s) for diagram number 8779
    // (none)
    // Amplitude(s) for diagram number 8779
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[171], w_fp[259], w_fp[576], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];

    // *** DIAGRAM 8780 OF 15495 ***
    // Wavefunction(s) for diagram number 8780
    // (none)
    // Amplitude(s) for diagram number 8780
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[468], w_fp[88], COUPs[0], 1.0, &amp_fp[0] );
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

    // *** DIAGRAM 8781 OF 15495 ***
    // Wavefunction(s) for diagram number 8781
    // (none)
    // Amplitude(s) for diagram number 8781
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[171], w_fp[460], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[157] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8782 OF 15495 ***
    // Wavefunction(s) for diagram number 8782
    // (none)
    // Amplitude(s) for diagram number 8782
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[200], w_fp[441], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[217] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8783 OF 15495 ***
    // Wavefunction(s) for diagram number 8783
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[86], w_fp[7], COUPs[2], 1.0, 0., 0., w_fp[550] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[86], w_fp[7], COUPs[2], 1.0, 0., 0., w_fp[480] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[86], w_fp[7], COUPs[2], 1.0, 0., 0., w_fp[558] );
    // Amplitude(s) for diagram number 8783
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[168], w_fp[259], w_fp[550], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[168], w_fp[259], w_fp[480], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[139] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[217] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[168], w_fp[259], w_fp[558], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[217] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8784 OF 15495 ***
    // Wavefunction(s) for diagram number 8784
    // (none)
    // Amplitude(s) for diagram number 8784
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[99], w_fp[449], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[133] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8785 OF 15495 ***
    // Wavefunction(s) for diagram number 8785
    // (none)
    // Amplitude(s) for diagram number 8785
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[170], w_fp[449], w_fp[102], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[125] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8786 OF 15495 ***
    // Wavefunction(s) for diagram number 8786
    // (none)
    // Amplitude(s) for diagram number 8786
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[168], w_fp[449], w_fp[104], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[125] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[136] += amp_sv[0];
    jamp_sv[139] -= amp_sv[0];

    // *** DIAGRAM 8787 OF 15495 ***
    // Wavefunction(s) for diagram number 8787
    // (none)
    // Amplitude(s) for diagram number 8787
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[527], w_fp[463], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8788 OF 15495 ***
    // Wavefunction(s) for diagram number 8788
    // (none)
    // Amplitude(s) for diagram number 8788
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[527], w_fp[439], w_fp[102], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8789 OF 15495 ***
    // Wavefunction(s) for diagram number 8789
    // (none)
    // Amplitude(s) for diagram number 8789
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[527], w_fp[259], w_fp[104], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[166] += amp_sv[0];
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[226] -= amp_sv[0];

    // *** DIAGRAM 8790 OF 15495 ***
    // Wavefunction(s) for diagram number 8790
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[102], COUPs[0], 1.0, 0., 0., w_fp[518] );
    // Amplitude(s) for diagram number 8790
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[518], w_fp[468], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[125] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[163] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[196] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8791 OF 15495 ***
    // Wavefunction(s) for diagram number 8791
    // (none)
    // Amplitude(s) for diagram number 8791
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[168], w_fp[439], w_fp[518], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[193] += amp_sv[0];
    jamp_sv[196] -= amp_sv[0];
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];

    // *** DIAGRAM 8792 OF 15495 ***
    // Wavefunction(s) for diagram number 8792
    // (none)
    // Amplitude(s) for diagram number 8792
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[170], w_fp[259], w_fp[518], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[125] += amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[163] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];

    // *** DIAGRAM 8793 OF 15495 ***
    // Wavefunction(s) for diagram number 8793
    // (none)
    // Amplitude(s) for diagram number 8793
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[468], w_fp[104], COUPs[0], 1.0, &amp_fp[0] );
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

    // *** DIAGRAM 8794 OF 15495 ***
    // Wavefunction(s) for diagram number 8794
    // (none)
    // Amplitude(s) for diagram number 8794
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[170], w_fp[463], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[163] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8795 OF 15495 ***
    // Wavefunction(s) for diagram number 8795
    // (none)
    // Amplitude(s) for diagram number 8795
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[99], w_fp[439], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[193] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[196] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8796 OF 15495 ***
    // Wavefunction(s) for diagram number 8796
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[102], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[486] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[102], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[547] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[102], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[529] );
    // Amplitude(s) for diagram number 8796
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[168], w_fp[259], w_fp[486], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[125] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[168], w_fp[259], w_fp[547], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[133] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[163] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[196] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[168], w_fp[259], w_fp[529], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[125] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[163] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[196] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8797 OF 15495 ***
    // Wavefunction(s) for diagram number 8797
    // (none)
    // Amplitude(s) for diagram number 8797
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[218], w_fp[449], w_fp[84], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[136] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8798 OF 15495 ***
    // Wavefunction(s) for diagram number 8798
    // (none)
    // Amplitude(s) for diagram number 8798
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[112], w_fp[449], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8799 OF 15495 ***
    // Wavefunction(s) for diagram number 8799
    // (none)
    // Amplitude(s) for diagram number 8799
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[168], w_fp[449], w_fp[130], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];

    // *** DIAGRAM 8800 OF 15495 ***
    // Wavefunction(s) for diagram number 8800
    // (none)
    // Amplitude(s) for diagram number 8800
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[527], w_fp[443], w_fp[84], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 437 );
    storeWf( wfs, w_cx, nevt, 444 );
    storeWf( wfs, w_cx, nevt, 450 );
    storeWf( wfs, w_cx, nevt, 451 );
    storeWf( wfs, w_cx, nevt, 452 );
    storeWf( wfs, w_cx, nevt, 474 );
    storeWf( wfs, w_cx, nevt, 479 );
    storeWf( wfs, w_cx, nevt, 480 );
    storeWf( wfs, w_cx, nevt, 481 );
    storeWf( wfs, w_cx, nevt, 486 );
    storeWf( wfs, w_cx, nevt, 488 );
    storeWf( wfs, w_cx, nevt, 505 );
    storeWf( wfs, w_cx, nevt, 509 );
    storeWf( wfs, w_cx, nevt, 516 );
    storeWf( wfs, w_cx, nevt, 518 );
    storeWf( wfs, w_cx, nevt, 520 );
    storeWf( wfs, w_cx, nevt, 522 );
    storeWf( wfs, w_cx, nevt, 523 );
    storeWf( wfs, w_cx, nevt, 527 );
    storeWf( wfs, w_cx, nevt, 528 );
    storeWf( wfs, w_cx, nevt, 529 );
    storeWf( wfs, w_cx, nevt, 535 );
    storeWf( wfs, w_cx, nevt, 538 );
    storeWf( wfs, w_cx, nevt, 543 );
    storeWf( wfs, w_cx, nevt, 546 );
    storeWf( wfs, w_cx, nevt, 547 );
    storeWf( wfs, w_cx, nevt, 549 );
    storeWf( wfs, w_cx, nevt, 550 );
    storeWf( wfs, w_cx, nevt, 558 );
    storeWf( wfs, w_cx, nevt, 559 );
    storeWf( wfs, w_cx, nevt, 571 );
    storeWf( wfs, w_cx, nevt, 576 );
    storeWf( wfs, w_cx, nevt, 580 );
#endif
#endif
  }

  //--------------------------------------------------------------------------
}
