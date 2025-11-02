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
  diagramgroup261( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 45 );
    retrieveWf( wfs, w_cx, nevt, 151 );
    retrieveWf( wfs, w_cx, nevt, 152 );
    retrieveWf( wfs, w_cx, nevt, 153 );
    retrieveWf( wfs, w_cx, nevt, 259 );
    retrieveWf( wfs, w_cx, nevt, 260 );
    retrieveWf( wfs, w_cx, nevt, 263 );
    retrieveWf( wfs, w_cx, nevt, 265 );
    retrieveWf( wfs, w_cx, nevt, 435 );
    retrieveWf( wfs, w_cx, nevt, 455 );
#endif
#endif

    // *** DIAGRAM 2601 OF 15495 ***
    // Wavefunction(s) for diagram number 2601
    // (none)
    // Amplitude(s) for diagram number 2601
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[4], w_fp[151], w_fp[455], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[155] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[185] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[4], w_fp[152], w_fp[455], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[159] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[165] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[209] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[4], w_fp[153], w_fp[455], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[159] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[165] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[185] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[209] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2602 OF 15495 ***
    // Wavefunction(s) for diagram number 2602
    // (none)
    // Amplitude(s) for diagram number 2602
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[435], w_fp[259], w_fp[151], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] += amp_sv[0];
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[435], w_fp[259], w_fp[152], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[209] += amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[233] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[435], w_fp[259], w_fp[153], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[209] += amp_sv[0];
    jamp_sv[233] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];

    // *** DIAGRAM 2603 OF 15495 ***
    // Wavefunction(s) for diagram number 2603
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[45], w_fp[151], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[455] );
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[45], w_fp[152], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[477] );
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[45], w_fp[153], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[478] );
    // Amplitude(s) for diagram number 2603
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[455], w_fp[259], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += amp_sv[0];
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[167] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[477], w_fp[259], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[159] += amp_sv[0];
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[478], w_fp[259], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[159] += amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[167] -= amp_sv[0];

    // *** DIAGRAM 2604 OF 15495 ***
    // Wavefunction(s) for diagram number 2604
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[45], w_fp[2], COUPs[1], 1.0, 0., 0., w_fp[479] );
    // Amplitude(s) for diagram number 2604
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[260], w_fp[6], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[155] += amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[393] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[417] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    jamp_sv[705] -= amp_sv[0];
    jamp_sv[711] += amp_sv[0];
    jamp_sv[717] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[260], w_fp[6], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[260], w_fp[6], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];

    // *** DIAGRAM 2605 OF 15495 ***
    // Wavefunction(s) for diagram number 2605
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[479], w_fp[6], COUPs[0], 1.0, 0., 0., w_fp[480] );
    // Amplitude(s) for diagram number 2605
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[260], w_fp[7], w_fp[480], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];

    // *** DIAGRAM 2606 OF 15495 ***
    // Wavefunction(s) for diagram number 2606
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[479], w_fp[7], COUPs[0], 1.0, 0., 0., w_fp[481] );
    // Amplitude(s) for diagram number 2606
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[260], w_fp[6], w_fp[481], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];

    // *** DIAGRAM 2607 OF 15495 ***
    // Wavefunction(s) for diagram number 2607
    // (none)
    // Amplitude(s) for diagram number 2607
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[263], w_fp[5], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] -= amp_sv[0];
    jamp_sv[161] += amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    jamp_sv[513] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[537] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[681] -= amp_sv[0];
    jamp_sv[687] += amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    jamp_sv[695] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[263], w_fp[5], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[263], w_fp[5], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];

    // *** DIAGRAM 2608 OF 15495 ***
    // Wavefunction(s) for diagram number 2608
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[479], w_fp[5], COUPs[0], 1.0, 0., 0., w_fp[482] );
    // Amplitude(s) for diagram number 2608
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[263], w_fp[7], w_fp[482], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];

    // *** DIAGRAM 2609 OF 15495 ***
    // Wavefunction(s) for diagram number 2609
    // (none)
    // Amplitude(s) for diagram number 2609
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[263], w_fp[5], w_fp[481], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];

    // *** DIAGRAM 2610 OF 15495 ***
    // Wavefunction(s) for diagram number 2610
    // (none)
    // Amplitude(s) for diagram number 2610
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[265], w_fp[5], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] -= amp_sv[0];
    jamp_sv[167] += amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[561] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[633] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[657] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[265], w_fp[5], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[167] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[265], w_fp[5], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] += amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[657] += amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 479 );
    storeWf( wfs, w_cx, nevt, 480 );
    storeWf( wfs, w_cx, nevt, 481 );
    storeWf( wfs, w_cx, nevt, 482 );
#endif
#endif
  }

  //--------------------------------------------------------------------------

#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup262( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 258 );
    retrieveWf( wfs, w_cx, nevt, 263 );
    retrieveWf( wfs, w_cx, nevt, 265 );
    retrieveWf( wfs, w_cx, nevt, 266 );
    retrieveWf( wfs, w_cx, nevt, 267 );
    retrieveWf( wfs, w_cx, nevt, 268 );
    retrieveWf( wfs, w_cx, nevt, 269 );
    retrieveWf( wfs, w_cx, nevt, 270 );
    retrieveWf( wfs, w_cx, nevt, 271 );
    retrieveWf( wfs, w_cx, nevt, 272 );
    retrieveWf( wfs, w_cx, nevt, 273 );
    retrieveWf( wfs, w_cx, nevt, 274 );
    retrieveWf( wfs, w_cx, nevt, 442 );
    retrieveWf( wfs, w_cx, nevt, 444 );
    retrieveWf( wfs, w_cx, nevt, 445 );
    retrieveWf( wfs, w_cx, nevt, 479 );
    retrieveWf( wfs, w_cx, nevt, 480 );
    retrieveWf( wfs, w_cx, nevt, 482 );
#endif
#endif

    // *** DIAGRAM 2611 OF 15495 ***
    // Wavefunction(s) for diagram number 2611
    // (none)
    // Amplitude(s) for diagram number 2611
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[265], w_fp[6], w_fp[482], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[167] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];

    // *** DIAGRAM 2612 OF 15495 ***
    // Wavefunction(s) for diagram number 2612
    // (none)
    // Amplitude(s) for diagram number 2612
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[265], w_fp[5], w_fp[480], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] += amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[657] += amp_sv[0];

    // *** DIAGRAM 2613 OF 15495 ***
    // Wavefunction(s) for diagram number 2613
    // (none)
    // Amplitude(s) for diagram number 2613
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[266], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[159] += amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[561] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[633] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[657] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[717] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[267], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[268], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];

    // *** DIAGRAM 2614 OF 15495 ***
    // Wavefunction(s) for diagram number 2614
    // (none)
    // Amplitude(s) for diagram number 2614
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[269], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[275] += amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    jamp_sv[513] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[537] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[597] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];
    jamp_sv[681] -= amp_sv[0];
    jamp_sv[687] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[270], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] += amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[271], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];

    // *** DIAGRAM 2615 OF 15495 ***
    // Wavefunction(s) for diagram number 2615
    // (none)
    // Amplitude(s) for diagram number 2615
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[272], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[167] += amp_sv[0];
    jamp_sv[281] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[393] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[417] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[477] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[705] -= amp_sv[0];
    jamp_sv[711] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[273], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[167] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[274], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];

    // *** DIAGRAM 2616 OF 15495 ***
    // Wavefunction(s) for diagram number 2616
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[258], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[483] );
    // Amplitude(s) for diagram number 2616
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[444], w_fp[483], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2617 OF 15495 ***
    // Wavefunction(s) for diagram number 2617
    // (none)
    // Amplitude(s) for diagram number 2617
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[445], w_fp[483], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2618 OF 15495 ***
    // Wavefunction(s) for diagram number 2618
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[442], w_fp[2], COUPs[1], 1.0, 0., 0., w_fp[484] );
    // Amplitude(s) for diagram number 2618
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[263], w_fp[7], w_fp[484], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2619 OF 15495 ***
    // Wavefunction(s) for diagram number 2619
    // (none)
    // Amplitude(s) for diagram number 2619
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[445], w_fp[2], w_fp[263], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];

    // *** DIAGRAM 2620 OF 15495 ***
    // Wavefunction(s) for diagram number 2620
    // (none)
    // Amplitude(s) for diagram number 2620
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[265], w_fp[6], w_fp[484], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] -= cxtype( 0, 1 ) * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 483 );
    storeWf( wfs, w_cx, nevt, 484 );
#endif
#endif
  }

  //--------------------------------------------------------------------------


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup263( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 260 );
    retrieveWf( wfs, w_cx, nevt, 265 );
    retrieveWf( wfs, w_cx, nevt, 269 );
    retrieveWf( wfs, w_cx, nevt, 270 );
    retrieveWf( wfs, w_cx, nevt, 271 );
    retrieveWf( wfs, w_cx, nevt, 272 );
    retrieveWf( wfs, w_cx, nevt, 273 );
    retrieveWf( wfs, w_cx, nevt, 274 );
    retrieveWf( wfs, w_cx, nevt, 442 );
    retrieveWf( wfs, w_cx, nevt, 444 );
    retrieveWf( wfs, w_cx, nevt, 447 );
    retrieveWf( wfs, w_cx, nevt, 448 );
    retrieveWf( wfs, w_cx, nevt, 449 );
    retrieveWf( wfs, w_cx, nevt, 452 );
    retrieveWf( wfs, w_cx, nevt, 483 );
#endif
#endif

    // *** DIAGRAM 2621 OF 15495 ***
    // Wavefunction(s) for diagram number 2621
    // (none)
    // Amplitude(s) for diagram number 2621
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[444], w_fp[2], w_fp[265], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[167] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];

    // *** DIAGRAM 2622 OF 15495 ***
    // Wavefunction(s) for diagram number 2622
    // (none)
    // Amplitude(s) for diagram number 2622
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[442], w_fp[2], w_fp[272], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[442], w_fp[2], w_fp[273], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[442], w_fp[2], w_fp[274], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2623 OF 15495 ***
    // Wavefunction(s) for diagram number 2623
    // (none)
    // Amplitude(s) for diagram number 2623
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[448], w_fp[483], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2624 OF 15495 ***
    // Wavefunction(s) for diagram number 2624
    // (none)
    // Amplitude(s) for diagram number 2624
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[449], w_fp[483], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2625 OF 15495 ***
    // Wavefunction(s) for diagram number 2625
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[447], w_fp[2], COUPs[1], 1.0, 0., 0., w_fp[485] );
    // Amplitude(s) for diagram number 2625
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[260], w_fp[7], w_fp[485], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[419] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[633] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[657] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[687] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2626 OF 15495 ***
    // Wavefunction(s) for diagram number 2626
    // (none)
    // Amplitude(s) for diagram number 2626
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[449], w_fp[2], w_fp[260], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];

    // *** DIAGRAM 2627 OF 15495 ***
    // Wavefunction(s) for diagram number 2627
    // (none)
    // Amplitude(s) for diagram number 2627
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[265], w_fp[5], w_fp[485], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[419] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[471] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[633] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[657] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2628 OF 15495 ***
    // Wavefunction(s) for diagram number 2628
    // (none)
    // Amplitude(s) for diagram number 2628
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[448], w_fp[2], w_fp[265], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] += amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[657] += amp_sv[0];

    // *** DIAGRAM 2629 OF 15495 ***
    // Wavefunction(s) for diagram number 2629
    // (none)
    // Amplitude(s) for diagram number 2629
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[447], w_fp[2], w_fp[269], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[165] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[471] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[687] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[447], w_fp[2], w_fp[270], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[419] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[471] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[633] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[657] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[447], w_fp[2], w_fp[271], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[419] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[633] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[657] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[687] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2630 OF 15495 ***
    // Wavefunction(s) for diagram number 2630
    // (none)
    // Amplitude(s) for diagram number 2630
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[452], w_fp[483], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 485 );
#endif
#endif
  }

  //--------------------------------------------------------------------------


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup264( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 45 );
    retrieveWf( wfs, w_cx, nevt, 169 );
    retrieveWf( wfs, w_cx, nevt, 191 );
    retrieveWf( wfs, w_cx, nevt, 193 );
    retrieveWf( wfs, w_cx, nevt, 258 );
    retrieveWf( wfs, w_cx, nevt, 260 );
    retrieveWf( wfs, w_cx, nevt, 263 );
    retrieveWf( wfs, w_cx, nevt, 265 );
    retrieveWf( wfs, w_cx, nevt, 266 );
    retrieveWf( wfs, w_cx, nevt, 267 );
    retrieveWf( wfs, w_cx, nevt, 268 );
    retrieveWf( wfs, w_cx, nevt, 451 );
    retrieveWf( wfs, w_cx, nevt, 452 );
    retrieveWf( wfs, w_cx, nevt, 453 );
    retrieveWf( wfs, w_cx, nevt, 483 );
#endif
#endif

    // *** DIAGRAM 2631 OF 15495 ***
    // Wavefunction(s) for diagram number 2631
    // (none)
    // Amplitude(s) for diagram number 2631
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[453], w_fp[483], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2632 OF 15495 ***
    // Wavefunction(s) for diagram number 2632
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[451], w_fp[2], COUPs[1], 1.0, 0., 0., w_fp[486] );
    // Amplitude(s) for diagram number 2632
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[260], w_fp[6], w_fp[486], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2633 OF 15495 ***
    // Wavefunction(s) for diagram number 2633
    // (none)
    // Amplitude(s) for diagram number 2633
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[453], w_fp[2], w_fp[260], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[417] += amp_sv[0];

    // *** DIAGRAM 2634 OF 15495 ***
    // Wavefunction(s) for diagram number 2634
    // (none)
    // Amplitude(s) for diagram number 2634
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[263], w_fp[5], w_fp[486], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2635 OF 15495 ***
    // Wavefunction(s) for diagram number 2635
    // (none)
    // Amplitude(s) for diagram number 2635
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[452], w_fp[2], w_fp[263], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[537] += amp_sv[0];

    // *** DIAGRAM 2636 OF 15495 ***
    // Wavefunction(s) for diagram number 2636
    // (none)
    // Amplitude(s) for diagram number 2636
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[2], w_fp[266], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[159] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[2], w_fp[267], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[2], w_fp[268], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2637 OF 15495 ***
    // Wavefunction(s) for diagram number 2637
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[45], w_fp[258], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[487] );
    // Amplitude(s) for diagram number 2637
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[487], w_fp[191], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[453] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2638 OF 15495 ***
    // Wavefunction(s) for diagram number 2638
    // (none)
    // Amplitude(s) for diagram number 2638
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[487], w_fp[193], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[477] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2639 OF 15495 ***
    // Wavefunction(s) for diagram number 2639
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[45], w_fp[169], COUPs[1], 1.0, 0., 0., w_fp[488] );
    // Amplitude(s) for diagram number 2639
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[263], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[471] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[477] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2640 OF 15495 ***
    // Wavefunction(s) for diagram number 2640
    // (none)
    // Amplitude(s) for diagram number 2640
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[265], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[395] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[419] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[453] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[471] -= cxtype( 0, 1 ) * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 486 );
    storeWf( wfs, w_cx, nevt, 487 );
    storeWf( wfs, w_cx, nevt, 488 );
#endif
#endif
  }

  //--------------------------------------------------------------------------


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup265( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 45 );
    retrieveWf( wfs, w_cx, nevt, 84 );
    retrieveWf( wfs, w_cx, nevt, 169 );
    retrieveWf( wfs, w_cx, nevt, 191 );
    retrieveWf( wfs, w_cx, nevt, 193 );
    retrieveWf( wfs, w_cx, nevt, 258 );
    retrieveWf( wfs, w_cx, nevt, 263 );
    retrieveWf( wfs, w_cx, nevt, 265 );
    retrieveWf( wfs, w_cx, nevt, 447 );
    retrieveWf( wfs, w_cx, nevt, 451 );
    retrieveWf( wfs, w_cx, nevt, 487 );
    retrieveWf( wfs, w_cx, nevt, 488 );
#endif
#endif

    // *** DIAGRAM 2641 OF 15495 ***
    // Wavefunction(s) for diagram number 2641
    // (none)
    // Amplitude(s) for diagram number 2641
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[258], w_fp[6], w_fp[7], w_fp[488], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[471] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[477] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[258], w_fp[6], w_fp[7], w_fp[488], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[395] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[419] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[453] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[471] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[258], w_fp[6], w_fp[7], w_fp[488], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[393] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[419] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[453] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[477] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2642 OF 15495 ***
    // Wavefunction(s) for diagram number 2642
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[169], w_fp[258], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[489] );
    // Amplitude(s) for diagram number 2642
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[447], w_fp[489], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[395] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[419] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2643 OF 15495 ***
    // Wavefunction(s) for diagram number 2643
    // (none)
    // Amplitude(s) for diagram number 2643
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[447], w_fp[169], w_fp[265], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[395] += amp_sv[0];
    jamp_sv[419] -= amp_sv[0];
    jamp_sv[465] -= amp_sv[0];
    jamp_sv[471] += amp_sv[0];

    // *** DIAGRAM 2644 OF 15495 ***
    // Wavefunction(s) for diagram number 2644
    // (none)
    // Amplitude(s) for diagram number 2644
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[447], w_fp[193], w_fp[258], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[465] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[471] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2645 OF 15495 ***
    // Wavefunction(s) for diagram number 2645
    // (none)
    // Amplitude(s) for diagram number 2645
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[489], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[393] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2646 OF 15495 ***
    // Wavefunction(s) for diagram number 2646
    // (none)
    // Amplitude(s) for diagram number 2646
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[169], w_fp[263], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[393] += amp_sv[0];
    jamp_sv[417] -= amp_sv[0];
    jamp_sv[441] -= amp_sv[0];
    jamp_sv[447] += amp_sv[0];

    // *** DIAGRAM 2647 OF 15495 ***
    // Wavefunction(s) for diagram number 2647
    // (none)
    // Amplitude(s) for diagram number 2647
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[191], w_fp[258], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2648 OF 15495 ***
    // Wavefunction(s) for diagram number 2648
    // (none)
    // Amplitude(s) for diagram number 2648
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[45], w_fp[193], w_fp[263], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[465] += amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];

    // *** DIAGRAM 2649 OF 15495 ***
    // Wavefunction(s) for diagram number 2649
    // (none)
    // Amplitude(s) for diagram number 2649
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[45], w_fp[191], w_fp[265], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[441] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];

    // *** DIAGRAM 2650 OF 15495 ***
    // Wavefunction(s) for diagram number 2650
    // (none)
    // Amplitude(s) for diagram number 2650
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[487], w_fp[169], w_fp[84], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[453] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 489 );
#endif
#endif
  }

  //--------------------------------------------------------------------------


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup266( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 45 );
    retrieveWf( wfs, w_cx, nevt, 84 );
    retrieveWf( wfs, w_cx, nevt, 169 );
    retrieveWf( wfs, w_cx, nevt, 197 );
    retrieveWf( wfs, w_cx, nevt, 212 );
    retrieveWf( wfs, w_cx, nevt, 213 );
    retrieveWf( wfs, w_cx, nevt, 258 );
    retrieveWf( wfs, w_cx, nevt, 260 );
    retrieveWf( wfs, w_cx, nevt, 265 );
    retrieveWf( wfs, w_cx, nevt, 442 );
    retrieveWf( wfs, w_cx, nevt, 458 );
    retrieveWf( wfs, w_cx, nevt, 487 );
    retrieveWf( wfs, w_cx, nevt, 488 );
#endif
#endif

    // *** DIAGRAM 2651 OF 15495 ***
    // Wavefunction(s) for diagram number 2651
    // (none)
    // Amplitude(s) for diagram number 2651
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[258], w_fp[84], w_fp[488], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[419] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[453] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[477] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2652 OF 15495 ***
    // Wavefunction(s) for diagram number 2652
    // (none)
    // Amplitude(s) for diagram number 2652
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[458], w_fp[169], w_fp[258], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[393] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[417] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];

    // *** DIAGRAM 2653 OF 15495 ***
    // Wavefunction(s) for diagram number 2653
    // (none)
    // Amplitude(s) for diagram number 2653
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[487], w_fp[212], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[573] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2654 OF 15495 ***
    // Wavefunction(s) for diagram number 2654
    // (none)
    // Amplitude(s) for diagram number 2654
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[487], w_fp[213], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[597] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2655 OF 15495 ***
    // Wavefunction(s) for diagram number 2655
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[45], w_fp[197], COUPs[1], 1.0, 0., 0., w_fp[490] );
    // Amplitude(s) for diagram number 2655
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[490], w_fp[260], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[513] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2656 OF 15495 ***
    // Wavefunction(s) for diagram number 2656
    // (none)
    // Amplitude(s) for diagram number 2656
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[490], w_fp[265], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[515] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2657 OF 15495 ***
    // Wavefunction(s) for diagram number 2657
    // (none)
    // Amplitude(s) for diagram number 2657
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[258], w_fp[5], w_fp[7], w_fp[490], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[513] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[258], w_fp[5], w_fp[7], w_fp[490], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[515] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[258], w_fp[5], w_fp[7], w_fp[490], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[513] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2658 OF 15495 ***
    // Wavefunction(s) for diagram number 2658
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[197], w_fp[258], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[491] );
    // Amplitude(s) for diagram number 2658
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[442], w_fp[491], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2659 OF 15495 ***
    // Wavefunction(s) for diagram number 2659
    // (none)
    // Amplitude(s) for diagram number 2659
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[442], w_fp[197], w_fp[265], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[515] += amp_sv[0];
    jamp_sv[539] -= amp_sv[0];
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[591] += amp_sv[0];

    // *** DIAGRAM 2660 OF 15495 ***
    // Wavefunction(s) for diagram number 2660
    // (none)
    // Amplitude(s) for diagram number 2660
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[442], w_fp[213], w_fp[258], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] -= cxtype( 0, 1 ) * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 490 );
    storeWf( wfs, w_cx, nevt, 491 );
#endif
#endif
  }

  //--------------------------------------------------------------------------


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup267( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 45 );
    retrieveWf( wfs, w_cx, nevt, 100 );
    retrieveWf( wfs, w_cx, nevt, 197 );
    retrieveWf( wfs, w_cx, nevt, 212 );
    retrieveWf( wfs, w_cx, nevt, 213 );
    retrieveWf( wfs, w_cx, nevt, 226 );
    retrieveWf( wfs, w_cx, nevt, 227 );
    retrieveWf( wfs, w_cx, nevt, 258 );
    retrieveWf( wfs, w_cx, nevt, 260 );
    retrieveWf( wfs, w_cx, nevt, 265 );
    retrieveWf( wfs, w_cx, nevt, 451 );
    retrieveWf( wfs, w_cx, nevt, 461 );
    retrieveWf( wfs, w_cx, nevt, 487 );
    retrieveWf( wfs, w_cx, nevt, 490 );
    retrieveWf( wfs, w_cx, nevt, 491 );
#endif
#endif

    // *** DIAGRAM 2661 OF 15495 ***
    // Wavefunction(s) for diagram number 2661
    // (none)
    // Amplitude(s) for diagram number 2661
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[491], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[513] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2662 OF 15495 ***
    // Wavefunction(s) for diagram number 2662
    // (none)
    // Amplitude(s) for diagram number 2662
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[197], w_fp[260], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[513] += amp_sv[0];
    jamp_sv[537] -= amp_sv[0];
    jamp_sv[561] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];

    // *** DIAGRAM 2663 OF 15495 ***
    // Wavefunction(s) for diagram number 2663
    // (none)
    // Amplitude(s) for diagram number 2663
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[212], w_fp[258], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2664 OF 15495 ***
    // Wavefunction(s) for diagram number 2664
    // (none)
    // Amplitude(s) for diagram number 2664
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[45], w_fp[213], w_fp[260], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[585] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];

    // *** DIAGRAM 2665 OF 15495 ***
    // Wavefunction(s) for diagram number 2665
    // (none)
    // Amplitude(s) for diagram number 2665
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[45], w_fp[212], w_fp[265], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[561] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];

    // *** DIAGRAM 2666 OF 15495 ***
    // Wavefunction(s) for diagram number 2666
    // (none)
    // Amplitude(s) for diagram number 2666
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[487], w_fp[197], w_fp[100], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[573] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];

    // *** DIAGRAM 2667 OF 15495 ***
    // Wavefunction(s) for diagram number 2667
    // (none)
    // Amplitude(s) for diagram number 2667
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[258], w_fp[100], w_fp[490], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[513] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2668 OF 15495 ***
    // Wavefunction(s) for diagram number 2668
    // (none)
    // Amplitude(s) for diagram number 2668
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[461], w_fp[197], w_fp[258], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[513] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[537] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];

    // *** DIAGRAM 2669 OF 15495 ***
    // Wavefunction(s) for diagram number 2669
    // (none)
    // Amplitude(s) for diagram number 2669
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[487], w_fp[226], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[693] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2670 OF 15495 ***
    // Wavefunction(s) for diagram number 2670
    // (none)
    // Amplitude(s) for diagram number 2670
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[487], w_fp[227], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[717] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] -= cxtype( 0, 1 ) * amp_sv[0];

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


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup268( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 45 );
    retrieveWf( wfs, w_cx, nevt, 215 );
    retrieveWf( wfs, w_cx, nevt, 226 );
    retrieveWf( wfs, w_cx, nevt, 227 );
    retrieveWf( wfs, w_cx, nevt, 258 );
    retrieveWf( wfs, w_cx, nevt, 260 );
    retrieveWf( wfs, w_cx, nevt, 263 );
    retrieveWf( wfs, w_cx, nevt, 442 );
    retrieveWf( wfs, w_cx, nevt, 447 );
#endif
#endif

    // *** DIAGRAM 2671 OF 15495 ***
    // Wavefunction(s) for diagram number 2671
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[45], w_fp[215], COUPs[1], 1.0, 0., 0., w_fp[492] );
    // Amplitude(s) for diagram number 2671
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[492], w_fp[260], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[633] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[657] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[687] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2672 OF 15495 ***
    // Wavefunction(s) for diagram number 2672
    // (none)
    // Amplitude(s) for diagram number 2672
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[492], w_fp[263], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[635] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[687] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2673 OF 15495 ***
    // Wavefunction(s) for diagram number 2673
    // (none)
    // Amplitude(s) for diagram number 2673
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[258], w_fp[5], w_fp[6], w_fp[492], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[633] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[657] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[687] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[258], w_fp[5], w_fp[6], w_fp[492], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[635] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[687] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[258], w_fp[5], w_fp[6], w_fp[492], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[633] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[657] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2674 OF 15495 ***
    // Wavefunction(s) for diagram number 2674
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[215], w_fp[258], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[493] );
    // Amplitude(s) for diagram number 2674
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[442], w_fp[493], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[635] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2675 OF 15495 ***
    // Wavefunction(s) for diagram number 2675
    // (none)
    // Amplitude(s) for diagram number 2675
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[442], w_fp[215], w_fp[263], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[635] += amp_sv[0];
    jamp_sv[659] -= amp_sv[0];
    jamp_sv[705] -= amp_sv[0];
    jamp_sv[711] += amp_sv[0];

    // *** DIAGRAM 2676 OF 15495 ***
    // Wavefunction(s) for diagram number 2676
    // (none)
    // Amplitude(s) for diagram number 2676
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[442], w_fp[227], w_fp[258], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2677 OF 15495 ***
    // Wavefunction(s) for diagram number 2677
    // (none)
    // Amplitude(s) for diagram number 2677
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[447], w_fp[493], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[633] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[657] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2678 OF 15495 ***
    // Wavefunction(s) for diagram number 2678
    // (none)
    // Amplitude(s) for diagram number 2678
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[447], w_fp[215], w_fp[260], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[633] += amp_sv[0];
    jamp_sv[657] -= amp_sv[0];
    jamp_sv[681] -= amp_sv[0];
    jamp_sv[687] += amp_sv[0];

    // *** DIAGRAM 2679 OF 15495 ***
    // Wavefunction(s) for diagram number 2679
    // (none)
    // Amplitude(s) for diagram number 2679
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[447], w_fp[226], w_fp[258], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[681] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[687] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2680 OF 15495 ***
    // Wavefunction(s) for diagram number 2680
    // (none)
    // Amplitude(s) for diagram number 2680
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[45], w_fp[227], w_fp[260], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[705] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 492 );
    storeWf( wfs, w_cx, nevt, 493 );
#endif
#endif
  }

  //--------------------------------------------------------------------------


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup269( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 45 );
    retrieveWf( wfs, w_cx, nevt, 113 );
    retrieveWf( wfs, w_cx, nevt, 116 );
    retrieveWf( wfs, w_cx, nevt, 215 );
    retrieveWf( wfs, w_cx, nevt, 226 );
    retrieveWf( wfs, w_cx, nevt, 244 );
    retrieveWf( wfs, w_cx, nevt, 258 );
    retrieveWf( wfs, w_cx, nevt, 263 );
    retrieveWf( wfs, w_cx, nevt, 265 );
    retrieveWf( wfs, w_cx, nevt, 286 );
    retrieveWf( wfs, w_cx, nevt, 464 );
    retrieveWf( wfs, w_cx, nevt, 479 );
    retrieveWf( wfs, w_cx, nevt, 487 );
    retrieveWf( wfs, w_cx, nevt, 492 );
#endif
#endif

    // *** DIAGRAM 2681 OF 15495 ***
    // Wavefunction(s) for diagram number 2681
    // (none)
    // Amplitude(s) for diagram number 2681
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[45], w_fp[226], w_fp[263], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[681] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];

    // *** DIAGRAM 2682 OF 15495 ***
    // Wavefunction(s) for diagram number 2682
    // (none)
    // Amplitude(s) for diagram number 2682
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[487], w_fp[215], w_fp[113], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[693] += amp_sv[0];
    jamp_sv[695] -= amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];

    // *** DIAGRAM 2683 OF 15495 ***
    // Wavefunction(s) for diagram number 2683
    // (none)
    // Amplitude(s) for diagram number 2683
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[258], w_fp[113], w_fp[492], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[633] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[657] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2684 OF 15495 ***
    // Wavefunction(s) for diagram number 2684
    // (none)
    // Amplitude(s) for diagram number 2684
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[464], w_fp[215], w_fp[258], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[633] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[657] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];

    // *** DIAGRAM 2685 OF 15495 ***
    // Wavefunction(s) for diagram number 2685
    // (none)
    // Amplitude(s) for diagram number 2685
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[487], w_fp[244], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[453] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];

    // *** DIAGRAM 2686 OF 15495 ***
    // Wavefunction(s) for diagram number 2686
    // (none)
    // Amplitude(s) for diagram number 2686
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[487], w_fp[2], w_fp[116], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[453] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2687 OF 15495 ***
    // Wavefunction(s) for diagram number 2687
    // (none)
    // Amplitude(s) for diagram number 2687
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[286], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[159] += amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[561] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[633] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[657] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[717] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];

    // *** DIAGRAM 2688 OF 15495 ***
    // Wavefunction(s) for diagram number 2688
    // (none)
    // Amplitude(s) for diagram number 2688
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[265], w_fp[113], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] -= amp_sv[0];
    jamp_sv[167] += amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[561] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[633] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[657] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];

    // *** DIAGRAM 2689 OF 15495 ***
    // Wavefunction(s) for diagram number 2689
    // (none)
    // Amplitude(s) for diagram number 2689
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[258], w_fp[116], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[159] += amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];
    jamp_sv[453] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[717] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];

    // *** DIAGRAM 2690 OF 15495 ***
    // Wavefunction(s) for diagram number 2690
    // (none)
    // Amplitude(s) for diagram number 2690
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[258], w_fp[113], w_fp[7], w_fp[479], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[159] += amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[561] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[633] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[657] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[717] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[258], w_fp[113], w_fp[7], w_fp[479], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] -= amp_sv[0];
    jamp_sv[167] += amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[561] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[633] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[657] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[258], w_fp[113], w_fp[7], w_fp[479], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += amp_sv[0];
    jamp_sv[159] -= amp_sv[0];
    jamp_sv[165] -= amp_sv[0];
    jamp_sv[167] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    jamp_sv[695] -= amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];

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


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup270( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 45 );
    retrieveWf( wfs, w_cx, nevt, 113 );
    retrieveWf( wfs, w_cx, nevt, 116 );
    retrieveWf( wfs, w_cx, nevt, 122 );
    retrieveWf( wfs, w_cx, nevt, 125 );
    retrieveWf( wfs, w_cx, nevt, 244 );
    retrieveWf( wfs, w_cx, nevt, 258 );
    retrieveWf( wfs, w_cx, nevt, 265 );
    retrieveWf( wfs, w_cx, nevt, 286 );
    retrieveWf( wfs, w_cx, nevt, 287 );
    retrieveWf( wfs, w_cx, nevt, 451 );
    retrieveWf( wfs, w_cx, nevt, 464 );
    retrieveWf( wfs, w_cx, nevt, 479 );
    retrieveWf( wfs, w_cx, nevt, 483 );
    retrieveWf( wfs, w_cx, nevt, 487 );
#endif
#endif

    // *** DIAGRAM 2691 OF 15495 ***
    // Wavefunction(s) for diagram number 2691
    // (none)
    // Amplitude(s) for diagram number 2691
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[464], w_fp[483], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] += amp_sv[0];
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];

    // *** DIAGRAM 2692 OF 15495 ***
    // Wavefunction(s) for diagram number 2692
    // (none)
    // Amplitude(s) for diagram number 2692
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[464], w_fp[2], w_fp[265], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[633] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[657] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2693 OF 15495 ***
    // Wavefunction(s) for diagram number 2693
    // (none)
    // Amplitude(s) for diagram number 2693
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[483], w_fp[113], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += amp_sv[0];
    jamp_sv[159] -= amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];

    // *** DIAGRAM 2694 OF 15495 ***
    // Wavefunction(s) for diagram number 2694
    // (none)
    // Amplitude(s) for diagram number 2694
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[2], w_fp[286], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[159] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2695 OF 15495 ***
    // Wavefunction(s) for diagram number 2695
    // (none)
    // Amplitude(s) for diagram number 2695
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[244], w_fp[258], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[441] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[561] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];

    // *** DIAGRAM 2696 OF 15495 ***
    // Wavefunction(s) for diagram number 2696
    // (none)
    // Amplitude(s) for diagram number 2696
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[45], w_fp[483], w_fp[116], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[159] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[165] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2697 OF 15495 ***
    // Wavefunction(s) for diagram number 2697
    // (none)
    // Amplitude(s) for diagram number 2697
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[45], w_fp[244], w_fp[265], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[441] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[453] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2698 OF 15495 ***
    // Wavefunction(s) for diagram number 2698
    // (none)
    // Amplitude(s) for diagram number 2698
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[487], w_fp[122], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[477] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];

    // *** DIAGRAM 2699 OF 15495 ***
    // Wavefunction(s) for diagram number 2699
    // (none)
    // Amplitude(s) for diagram number 2699
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[487], w_fp[2], w_fp[125], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[477] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2700 OF 15495 ***
    // Wavefunction(s) for diagram number 2700
    // (none)
    // Amplitude(s) for diagram number 2700
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[287], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[275] += amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    jamp_sv[513] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[537] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[597] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];
    jamp_sv[681] -= amp_sv[0];
    jamp_sv[687] += amp_sv[0];

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
