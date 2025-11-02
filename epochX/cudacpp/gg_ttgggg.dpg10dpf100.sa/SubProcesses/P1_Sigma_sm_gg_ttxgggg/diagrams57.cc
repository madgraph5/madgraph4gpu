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
  diagramgroup561( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 9 );
    retrieveWf( wfs, w_cx, nevt, 16 );
    retrieveWf( wfs, w_cx, nevt, 27 );
    retrieveWf( wfs, w_cx, nevt, 435 );
    retrieveWf( wfs, w_cx, nevt, 516 );
    retrieveWf( wfs, w_cx, nevt, 547 );
    retrieveWf( wfs, w_cx, nevt, 565 );
    retrieveWf( wfs, w_cx, nevt, 568 );
    retrieveWf( wfs, w_cx, nevt, 569 );
    retrieveWf( wfs, w_cx, nevt, 570 );
#endif
#endif

    // *** DIAGRAM 5601 OF 15495 ***
    // Wavefunction(s) for diagram number 5601
    // (none)
    // Amplitude(s) for diagram number 5601
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[565], w_fp[7], w_fp[435], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5602 OF 15495 ***
    // Wavefunction(s) for diagram number 5602
    // (none)
    // Amplitude(s) for diagram number 5602
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[565], w_fp[4], w_fp[516], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5603 OF 15495 ***
    // Wavefunction(s) for diagram number 5603
    // (none)
    // Amplitude(s) for diagram number 5603
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[547], w_fp[1], w_fp[27], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[61] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[547], w_fp[1], w_fp[27], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[61] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[64] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[547], w_fp[1], w_fp[27], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5604 OF 15495 ***
    // Wavefunction(s) for diagram number 5604
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[547], w_fp[1], COUPs[0], 1.0, 0., 0., w_fp[585] );
    // Amplitude(s) for diagram number 5604
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[27], w_fp[7], w_fp[585], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[61] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5605 OF 15495 ***
    // Wavefunction(s) for diagram number 5605
    // (none)
    // Amplitude(s) for diagram number 5605
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[27], w_fp[516], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5606 OF 15495 ***
    // Wavefunction(s) for diagram number 5606
    // (none)
    // Amplitude(s) for diagram number 5606
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[547], w_fp[1], w_fp[16], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[547], w_fp[1], w_fp[16], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[547], w_fp[1], w_fp[16], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5607 OF 15495 ***
    // Wavefunction(s) for diagram number 5607
    // (none)
    // Amplitude(s) for diagram number 5607
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[4], w_fp[585], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5608 OF 15495 ***
    // Wavefunction(s) for diagram number 5608
    // (none)
    // Amplitude(s) for diagram number 5608
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[16], w_fp[435], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5609 OF 15495 ***
    // Wavefunction(s) for diagram number 5609
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[4], COUPs[2], 1.0, 0., 0., w_fp[490] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[4], COUPs[2], 1.0, 0., 0., w_fp[564] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[4], COUPs[2], 1.0, 0., 0., w_fp[563] );
    // Amplitude(s) for diagram number 5609
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[547], w_fp[490], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[547], w_fp[564], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[547], w_fp[563], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5610 OF 15495 ***
    // Wavefunction(s) for diagram number 5610
    // (none)
    // Amplitude(s) for diagram number 5610
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[547], w_fp[568], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[62] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[547], w_fp[569], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[62] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[547], w_fp[570], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 563 );
    storeWf( wfs, w_cx, nevt, 564 );
    storeWf( wfs, w_cx, nevt, 585 );
#endif
#endif
  }

  //--------------------------------------------------------------------------

#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup562( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 14 );
    retrieveWf( wfs, w_cx, nevt, 27 );
    retrieveWf( wfs, w_cx, nevt, 34 );
    retrieveWf( wfs, w_cx, nevt, 35 );
    retrieveWf( wfs, w_cx, nevt, 36 );
    retrieveWf( wfs, w_cx, nevt, 488 );
    retrieveWf( wfs, w_cx, nevt, 543 );
    retrieveWf( wfs, w_cx, nevt, 547 );
    retrieveWf( wfs, w_cx, nevt, 565 );
    retrieveWf( wfs, w_cx, nevt, 578 );
#endif
#endif

    // *** DIAGRAM 5611 OF 15495 ***
    // Wavefunction(s) for diagram number 5611
    // (none)
    // Amplitude(s) for diagram number 5611
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[547], w_fp[1], w_fp[34], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[60] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[547], w_fp[1], w_fp[35], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[60] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[547], w_fp[1], w_fp[36], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5612 OF 15495 ***
    // Wavefunction(s) for diagram number 5612
    // (none)
    // Amplitude(s) for diagram number 5612
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[565], w_fp[4], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[565], w_fp[4], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[565], w_fp[4], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[71] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5613 OF 15495 ***
    // Wavefunction(s) for diagram number 5613
    // (none)
    // Amplitude(s) for diagram number 5613
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[565], w_fp[6], w_fp[578], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5614 OF 15495 ***
    // Wavefunction(s) for diagram number 5614
    // (none)
    // Amplitude(s) for diagram number 5614
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[565], w_fp[4], w_fp[543], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[71] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5615 OF 15495 ***
    // Wavefunction(s) for diagram number 5615
    // (none)
    // Amplitude(s) for diagram number 5615
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[1], w_fp[27], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[1], w_fp[27], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[1], w_fp[27], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5616 OF 15495 ***
    // Wavefunction(s) for diagram number 5616
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[488], w_fp[1], COUPs[0], 1.0, 0., 0., w_fp[562] );
    // Amplitude(s) for diagram number 5616
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[27], w_fp[6], w_fp[562], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5617 OF 15495 ***
    // Wavefunction(s) for diagram number 5617
    // (none)
    // Amplitude(s) for diagram number 5617
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[27], w_fp[543], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5618 OF 15495 ***
    // Wavefunction(s) for diagram number 5618
    // (none)
    // Amplitude(s) for diagram number 5618
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[1], w_fp[14], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[1], w_fp[14], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[1], w_fp[14], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5619 OF 15495 ***
    // Wavefunction(s) for diagram number 5619
    // (none)
    // Amplitude(s) for diagram number 5619
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[14], w_fp[4], w_fp[562], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5620 OF 15495 ***
    // Wavefunction(s) for diagram number 5620
    // (none)
    // Amplitude(s) for diagram number 5620
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[14], w_fp[578], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 562 );
#endif
#endif
  }

  //--------------------------------------------------------------------------


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup563( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 9 );
    retrieveWf( wfs, w_cx, nevt, 14 );
    retrieveWf( wfs, w_cx, nevt, 16 );
    retrieveWf( wfs, w_cx, nevt, 27 );
    retrieveWf( wfs, w_cx, nevt, 31 );
    retrieveWf( wfs, w_cx, nevt, 32 );
    retrieveWf( wfs, w_cx, nevt, 33 );
    retrieveWf( wfs, w_cx, nevt, 479 );
    retrieveWf( wfs, w_cx, nevt, 488 );
    retrieveWf( wfs, w_cx, nevt, 490 );
    retrieveWf( wfs, w_cx, nevt, 552 );
    retrieveWf( wfs, w_cx, nevt, 563 );
    retrieveWf( wfs, w_cx, nevt, 564 );
    retrieveWf( wfs, w_cx, nevt, 566 );
    retrieveWf( wfs, w_cx, nevt, 567 );
#endif
#endif

    // *** DIAGRAM 5621 OF 15495 ***
    // Wavefunction(s) for diagram number 5621
    // (none)
    // Amplitude(s) for diagram number 5621
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[490], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[564], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[563], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[71] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5622 OF 15495 ***
    // Wavefunction(s) for diagram number 5622
    // (none)
    // Amplitude(s) for diagram number 5622
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[552], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[566], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[567], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5623 OF 15495 ***
    // Wavefunction(s) for diagram number 5623
    // (none)
    // Amplitude(s) for diagram number 5623
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[1], w_fp[31], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[1], w_fp[32], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[1], w_fp[33], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[67] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5624 OF 15495 ***
    // Wavefunction(s) for diagram number 5624
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[4], COUPs[2], 1.0, 0., 0., w_fp[524] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[4], COUPs[2], 1.0, 0., 0., w_fp[523] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[4], COUPs[2], 1.0, 0., 0., w_fp[522] );
    // Amplitude(s) for diagram number 5624
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[524], w_fp[14], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[523], w_fp[14], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
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
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[522], w_fp[14], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
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
    jamp_sv[492] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5625 OF 15495 ***
    // Wavefunction(s) for diagram number 5625
    // (none)
    // Amplitude(s) for diagram number 5625
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[524], w_fp[16], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[523], w_fp[16], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[522], w_fp[16], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5626 OF 15495 ***
    // Wavefunction(s) for diagram number 5626
    // (none)
    // Amplitude(s) for diagram number 5626
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[6], w_fp[7], w_fp[524], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[6], w_fp[7], w_fp[524], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[6], w_fp[7], w_fp[524], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[6], w_fp[7], w_fp[523], COUPs[2], 1.0, &amp_fp[0] );
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
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[6], w_fp[7], w_fp[523], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[6], w_fp[7], w_fp[523], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[6], w_fp[7], w_fp[522], COUPs[2], 1.0, &amp_fp[0] );
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
    jamp_sv[492] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[6], w_fp[7], w_fp[522], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[6], w_fp[7], w_fp[522], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5627 OF 15495 ***
    // Wavefunction(s) for diagram number 5627
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[557] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[556] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[555] );
    // Amplitude(s) for diagram number 5627
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[557], w_fp[27], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[556], w_fp[27], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[555], w_fp[27], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5628 OF 15495 ***
    // Wavefunction(s) for diagram number 5628
    // (none)
    // Amplitude(s) for diagram number 5628
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[557], w_fp[16], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[556], w_fp[16], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[60] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[555], w_fp[16], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5629 OF 15495 ***
    // Wavefunction(s) for diagram number 5629
    // (none)
    // Amplitude(s) for diagram number 5629
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[4], w_fp[7], w_fp[557], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[4], w_fp[7], w_fp[557], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[4], w_fp[7], w_fp[557], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[4], w_fp[7], w_fp[556], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[4], w_fp[7], w_fp[556], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[60] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[4], w_fp[7], w_fp[556], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[60] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[4], w_fp[7], w_fp[555], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[4], w_fp[7], w_fp[555], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[4], w_fp[7], w_fp[555], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5630 OF 15495 ***
    // Wavefunction(s) for diagram number 5630
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[7], COUPs[2], 1.0, 0., 0., w_fp[560] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[7], COUPs[2], 1.0, 0., 0., w_fp[559] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[7], COUPs[2], 1.0, 0., 0., w_fp[558] );
    // Amplitude(s) for diagram number 5630
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[560], w_fp[27], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[559], w_fp[27], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[67] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[558], w_fp[27], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[53] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 522 );
    storeWf( wfs, w_cx, nevt, 523 );
    storeWf( wfs, w_cx, nevt, 524 );
    storeWf( wfs, w_cx, nevt, 555 );
    storeWf( wfs, w_cx, nevt, 556 );
    storeWf( wfs, w_cx, nevt, 557 );
    storeWf( wfs, w_cx, nevt, 558 );
    storeWf( wfs, w_cx, nevt, 559 );
    storeWf( wfs, w_cx, nevt, 560 );
#endif
#endif
  }

  //--------------------------------------------------------------------------


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup564( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 9 );
    retrieveWf( wfs, w_cx, nevt, 14 );
    retrieveWf( wfs, w_cx, nevt, 16 );
    retrieveWf( wfs, w_cx, nevt, 27 );
    retrieveWf( wfs, w_cx, nevt, 447 );
    retrieveWf( wfs, w_cx, nevt, 495 );
    retrieveWf( wfs, w_cx, nevt, 527 );
    retrieveWf( wfs, w_cx, nevt, 528 );
    retrieveWf( wfs, w_cx, nevt, 533 );
    retrieveWf( wfs, w_cx, nevt, 535 );
    retrieveWf( wfs, w_cx, nevt, 538 );
    retrieveWf( wfs, w_cx, nevt, 545 );
    retrieveWf( wfs, w_cx, nevt, 553 );
    retrieveWf( wfs, w_cx, nevt, 558 );
    retrieveWf( wfs, w_cx, nevt, 559 );
    retrieveWf( wfs, w_cx, nevt, 560 );
    retrieveWf( wfs, w_cx, nevt, 565 );
#endif
#endif

    // *** DIAGRAM 5631 OF 15495 ***
    // Wavefunction(s) for diagram number 5631
    // (none)
    // Amplitude(s) for diagram number 5631
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[560], w_fp[14], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[559], w_fp[14], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[558], w_fp[14], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[52] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5632 OF 15495 ***
    // Wavefunction(s) for diagram number 5632
    // (none)
    // Amplitude(s) for diagram number 5632
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[4], w_fp[6], w_fp[560], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[4], w_fp[6], w_fp[560], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[4], w_fp[6], w_fp[560], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[4], w_fp[6], w_fp[559], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[67] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[4], w_fp[6], w_fp[559], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[4], w_fp[6], w_fp[559], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[4], w_fp[6], w_fp[558], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[53] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[4], w_fp[6], w_fp[558], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[52] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[4], w_fp[6], w_fp[558], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[52] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5633 OF 15495 ***
    // Wavefunction(s) for diagram number 5633
    // (none)
    // Amplitude(s) for diagram number 5633
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[535], w_fp[565], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[553], w_fp[565], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[495], w_fp[565], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5634 OF 15495 ***
    // Wavefunction(s) for diagram number 5634
    // (none)
    // Amplitude(s) for diagram number 5634
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[535], w_fp[1], w_fp[16], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[553], w_fp[1], w_fp[16], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[62] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[495], w_fp[1], w_fp[16], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5635 OF 15495 ***
    // Wavefunction(s) for diagram number 5635
    // (none)
    // Amplitude(s) for diagram number 5635
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[7], w_fp[535], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[7], w_fp[535], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[7], w_fp[535], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[7], w_fp[553], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[7], w_fp[553], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[7], w_fp[553], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[7], w_fp[495], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[7], w_fp[495], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[7], w_fp[495], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5636 OF 15495 ***
    // Wavefunction(s) for diagram number 5636
    // (none)
    // Amplitude(s) for diagram number 5636
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[545], w_fp[565], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[59] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[528], w_fp[565], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[447], w_fp[565], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5637 OF 15495 ***
    // Wavefunction(s) for diagram number 5637
    // (none)
    // Amplitude(s) for diagram number 5637
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[545], w_fp[1], w_fp[14], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[528], w_fp[1], w_fp[14], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[447], w_fp[1], w_fp[14], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5638 OF 15495 ***
    // Wavefunction(s) for diagram number 5638
    // (none)
    // Amplitude(s) for diagram number 5638
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[6], w_fp[545], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[59] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[6], w_fp[545], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[6], w_fp[545], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[6], w_fp[528], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[6], w_fp[528], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[6], w_fp[528], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[6], w_fp[447], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[6], w_fp[447], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[6], w_fp[447], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5639 OF 15495 ***
    // Wavefunction(s) for diagram number 5639
    // (none)
    // Amplitude(s) for diagram number 5639
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[538], w_fp[565], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[527], w_fp[565], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[533], w_fp[565], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5640 OF 15495 ***
    // Wavefunction(s) for diagram number 5640
    // (none)
    // Amplitude(s) for diagram number 5640
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[538], w_fp[1], w_fp[27], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[527], w_fp[1], w_fp[27], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[533], w_fp[1], w_fp[27], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] -= cxtype( 0, 1 ) * amp_sv[0];

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
  diagramgroup565( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 9 );
    retrieveWf( wfs, w_cx, nevt, 16 );
    retrieveWf( wfs, w_cx, nevt, 86 );
    retrieveWf( wfs, w_cx, nevt, 88 );
    retrieveWf( wfs, w_cx, nevt, 90 );
    retrieveWf( wfs, w_cx, nevt, 362 );
    retrieveWf( wfs, w_cx, nevt, 474 );
    retrieveWf( wfs, w_cx, nevt, 521 );
    retrieveWf( wfs, w_cx, nevt, 527 );
    retrieveWf( wfs, w_cx, nevt, 533 );
    retrieveWf( wfs, w_cx, nevt, 538 );
    retrieveWf( wfs, w_cx, nevt, 565 );
    retrieveWf( wfs, w_cx, nevt, 571 );
#endif
#endif

    // *** DIAGRAM 5641 OF 15495 ***
    // Wavefunction(s) for diagram number 5641
    // (none)
    // Amplitude(s) for diagram number 5641
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[4], w_fp[538], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[4], w_fp[538], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[4], w_fp[538], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[4], w_fp[527], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[4], w_fp[527], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[4], w_fp[527], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[4], w_fp[533], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[4], w_fp[533], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[4], w_fp[533], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5642 OF 15495 ***
    // Wavefunction(s) for diagram number 5642
    // (none)
    // Amplitude(s) for diagram number 5642
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[521], w_fp[90], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5643 OF 15495 ***
    // Wavefunction(s) for diagram number 5643
    // (none)
    // Amplitude(s) for diagram number 5643
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[521], w_fp[16], w_fp[86], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5644 OF 15495 ***
    // Wavefunction(s) for diagram number 5644
    // (none)
    // Amplitude(s) for diagram number 5644
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[521], w_fp[9], w_fp[88], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5645 OF 15495 ***
    // Wavefunction(s) for diagram number 5645
    // (none)
    // Amplitude(s) for diagram number 5645
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[86], w_fp[7], w_fp[521], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[86], w_fp[7], w_fp[521], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[86], w_fp[7], w_fp[521], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5646 OF 15495 ***
    // Wavefunction(s) for diagram number 5646
    // (none)
    // Amplitude(s) for diagram number 5646
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[474], w_fp[362], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5647 OF 15495 ***
    // Wavefunction(s) for diagram number 5647
    // (none)
    // Amplitude(s) for diagram number 5647
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[474], w_fp[1], w_fp[88], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5648 OF 15495 ***
    // Wavefunction(s) for diagram number 5648
    // (none)
    // Amplitude(s) for diagram number 5648
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[86], w_fp[7], w_fp[474], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[86], w_fp[7], w_fp[474], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[52] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[86], w_fp[7], w_fp[474], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5649 OF 15495 ***
    // Wavefunction(s) for diagram number 5649
    // (none)
    // Amplitude(s) for diagram number 5649
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[571], w_fp[565], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5650 OF 15495 ***
    // Wavefunction(s) for diagram number 5650
    // (none)
    // Amplitude(s) for diagram number 5650
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[571], w_fp[1], w_fp[16], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];

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
  diagramgroup566( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 9 );
    retrieveWf( wfs, w_cx, nevt, 16 );
    retrieveWf( wfs, w_cx, nevt, 86 );
    retrieveWf( wfs, w_cx, nevt, 88 );
    retrieveWf( wfs, w_cx, nevt, 90 );
    retrieveWf( wfs, w_cx, nevt, 362 );
    retrieveWf( wfs, w_cx, nevt, 479 );
    retrieveWf( wfs, w_cx, nevt, 488 );
    retrieveWf( wfs, w_cx, nevt, 558 );
    retrieveWf( wfs, w_cx, nevt, 559 );
    retrieveWf( wfs, w_cx, nevt, 560 );
    retrieveWf( wfs, w_cx, nevt, 565 );
    retrieveWf( wfs, w_cx, nevt, 571 );
#endif
#endif

    // *** DIAGRAM 5651 OF 15495 ***
    // Wavefunction(s) for diagram number 5651
    // (none)
    // Amplitude(s) for diagram number 5651
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[7], w_fp[571], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[7], w_fp[571], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[7], w_fp[571], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5652 OF 15495 ***
    // Wavefunction(s) for diagram number 5652
    // (none)
    // Amplitude(s) for diagram number 5652
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[565], w_fp[86], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5653 OF 15495 ***
    // Wavefunction(s) for diagram number 5653
    // (none)
    // Amplitude(s) for diagram number 5653
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[362], w_fp[9], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5654 OF 15495 ***
    // Wavefunction(s) for diagram number 5654
    // (none)
    // Amplitude(s) for diagram number 5654
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[1], w_fp[90], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5655 OF 15495 ***
    // Wavefunction(s) for diagram number 5655
    // (none)
    // Amplitude(s) for diagram number 5655
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[86], w_fp[488], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[86], w_fp[488], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[86], w_fp[488], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5656 OF 15495 ***
    // Wavefunction(s) for diagram number 5656
    // (none)
    // Amplitude(s) for diagram number 5656
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[565], w_fp[88], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5657 OF 15495 ***
    // Wavefunction(s) for diagram number 5657
    // (none)
    // Amplitude(s) for diagram number 5657
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[362], w_fp[16], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5658 OF 15495 ***
    // Wavefunction(s) for diagram number 5658
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[9], COUPs[2], 1.0, 0., 0., w_fp[591] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[9], COUPs[2], 1.0, 0., 0., w_fp[592] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[9], COUPs[2], 1.0, 0., 0., w_fp[593] );
    // Amplitude(s) for diagram number 5658
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[86], w_fp[7], w_fp[591], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[86], w_fp[7], w_fp[592], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[86], w_fp[7], w_fp[593], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5659 OF 15495 ***
    // Wavefunction(s) for diagram number 5659
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[86], COUPs[2], 1.0, 0., 0., w_fp[594] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[86], COUPs[2], 1.0, 0., 0., w_fp[595] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[86], COUPs[2], 1.0, 0., 0., w_fp[596] );
    // Amplitude(s) for diagram number 5659
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[7], w_fp[594], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[7], w_fp[595], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[7], w_fp[596], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5660 OF 15495 ***
    // Wavefunction(s) for diagram number 5660
    // (none)
    // Amplitude(s) for diagram number 5660
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[86], w_fp[560], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[52] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[86], w_fp[559], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[86], w_fp[558], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 591 );
    storeWf( wfs, w_cx, nevt, 592 );
    storeWf( wfs, w_cx, nevt, 593 );
    storeWf( wfs, w_cx, nevt, 594 );
    storeWf( wfs, w_cx, nevt, 595 );
    storeWf( wfs, w_cx, nevt, 596 );
#endif
#endif
  }

  //--------------------------------------------------------------------------


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup567( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 9 );
    retrieveWf( wfs, w_cx, nevt, 14 );
    retrieveWf( wfs, w_cx, nevt, 86 );
    retrieveWf( wfs, w_cx, nevt, 102 );
    retrieveWf( wfs, w_cx, nevt, 104 );
    retrieveWf( wfs, w_cx, nevt, 106 );
    retrieveWf( wfs, w_cx, nevt, 363 );
    retrieveWf( wfs, w_cx, nevt, 444 );
    retrieveWf( wfs, w_cx, nevt, 450 );
    retrieveWf( wfs, w_cx, nevt, 452 );
    retrieveWf( wfs, w_cx, nevt, 474 );
    retrieveWf( wfs, w_cx, nevt, 505 );
    retrieveWf( wfs, w_cx, nevt, 509 );
    retrieveWf( wfs, w_cx, nevt, 521 );
    retrieveWf( wfs, w_cx, nevt, 531 );
    retrieveWf( wfs, w_cx, nevt, 565 );
    retrieveWf( wfs, w_cx, nevt, 579 );
#endif
#endif

    // *** DIAGRAM 5661 OF 15495 ***
    // Wavefunction(s) for diagram number 5661
    // (none)
    // Amplitude(s) for diagram number 5661
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[86], w_fp[452], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[86], w_fp[444], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[86], w_fp[509], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5662 OF 15495 ***
    // Wavefunction(s) for diagram number 5662
    // (none)
    // Amplitude(s) for diagram number 5662
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[450], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[505], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[531], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5663 OF 15495 ***
    // Wavefunction(s) for diagram number 5663
    // (none)
    // Amplitude(s) for diagram number 5663
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[521], w_fp[106], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5664 OF 15495 ***
    // Wavefunction(s) for diagram number 5664
    // (none)
    // Amplitude(s) for diagram number 5664
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[521], w_fp[14], w_fp[102], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5665 OF 15495 ***
    // Wavefunction(s) for diagram number 5665
    // (none)
    // Amplitude(s) for diagram number 5665
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[521], w_fp[9], w_fp[104], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5666 OF 15495 ***
    // Wavefunction(s) for diagram number 5666
    // (none)
    // Amplitude(s) for diagram number 5666
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[102], w_fp[6], w_fp[521], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[102], w_fp[6], w_fp[521], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[102], w_fp[6], w_fp[521], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5667 OF 15495 ***
    // Wavefunction(s) for diagram number 5667
    // (none)
    // Amplitude(s) for diagram number 5667
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[474], w_fp[363], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5668 OF 15495 ***
    // Wavefunction(s) for diagram number 5668
    // (none)
    // Amplitude(s) for diagram number 5668
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[474], w_fp[1], w_fp[104], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5669 OF 15495 ***
    // Wavefunction(s) for diagram number 5669
    // (none)
    // Amplitude(s) for diagram number 5669
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[102], w_fp[6], w_fp[474], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[102], w_fp[6], w_fp[474], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[102], w_fp[6], w_fp[474], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5670 OF 15495 ***
    // Wavefunction(s) for diagram number 5670
    // (none)
    // Amplitude(s) for diagram number 5670
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[579], w_fp[565], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[59] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];

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
  diagramgroup568( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 9 );
    retrieveWf( wfs, w_cx, nevt, 14 );
    retrieveWf( wfs, w_cx, nevt, 102 );
    retrieveWf( wfs, w_cx, nevt, 104 );
    retrieveWf( wfs, w_cx, nevt, 106 );
    retrieveWf( wfs, w_cx, nevt, 363 );
    retrieveWf( wfs, w_cx, nevt, 479 );
    retrieveWf( wfs, w_cx, nevt, 547 );
    retrieveWf( wfs, w_cx, nevt, 565 );
    retrieveWf( wfs, w_cx, nevt, 579 );
    retrieveWf( wfs, w_cx, nevt, 591 );
    retrieveWf( wfs, w_cx, nevt, 592 );
    retrieveWf( wfs, w_cx, nevt, 593 );
#endif
#endif

    // *** DIAGRAM 5671 OF 15495 ***
    // Wavefunction(s) for diagram number 5671
    // (none)
    // Amplitude(s) for diagram number 5671
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[579], w_fp[1], w_fp[14], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5672 OF 15495 ***
    // Wavefunction(s) for diagram number 5672
    // (none)
    // Amplitude(s) for diagram number 5672
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[6], w_fp[579], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[59] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[6], w_fp[579], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[6], w_fp[579], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5673 OF 15495 ***
    // Wavefunction(s) for diagram number 5673
    // (none)
    // Amplitude(s) for diagram number 5673
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[547], w_fp[565], w_fp[102], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5674 OF 15495 ***
    // Wavefunction(s) for diagram number 5674
    // (none)
    // Amplitude(s) for diagram number 5674
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[547], w_fp[363], w_fp[9], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5675 OF 15495 ***
    // Wavefunction(s) for diagram number 5675
    // (none)
    // Amplitude(s) for diagram number 5675
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[547], w_fp[1], w_fp[106], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[60] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5676 OF 15495 ***
    // Wavefunction(s) for diagram number 5676
    // (none)
    // Amplitude(s) for diagram number 5676
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[102], w_fp[547], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[102], w_fp[547], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[102], w_fp[547], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5677 OF 15495 ***
    // Wavefunction(s) for diagram number 5677
    // (none)
    // Amplitude(s) for diagram number 5677
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[565], w_fp[104], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[59] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5678 OF 15495 ***
    // Wavefunction(s) for diagram number 5678
    // (none)
    // Amplitude(s) for diagram number 5678
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[363], w_fp[14], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5679 OF 15495 ***
    // Wavefunction(s) for diagram number 5679
    // (none)
    // Amplitude(s) for diagram number 5679
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[102], w_fp[6], w_fp[591], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[59] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[102], w_fp[6], w_fp[592], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[102], w_fp[6], w_fp[593], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5680 OF 15495 ***
    // Wavefunction(s) for diagram number 5680
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[102], COUPs[2], 1.0, 0., 0., w_fp[531] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[102], COUPs[2], 1.0, 0., 0., w_fp[505] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[102], COUPs[2], 1.0, 0., 0., w_fp[450] );
    // Amplitude(s) for diagram number 5680
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[6], w_fp[531], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[6], w_fp[505], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[6], w_fp[450], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 450 );
    storeWf( wfs, w_cx, nevt, 505 );
    storeWf( wfs, w_cx, nevt, 531 );
#endif
#endif
  }

  //--------------------------------------------------------------------------


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup569( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 9 );
    retrieveWf( wfs, w_cx, nevt, 27 );
    retrieveWf( wfs, w_cx, nevt, 84 );
    retrieveWf( wfs, w_cx, nevt, 102 );
    retrieveWf( wfs, w_cx, nevt, 130 );
    retrieveWf( wfs, w_cx, nevt, 132 );
    retrieveWf( wfs, w_cx, nevt, 361 );
    retrieveWf( wfs, w_cx, nevt, 445 );
    retrieveWf( wfs, w_cx, nevt, 472 );
    retrieveWf( wfs, w_cx, nevt, 474 );
    retrieveWf( wfs, w_cx, nevt, 521 );
    retrieveWf( wfs, w_cx, nevt, 555 );
    retrieveWf( wfs, w_cx, nevt, 556 );
    retrieveWf( wfs, w_cx, nevt, 557 );
    retrieveWf( wfs, w_cx, nevt, 561 );
    retrieveWf( wfs, w_cx, nevt, 583 );
    retrieveWf( wfs, w_cx, nevt, 584 );
    retrieveWf( wfs, w_cx, nevt, 587 );
#endif
#endif

    // *** DIAGRAM 5681 OF 15495 ***
    // Wavefunction(s) for diagram number 5681
    // (none)
    // Amplitude(s) for diagram number 5681
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[102], w_fp[557], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[102], w_fp[556], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[102], w_fp[555], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5682 OF 15495 ***
    // Wavefunction(s) for diagram number 5682
    // (none)
    // Amplitude(s) for diagram number 5682
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[102], w_fp[587], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[102], w_fp[472], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[102], w_fp[445], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5683 OF 15495 ***
    // Wavefunction(s) for diagram number 5683
    // (none)
    // Amplitude(s) for diagram number 5683
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[561], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[584], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[583], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[59] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5684 OF 15495 ***
    // Wavefunction(s) for diagram number 5684
    // (none)
    // Amplitude(s) for diagram number 5684
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[521], w_fp[27], w_fp[84], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5685 OF 15495 ***
    // Wavefunction(s) for diagram number 5685
    // (none)
    // Amplitude(s) for diagram number 5685
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[521], w_fp[132], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5686 OF 15495 ***
    // Wavefunction(s) for diagram number 5686
    // (none)
    // Amplitude(s) for diagram number 5686
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[521], w_fp[9], w_fp[130], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5687 OF 15495 ***
    // Wavefunction(s) for diagram number 5687
    // (none)
    // Amplitude(s) for diagram number 5687
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[4], w_fp[84], w_fp[521], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[4], w_fp[84], w_fp[521], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[4], w_fp[84], w_fp[521], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5688 OF 15495 ***
    // Wavefunction(s) for diagram number 5688
    // (none)
    // Amplitude(s) for diagram number 5688
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[474], w_fp[361], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5689 OF 15495 ***
    // Wavefunction(s) for diagram number 5689
    // (none)
    // Amplitude(s) for diagram number 5689
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[474], w_fp[1], w_fp[130], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5690 OF 15495 ***
    // Wavefunction(s) for diagram number 5690
    // (none)
    // Amplitude(s) for diagram number 5690
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[4], w_fp[84], w_fp[474], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[4], w_fp[84], w_fp[474], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[4], w_fp[84], w_fp[474], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] -= cxtype( 0, 1 ) * amp_sv[0];

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
  diagramgroup570( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 9 );
    retrieveWf( wfs, w_cx, nevt, 27 );
    retrieveWf( wfs, w_cx, nevt, 84 );
    retrieveWf( wfs, w_cx, nevt, 130 );
    retrieveWf( wfs, w_cx, nevt, 132 );
    retrieveWf( wfs, w_cx, nevt, 361 );
    retrieveWf( wfs, w_cx, nevt, 479 );
    retrieveWf( wfs, w_cx, nevt, 565 );
    retrieveWf( wfs, w_cx, nevt, 577 );
    retrieveWf( wfs, w_cx, nevt, 581 );
    retrieveWf( wfs, w_cx, nevt, 591 );
    retrieveWf( wfs, w_cx, nevt, 592 );
    retrieveWf( wfs, w_cx, nevt, 593 );
#endif
#endif

    // *** DIAGRAM 5691 OF 15495 ***
    // Wavefunction(s) for diagram number 5691
    // (none)
    // Amplitude(s) for diagram number 5691
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[577], w_fp[565], w_fp[84], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5692 OF 15495 ***
    // Wavefunction(s) for diagram number 5692
    // (none)
    // Amplitude(s) for diagram number 5692
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[577], w_fp[361], w_fp[9], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5693 OF 15495 ***
    // Wavefunction(s) for diagram number 5693
    // (none)
    // Amplitude(s) for diagram number 5693
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[577], w_fp[1], w_fp[132], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5694 OF 15495 ***
    // Wavefunction(s) for diagram number 5694
    // (none)
    // Amplitude(s) for diagram number 5694
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[84], w_fp[577], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[84], w_fp[577], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[84], w_fp[577], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5695 OF 15495 ***
    // Wavefunction(s) for diagram number 5695
    // (none)
    // Amplitude(s) for diagram number 5695
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[581], w_fp[565], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5696 OF 15495 ***
    // Wavefunction(s) for diagram number 5696
    // (none)
    // Amplitude(s) for diagram number 5696
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[581], w_fp[1], w_fp[27], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5697 OF 15495 ***
    // Wavefunction(s) for diagram number 5697
    // (none)
    // Amplitude(s) for diagram number 5697
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[4], w_fp[581], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[4], w_fp[581], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], w_fp[4], w_fp[581], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5698 OF 15495 ***
    // Wavefunction(s) for diagram number 5698
    // (none)
    // Amplitude(s) for diagram number 5698
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[565], w_fp[130], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5699 OF 15495 ***
    // Wavefunction(s) for diagram number 5699
    // (none)
    // Amplitude(s) for diagram number 5699
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[361], w_fp[27], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5700 OF 15495 ***
    // Wavefunction(s) for diagram number 5700
    // (none)
    // Amplitude(s) for diagram number 5700
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[4], w_fp[84], w_fp[591], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[4], w_fp[84], w_fp[592], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[4], w_fp[84], w_fp[593], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];

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
