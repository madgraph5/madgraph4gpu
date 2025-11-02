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
  diagramgroup701( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 45 );
    retrieveWf( wfs, w_cx, nevt, 100 );
    retrieveWf( wfs, w_cx, nevt, 156 );
    retrieveWf( wfs, w_cx, nevt, 158 );
    retrieveWf( wfs, w_cx, nevt, 164 );
    retrieveWf( wfs, w_cx, nevt, 179 );
    retrieveWf( wfs, w_cx, nevt, 180 );
    retrieveWf( wfs, w_cx, nevt, 183 );
    retrieveWf( wfs, w_cx, nevt, 186 );
    retrieveWf( wfs, w_cx, nevt, 187 );
    retrieveWf( wfs, w_cx, nevt, 252 );
    retrieveWf( wfs, w_cx, nevt, 450 );
    retrieveWf( wfs, w_cx, nevt, 474 );
    retrieveWf( wfs, w_cx, nevt, 512 );
    retrieveWf( wfs, w_cx, nevt, 525 );
    retrieveWf( wfs, w_cx, nevt, 526 );
    retrieveWf( wfs, w_cx, nevt, 530 );
    retrieveWf( wfs, w_cx, nevt, 557 );
#endif
#endif

    // *** DIAGRAM 7001 OF 15495 ***
    // Wavefunction(s) for diagram number 7001
    // (none)
    // Amplitude(s) for diagram number 7001
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[557], w_fp[158], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[294] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7002 OF 15495 ***
    // Wavefunction(s) for diagram number 7002
    // (none)
    // Amplitude(s) for diagram number 7002
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[179], w_fp[512], w_fp[474], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[266] += amp_sv[0];
    jamp_sv[270] -= amp_sv[0];
    jamp_sv[272] += amp_sv[0];
    jamp_sv[276] -= amp_sv[0];

    // *** DIAGRAM 7003 OF 15495 ***
    // Wavefunction(s) for diagram number 7003
    // (none)
    // Amplitude(s) for diagram number 7003
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[252], w_fp[156], w_fp[474], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[254] += amp_sv[0];
    jamp_sv[290] -= amp_sv[0];
    jamp_sv[300] += amp_sv[0];
    jamp_sv[314] -= amp_sv[0];

    // *** DIAGRAM 7004 OF 15495 ***
    // Wavefunction(s) for diagram number 7004
    // (none)
    // Amplitude(s) for diagram number 7004
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[474], w_fp[1], w_fp[183], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[254] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[266] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[300] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7005 OF 15495 ***
    // Wavefunction(s) for diagram number 7005
    // (none)
    // Amplitude(s) for diagram number 7005
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[180], w_fp[512], w_fp[530], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[266] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7006 OF 15495 ***
    // Wavefunction(s) for diagram number 7006
    // (none)
    // Amplitude(s) for diagram number 7006
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[252], w_fp[158], w_fp[530], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[290] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[300] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7007 OF 15495 ***
    // Wavefunction(s) for diagram number 7007
    // (none)
    // Amplitude(s) for diagram number 7007
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[179], w_fp[156], w_fp[525], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[252] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[254] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[294] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[312] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[179], w_fp[156], w_fp[526], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[254] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[266] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[300] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[179], w_fp[156], w_fp[45], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[252] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[266] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[294] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[300] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[312] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7008 OF 15495 ***
    // Wavefunction(s) for diagram number 7008
    // (none)
    // Amplitude(s) for diagram number 7008
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[450], w_fp[164], w_fp[100], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[252] += amp_sv[0];
    jamp_sv[253] -= amp_sv[0];
    jamp_sv[266] -= amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[276] += amp_sv[0];
    jamp_sv[277] -= amp_sv[0];
    jamp_sv[307] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[312] -= amp_sv[0];
    jamp_sv[313] += amp_sv[0];
    jamp_sv[349] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];
    jamp_sv[351] += amp_sv[0];
    jamp_sv[352] -= amp_sv[0];

    // *** DIAGRAM 7009 OF 15495 ***
    // Wavefunction(s) for diagram number 7009
    // (none)
    // Amplitude(s) for diagram number 7009
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[187], w_fp[450], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[307] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[309] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[310] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[350] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[351] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7010 OF 15495 ***
    // Wavefunction(s) for diagram number 7010
    // (none)
    // Amplitude(s) for diagram number 7010
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[186], w_fp[156], w_fp[450], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[252] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[253] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[266] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[267] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[277] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[312] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[313] -= cxtype( 0, 1 ) * amp_sv[0];

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
  diagramgroup702( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 100 );
    retrieveWf( wfs, w_cx, nevt, 156 );
    retrieveWf( wfs, w_cx, nevt, 164 );
    retrieveWf( wfs, w_cx, nevt, 186 );
    retrieveWf( wfs, w_cx, nevt, 187 );
    retrieveWf( wfs, w_cx, nevt, 360 );
    retrieveWf( wfs, w_cx, nevt, 445 );
    retrieveWf( wfs, w_cx, nevt, 472 );
    retrieveWf( wfs, w_cx, nevt, 510 );
    retrieveWf( wfs, w_cx, nevt, 512 );
    retrieveWf( wfs, w_cx, nevt, 530 );
    retrieveWf( wfs, w_cx, nevt, 532 );
    retrieveWf( wfs, w_cx, nevt, 553 );
    retrieveWf( wfs, w_cx, nevt, 587 );
#endif
#endif

    // *** DIAGRAM 7011 OF 15495 ***
    // Wavefunction(s) for diagram number 7011
    // (none)
    // Amplitude(s) for diagram number 7011
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[553], w_fp[360], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[252] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[253] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[312] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[313] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7012 OF 15495 ***
    // Wavefunction(s) for diagram number 7012
    // (none)
    // Amplitude(s) for diagram number 7012
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[186], w_fp[553], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[252] += amp_sv[0];
    jamp_sv[253] -= amp_sv[0];
    jamp_sv[312] -= amp_sv[0];
    jamp_sv[313] += amp_sv[0];

    // *** DIAGRAM 7013 OF 15495 ***
    // Wavefunction(s) for diagram number 7013
    // (none)
    // Amplitude(s) for diagram number 7013
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[532], w_fp[512], w_fp[100], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[274] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[284] -= amp_sv[0];
    jamp_sv[285] += amp_sv[0];

    // *** DIAGRAM 7014 OF 15495 ***
    // Wavefunction(s) for diagram number 7014
    // (none)
    // Amplitude(s) for diagram number 7014
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[532], w_fp[156], w_fp[360], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[274] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[309] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[350] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[351] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7015 OF 15495 ***
    // Wavefunction(s) for diagram number 7015
    // (none)
    // Amplitude(s) for diagram number 7015
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[532], w_fp[187], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[308] += amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[350] -= amp_sv[0];
    jamp_sv[351] += amp_sv[0];

    // *** DIAGRAM 7016 OF 15495 ***
    // Wavefunction(s) for diagram number 7016
    // (none)
    // Amplitude(s) for diagram number 7016
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[512], w_fp[510], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[266] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[267] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[277] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7017 OF 15495 ***
    // Wavefunction(s) for diagram number 7017
    // (none)
    // Amplitude(s) for diagram number 7017
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[510], w_fp[1], w_fp[164], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[266] += amp_sv[0];
    jamp_sv[267] -= amp_sv[0];
    jamp_sv[274] -= amp_sv[0];
    jamp_sv[275] += amp_sv[0];
    jamp_sv[276] -= amp_sv[0];
    jamp_sv[277] += amp_sv[0];
    jamp_sv[284] += amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[307] += amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[349] -= amp_sv[0];
    jamp_sv[352] += amp_sv[0];

    // *** DIAGRAM 7018 OF 15495 ***
    // Wavefunction(s) for diagram number 7018
    // (none)
    // Amplitude(s) for diagram number 7018
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[186], w_fp[512], w_fp[530], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[266] += amp_sv[0];
    jamp_sv[267] -= amp_sv[0];
    jamp_sv[276] -= amp_sv[0];
    jamp_sv[277] += amp_sv[0];

    // *** DIAGRAM 7019 OF 15495 ***
    // Wavefunction(s) for diagram number 7019
    // (none)
    // Amplitude(s) for diagram number 7019
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[530], w_fp[360], w_fp[164], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[252] -= amp_sv[0];
    jamp_sv[253] += amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[284] -= amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[309] += amp_sv[0];
    jamp_sv[312] += amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];

    // *** DIAGRAM 7020 OF 15495 ***
    // Wavefunction(s) for diagram number 7020
    // (none)
    // Amplitude(s) for diagram number 7020
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[156], w_fp[587], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[252] -= amp_sv[0];
    jamp_sv[253] += amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[284] -= amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[309] += amp_sv[0];
    jamp_sv[312] += amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[156], w_fp[472], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[255] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[266] -= amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[276] += amp_sv[0];
    jamp_sv[277] -= amp_sv[0];
    jamp_sv[284] -= amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[307] -= amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[349] += amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[156], w_fp[445], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[252] += amp_sv[0];
    jamp_sv[253] -= amp_sv[0];
    jamp_sv[266] -= amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[276] += amp_sv[0];
    jamp_sv[277] -= amp_sv[0];
    jamp_sv[307] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[312] -= amp_sv[0];
    jamp_sv[313] += amp_sv[0];
    jamp_sv[349] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];
    jamp_sv[351] += amp_sv[0];
    jamp_sv[352] -= amp_sv[0];

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
  diagramgroup703( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 190 );
    retrieveWf( wfs, w_cx, nevt, 193 );
    retrieveWf( wfs, w_cx, nevt, 194 );
    retrieveWf( wfs, w_cx, nevt, 437 );
    retrieveWf( wfs, w_cx, nevt, 450 );
    retrieveWf( wfs, w_cx, nevt, 453 );
    retrieveWf( wfs, w_cx, nevt, 477 );
    retrieveWf( wfs, w_cx, nevt, 524 );
    retrieveWf( wfs, w_cx, nevt, 531 );
    retrieveWf( wfs, w_cx, nevt, 562 );
    retrieveWf( wfs, w_cx, nevt, 579 );
#endif
#endif

    // *** DIAGRAM 7021 OF 15495 ***
    // Wavefunction(s) for diagram number 7021
    // (none)
    // Amplitude(s) for diagram number 7021
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[450], w_fp[194], w_fp[4], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[372] += amp_sv[0];
    jamp_sv[373] -= amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[396] += amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[427] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[429] -= amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[432] -= amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[469] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[471] += amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[450], w_fp[194], w_fp[4], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[372] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[396] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[414] += amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[432] -= amp_sv[0];
    jamp_sv[460] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[469] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[471] += amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[450], w_fp[194], w_fp[4], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[373] += amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[397] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[414] += amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[427] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[429] += amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[460] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[474] += amp_sv[0];

    // *** DIAGRAM 7022 OF 15495 ***
    // Wavefunction(s) for diagram number 7022
    // (none)
    // Amplitude(s) for diagram number 7022
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[194], w_fp[7], w_fp[531], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[372] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[396] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[414] += amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[432] -= amp_sv[0];
    jamp_sv[460] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[469] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[471] += amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[474] += amp_sv[0];

    // *** DIAGRAM 7023 OF 15495 ***
    // Wavefunction(s) for diagram number 7023
    // (none)
    // Amplitude(s) for diagram number 7023
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[194], w_fp[4], w_fp[579], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[373] += amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[397] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[414] += amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[427] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[429] += amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[460] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[474] += amp_sv[0];

    // *** DIAGRAM 7024 OF 15495 ***
    // Wavefunction(s) for diagram number 7024
    // (none)
    // Amplitude(s) for diagram number 7024
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[437], w_fp[190], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[427] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[429] += amp_sv[0];
    jamp_sv[430] -= amp_sv[0];

    // *** DIAGRAM 7025 OF 15495 ***
    // Wavefunction(s) for diagram number 7025
    // (none)
    // Amplitude(s) for diagram number 7025
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[190], w_fp[579], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[427] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[429] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7026 OF 15495 ***
    // Wavefunction(s) for diagram number 7026
    // (none)
    // Amplitude(s) for diagram number 7026
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[437], w_fp[193], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[469] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[471] += amp_sv[0];
    jamp_sv[472] -= amp_sv[0];

    // *** DIAGRAM 7027 OF 15495 ***
    // Wavefunction(s) for diagram number 7027
    // (none)
    // Amplitude(s) for diagram number 7027
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[193], w_fp[531], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[460] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[469] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[471] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7028 OF 15495 ***
    // Wavefunction(s) for diagram number 7028
    // (none)
    // Amplitude(s) for diagram number 7028
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[524], w_fp[477], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[404] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7029 OF 15495 ***
    // Wavefunction(s) for diagram number 7029
    // (none)
    // Amplitude(s) for diagram number 7029
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[453], w_fp[477], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[394] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7030 OF 15495 ***
    // Wavefunction(s) for diagram number 7030
    // (none)
    // Amplitude(s) for diagram number 7030
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[562], w_fp[190], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[428] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[429] -= cxtype( 0, 1 ) * amp_sv[0];

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
  diagramgroup704( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 190 );
    retrieveWf( wfs, w_cx, nevt, 193 );
    retrieveWf( wfs, w_cx, nevt, 194 );
    retrieveWf( wfs, w_cx, nevt, 453 );
    retrieveWf( wfs, w_cx, nevt, 477 );
    retrieveWf( wfs, w_cx, nevt, 515 );
    retrieveWf( wfs, w_cx, nevt, 523 );
    retrieveWf( wfs, w_cx, nevt, 524 );
    retrieveWf( wfs, w_cx, nevt, 562 );
    retrieveWf( wfs, w_cx, nevt, 577 );
    retrieveWf( wfs, w_cx, nevt, 578 );
#endif
#endif

    // *** DIAGRAM 7031 OF 15495 ***
    // Wavefunction(s) for diagram number 7031
    // (none)
    // Amplitude(s) for diagram number 7031
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[453], w_fp[190], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[418] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[419] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7032 OF 15495 ***
    // Wavefunction(s) for diagram number 7032
    // (none)
    // Amplitude(s) for diagram number 7032
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[562], w_fp[193], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[470] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[471] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7033 OF 15495 ***
    // Wavefunction(s) for diagram number 7033
    // (none)
    // Amplitude(s) for diagram number 7033
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[524], w_fp[193], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[464] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7034 OF 15495 ***
    // Wavefunction(s) for diagram number 7034
    // (none)
    // Amplitude(s) for diagram number 7034
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[577], w_fp[477], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[403] += amp_sv[0];
    jamp_sv[404] -= amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[406] -= amp_sv[0];

    // *** DIAGRAM 7035 OF 15495 ***
    // Wavefunction(s) for diagram number 7035
    // (none)
    // Amplitude(s) for diagram number 7035
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[477], w_fp[515], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[386] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[390] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7036 OF 15495 ***
    // Wavefunction(s) for diagram number 7036
    // (none)
    // Amplitude(s) for diagram number 7036
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[523], w_fp[1], w_fp[194], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[374] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[390] += amp_sv[0];
    jamp_sv[392] -= amp_sv[0];
    jamp_sv[396] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[434] -= amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[469] += amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[475] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[523], w_fp[1], w_fp[194], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[374] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[434] -= amp_sv[0];
    jamp_sv[435] += amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[523], w_fp[1], w_fp[194], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[396] -= amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[435] += amp_sv[0];
    jamp_sv[461] += amp_sv[0];
    jamp_sv[469] -= amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[475] -= amp_sv[0];

    // *** DIAGRAM 7037 OF 15495 ***
    // Wavefunction(s) for diagram number 7037
    // (none)
    // Amplitude(s) for diagram number 7037
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[194], w_fp[7], w_fp[578], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[374] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[390] += amp_sv[0];
    jamp_sv[392] -= amp_sv[0];
    jamp_sv[396] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[434] -= amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[469] += amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[475] += amp_sv[0];

    // *** DIAGRAM 7038 OF 15495 ***
    // Wavefunction(s) for diagram number 7038
    // (none)
    // Amplitude(s) for diagram number 7038
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[194], w_fp[515], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[396] -= amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[435] += amp_sv[0];
    jamp_sv[461] += amp_sv[0];
    jamp_sv[469] -= amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[475] -= amp_sv[0];

    // *** DIAGRAM 7039 OF 15495 ***
    // Wavefunction(s) for diagram number 7039
    // (none)
    // Amplitude(s) for diagram number 7039
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[193], w_fp[578], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[461] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[464] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[469] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[475] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7040 OF 15495 ***
    // Wavefunction(s) for diagram number 7040
    // (none)
    // Amplitude(s) for diagram number 7040
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[577], w_fp[193], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[463] += amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[466] -= amp_sv[0];

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
  diagramgroup705( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 190 );
    retrieveWf( wfs, w_cx, nevt, 193 );
    retrieveWf( wfs, w_cx, nevt, 194 );
    retrieveWf( wfs, w_cx, nevt, 477 );
    retrieveWf( wfs, w_cx, nevt, 488 );
    retrieveWf( wfs, w_cx, nevt, 516 );
    retrieveWf( wfs, w_cx, nevt, 522 );
    retrieveWf( wfs, w_cx, nevt, 528 );
    retrieveWf( wfs, w_cx, nevt, 537 );
    retrieveWf( wfs, w_cx, nevt, 547 );
    retrieveWf( wfs, w_cx, nevt, 561 );
    retrieveWf( wfs, w_cx, nevt, 582 );
    retrieveWf( wfs, w_cx, nevt, 583 );
    retrieveWf( wfs, w_cx, nevt, 584 );
#endif
#endif

    // *** DIAGRAM 7041 OF 15495 ***
    // Wavefunction(s) for diagram number 7041
    // (none)
    // Amplitude(s) for diagram number 7041
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[477], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[390] += amp_sv[0];
    jamp_sv[392] -= amp_sv[0];
    jamp_sv[394] -= amp_sv[0];
    jamp_sv[395] += amp_sv[0];

    // *** DIAGRAM 7042 OF 15495 ***
    // Wavefunction(s) for diagram number 7042
    // (none)
    // Amplitude(s) for diagram number 7042
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[477], w_fp[522], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[387] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[390] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7043 OF 15495 ***
    // Wavefunction(s) for diagram number 7043
    // (none)
    // Amplitude(s) for diagram number 7043
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[547], w_fp[1], w_fp[194], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[376] += amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[397] += amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[414] += amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[427] += amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[460] -= amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[547], w_fp[1], w_fp[194], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[376] += amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[414] += amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[460] -= amp_sv[0];
    jamp_sv[461] += amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[475] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[547], w_fp[1], w_fp[194], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[427] -= amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[461] += amp_sv[0];
    jamp_sv[475] -= amp_sv[0];

    // *** DIAGRAM 7044 OF 15495 ***
    // Wavefunction(s) for diagram number 7044
    // (none)
    // Amplitude(s) for diagram number 7044
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[194], w_fp[4], w_fp[528], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[376] += amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[397] += amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[414] += amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[427] += amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[460] -= amp_sv[0];
    jamp_sv[474] += amp_sv[0];

    // *** DIAGRAM 7045 OF 15495 ***
    // Wavefunction(s) for diagram number 7045
    // (none)
    // Amplitude(s) for diagram number 7045
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[194], w_fp[522], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[427] -= amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[461] += amp_sv[0];
    jamp_sv[475] -= amp_sv[0];

    // *** DIAGRAM 7046 OF 15495 ***
    // Wavefunction(s) for diagram number 7046
    // (none)
    // Amplitude(s) for diagram number 7046
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[190], w_fp[528], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[419] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[427] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7047 OF 15495 ***
    // Wavefunction(s) for diagram number 7047
    // (none)
    // Amplitude(s) for diagram number 7047
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[190], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[414] += amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];

    // *** DIAGRAM 7048 OF 15495 ***
    // Wavefunction(s) for diagram number 7048
    // (none)
    // Amplitude(s) for diagram number 7048
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[516], w_fp[194], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[372] += amp_sv[0];
    jamp_sv[374] -= amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[414] += amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[432] -= amp_sv[0];
    jamp_sv[434] += amp_sv[0];
    jamp_sv[460] -= amp_sv[0];
    jamp_sv[461] += amp_sv[0];
    jamp_sv[464] += amp_sv[0];
    jamp_sv[465] -= amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[471] += amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[475] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[537], w_fp[194], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[374] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[396] -= amp_sv[0];
    jamp_sv[410] += amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[434] += amp_sv[0];
    jamp_sv[461] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[464] += amp_sv[0];
    jamp_sv[465] -= amp_sv[0];
    jamp_sv[466] += amp_sv[0];
    jamp_sv[469] -= amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[475] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[582], w_fp[194], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[372] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[396] -= amp_sv[0];
    jamp_sv[410] += amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[416] += amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[432] += amp_sv[0];
    jamp_sv[460] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[466] += amp_sv[0];
    jamp_sv[469] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[474] -= amp_sv[0];

    // *** DIAGRAM 7049 OF 15495 ***
    // Wavefunction(s) for diagram number 7049
    // (none)
    // Amplitude(s) for diagram number 7049
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[193], w_fp[516], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[460] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[464] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[471] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[475] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[193], w_fp[537], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[461] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[464] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[469] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[475] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[193], w_fp[582], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[460] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[469] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[471] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7050 OF 15495 ***
    // Wavefunction(s) for diagram number 7050
    // (none)
    // Amplitude(s) for diagram number 7050
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[583], w_fp[194], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[373] += amp_sv[0];
    jamp_sv[376] -= amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[418] += amp_sv[0];
    jamp_sv[419] -= amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[429] += amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[436] += amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[584], w_fp[194], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[376] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[416] += amp_sv[0];
    jamp_sv[418] += amp_sv[0];
    jamp_sv[419] -= amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[427] -= amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[436] += amp_sv[0];
    jamp_sv[460] += amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[561], w_fp[194], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[373] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[410] += amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[416] += amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[427] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[429] -= amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[460] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[466] += amp_sv[0];
    jamp_sv[474] -= amp_sv[0];

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
  diagramgroup706( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 169 );
    retrieveWf( wfs, w_cx, nevt, 190 );
    retrieveWf( wfs, w_cx, nevt, 193 );
    retrieveWf( wfs, w_cx, nevt, 194 );
    retrieveWf( wfs, w_cx, nevt, 196 );
    retrieveWf( wfs, w_cx, nevt, 199 );
    retrieveWf( wfs, w_cx, nevt, 201 );
    retrieveWf( wfs, w_cx, nevt, 355 );
    retrieveWf( wfs, w_cx, nevt, 450 );
    retrieveWf( wfs, w_cx, nevt, 477 );
    retrieveWf( wfs, w_cx, nevt, 486 );
    retrieveWf( wfs, w_cx, nevt, 520 );
    retrieveWf( wfs, w_cx, nevt, 529 );
    retrieveWf( wfs, w_cx, nevt, 535 );
    retrieveWf( wfs, w_cx, nevt, 561 );
    retrieveWf( wfs, w_cx, nevt, 576 );
    retrieveWf( wfs, w_cx, nevt, 583 );
    retrieveWf( wfs, w_cx, nevt, 584 );
#endif
#endif

    // *** DIAGRAM 7051 OF 15495 ***
    // Wavefunction(s) for diagram number 7051
    // (none)
    // Amplitude(s) for diagram number 7051
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[190], w_fp[583], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[419] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[429] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[190], w_fp[584], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[411] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[419] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[427] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[190], w_fp[561], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[410] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[427] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[429] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7052 OF 15495 ***
    // Wavefunction(s) for diagram number 7052
    // (none)
    // Amplitude(s) for diagram number 7052
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[477], w_fp[520], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[386] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[477], w_fp[486], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[387] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[390] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[477], w_fp[576], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[386] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[390] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7053 OF 15495 ***
    // Wavefunction(s) for diagram number 7053
    // (none)
    // Amplitude(s) for diagram number 7053
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[520], w_fp[1], w_fp[194], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[377] += amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[394] -= amp_sv[0];
    jamp_sv[395] += amp_sv[0];
    jamp_sv[396] -= amp_sv[0];
    jamp_sv[397] += amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[427] += amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[435] += amp_sv[0];
    jamp_sv[437] -= amp_sv[0];
    jamp_sv[469] -= amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[486], w_fp[1], w_fp[194], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[377] += amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[390] += amp_sv[0];
    jamp_sv[392] -= amp_sv[0];
    jamp_sv[394] -= amp_sv[0];
    jamp_sv[395] += amp_sv[0];
    jamp_sv[397] += amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[427] += amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[437] -= amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    jamp_sv[475] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[576], w_fp[1], w_fp[194], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[375] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[390] += amp_sv[0];
    jamp_sv[392] -= amp_sv[0];
    jamp_sv[396] += amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[404] -= amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    jamp_sv[469] += amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[475] += amp_sv[0];

    // *** DIAGRAM 7054 OF 15495 ***
    // Wavefunction(s) for diagram number 7054
    // (none)
    // Amplitude(s) for diagram number 7054
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[450], w_fp[201], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[433] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[460] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7055 OF 15495 ***
    // Wavefunction(s) for diagram number 7055
    // (none)
    // Amplitude(s) for diagram number 7055
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[193], w_fp[450], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[460] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[466] += amp_sv[0];
    jamp_sv[474] -= amp_sv[0];

    // *** DIAGRAM 7056 OF 15495 ***
    // Wavefunction(s) for diagram number 7056
    // (none)
    // Amplitude(s) for diagram number 7056
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[199], w_fp[169], w_fp[450], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[373] += amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[397] += amp_sv[0];
    jamp_sv[433] -= amp_sv[0];

    // *** DIAGRAM 7057 OF 15495 ***
    // Wavefunction(s) for diagram number 7057
    // (none)
    // Amplitude(s) for diagram number 7057
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[355], w_fp[535], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[376] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[436] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7058 OF 15495 ***
    // Wavefunction(s) for diagram number 7058
    // (none)
    // Amplitude(s) for diagram number 7058
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[199], w_fp[535], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[433] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7059 OF 15495 ***
    // Wavefunction(s) for diagram number 7059
    // (none)
    // Amplitude(s) for diagram number 7059
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[529], w_fp[477], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[403] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7060 OF 15495 ***
    // Wavefunction(s) for diagram number 7060
    // (none)
    // Amplitude(s) for diagram number 7060
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[529], w_fp[193], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[463] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] -= cxtype( 0, 1 ) * amp_sv[0];

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
  diagramgroup707( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 169 );
    retrieveWf( wfs, w_cx, nevt, 179 );
    retrieveWf( wfs, w_cx, nevt, 190 );
    retrieveWf( wfs, w_cx, nevt, 193 );
    retrieveWf( wfs, w_cx, nevt, 196 );
    retrieveWf( wfs, w_cx, nevt, 199 );
    retrieveWf( wfs, w_cx, nevt, 201 );
    retrieveWf( wfs, w_cx, nevt, 204 );
    retrieveWf( wfs, w_cx, nevt, 205 );
    retrieveWf( wfs, w_cx, nevt, 252 );
    retrieveWf( wfs, w_cx, nevt, 355 );
    retrieveWf( wfs, w_cx, nevt, 450 );
    retrieveWf( wfs, w_cx, nevt, 477 );
    retrieveWf( wfs, w_cx, nevt, 530 );
    retrieveWf( wfs, w_cx, nevt, 535 );
    retrieveWf( wfs, w_cx, nevt, 547 );
    retrieveWf( wfs, w_cx, nevt, 561 );
    retrieveWf( wfs, w_cx, nevt, 583 );
    retrieveWf( wfs, w_cx, nevt, 584 );
#endif
#endif

    // *** DIAGRAM 7061 OF 15495 ***
    // Wavefunction(s) for diagram number 7061
    // (none)
    // Amplitude(s) for diagram number 7061
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[477], w_fp[547], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[387] += amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[406] += amp_sv[0];

    // *** DIAGRAM 7062 OF 15495 ***
    // Wavefunction(s) for diagram number 7062
    // (none)
    // Amplitude(s) for diagram number 7062
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[355], w_fp[169], w_fp[547], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[376] += amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[460] -= amp_sv[0];
    jamp_sv[474] += amp_sv[0];

    // *** DIAGRAM 7063 OF 15495 ***
    // Wavefunction(s) for diagram number 7063
    // (none)
    // Amplitude(s) for diagram number 7063
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[547], w_fp[1], w_fp[201], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[376] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[436] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[460] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7064 OF 15495 ***
    // Wavefunction(s) for diagram number 7064
    // (none)
    // Amplitude(s) for diagram number 7064
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[199], w_fp[477], w_fp[530], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[387] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7065 OF 15495 ***
    // Wavefunction(s) for diagram number 7065
    // (none)
    // Amplitude(s) for diagram number 7065
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[355], w_fp[193], w_fp[530], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[460] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7066 OF 15495 ***
    // Wavefunction(s) for diagram number 7066
    // (none)
    // Amplitude(s) for diagram number 7066
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[169], w_fp[583], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[373] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[433] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[436] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[169], w_fp[584], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[376] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[436] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[460] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[169], w_fp[561], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[433] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[460] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7067 OF 15495 ***
    // Wavefunction(s) for diagram number 7067
    // (none)
    // Amplitude(s) for diagram number 7067
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[450], w_fp[205], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[372] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[432] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7068 OF 15495 ***
    // Wavefunction(s) for diagram number 7068
    // (none)
    // Amplitude(s) for diagram number 7068
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[179], w_fp[190], w_fp[450], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[410] += amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[416] += amp_sv[0];
    jamp_sv[420] -= amp_sv[0];

    // *** DIAGRAM 7069 OF 15495 ***
    // Wavefunction(s) for diagram number 7069
    // (none)
    // Amplitude(s) for diagram number 7069
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[204], w_fp[169], w_fp[450], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[372] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[396] += amp_sv[0];
    jamp_sv[432] -= amp_sv[0];

    // *** DIAGRAM 7070 OF 15495 ***
    // Wavefunction(s) for diagram number 7070
    // (none)
    // Amplitude(s) for diagram number 7070
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[252], w_fp[535], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[374] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[434] -= cxtype( 0, 1 ) * amp_sv[0];

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
  diagramgroup708( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 102 );
    retrieveWf( wfs, w_cx, nevt, 169 );
    retrieveWf( wfs, w_cx, nevt, 179 );
    retrieveWf( wfs, w_cx, nevt, 190 );
    retrieveWf( wfs, w_cx, nevt, 194 );
    retrieveWf( wfs, w_cx, nevt, 204 );
    retrieveWf( wfs, w_cx, nevt, 205 );
    retrieveWf( wfs, w_cx, nevt, 252 );
    retrieveWf( wfs, w_cx, nevt, 450 );
    retrieveWf( wfs, w_cx, nevt, 477 );
    retrieveWf( wfs, w_cx, nevt, 516 );
    retrieveWf( wfs, w_cx, nevt, 523 );
    retrieveWf( wfs, w_cx, nevt, 530 );
    retrieveWf( wfs, w_cx, nevt, 535 );
    retrieveWf( wfs, w_cx, nevt, 537 );
    retrieveWf( wfs, w_cx, nevt, 557 );
    retrieveWf( wfs, w_cx, nevt, 582 );
#endif
#endif

    // *** DIAGRAM 7071 OF 15495 ***
    // Wavefunction(s) for diagram number 7071
    // (none)
    // Amplitude(s) for diagram number 7071
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[204], w_fp[535], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[372] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[432] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7072 OF 15495 ***
    // Wavefunction(s) for diagram number 7072
    // (none)
    // Amplitude(s) for diagram number 7072
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[557], w_fp[477], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[390] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7073 OF 15495 ***
    // Wavefunction(s) for diagram number 7073
    // (none)
    // Amplitude(s) for diagram number 7073
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[557], w_fp[190], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[414] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7074 OF 15495 ***
    // Wavefunction(s) for diagram number 7074
    // (none)
    // Amplitude(s) for diagram number 7074
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[179], w_fp[477], w_fp[523], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[386] += amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[396] -= amp_sv[0];

    // *** DIAGRAM 7075 OF 15495 ***
    // Wavefunction(s) for diagram number 7075
    // (none)
    // Amplitude(s) for diagram number 7075
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[252], w_fp[169], w_fp[523], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[374] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[434] -= amp_sv[0];

    // *** DIAGRAM 7076 OF 15495 ***
    // Wavefunction(s) for diagram number 7076
    // (none)
    // Amplitude(s) for diagram number 7076
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[523], w_fp[1], w_fp[205], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[374] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[390] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[434] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7077 OF 15495 ***
    // Wavefunction(s) for diagram number 7077
    // (none)
    // Amplitude(s) for diagram number 7077
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[204], w_fp[477], w_fp[530], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[386] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7078 OF 15495 ***
    // Wavefunction(s) for diagram number 7078
    // (none)
    // Amplitude(s) for diagram number 7078
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[252], w_fp[190], w_fp[530], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[410] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7079 OF 15495 ***
    // Wavefunction(s) for diagram number 7079
    // (none)
    // Amplitude(s) for diagram number 7079
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[179], w_fp[169], w_fp[516], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[372] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[390] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[432] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[434] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[179], w_fp[169], w_fp[537], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[374] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[390] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[434] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[179], w_fp[169], w_fp[582], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[372] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[432] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7080 OF 15495 ***
    // Wavefunction(s) for diagram number 7080
    // (none)
    // Amplitude(s) for diagram number 7080
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[450], w_fp[194], w_fp[102], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[372] += amp_sv[0];
    jamp_sv[373] -= amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[396] += amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[427] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[429] -= amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[432] -= amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[469] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[471] += amp_sv[0];
    jamp_sv[472] -= amp_sv[0];

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
  diagramgroup709( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 102 );
    retrieveWf( wfs, w_cx, nevt, 169 );
    retrieveWf( wfs, w_cx, nevt, 194 );
    retrieveWf( wfs, w_cx, nevt, 208 );
    retrieveWf( wfs, w_cx, nevt, 209 );
    retrieveWf( wfs, w_cx, nevt, 363 );
    retrieveWf( wfs, w_cx, nevt, 450 );
    retrieveWf( wfs, w_cx, nevt, 477 );
    retrieveWf( wfs, w_cx, nevt, 530 );
    retrieveWf( wfs, w_cx, nevt, 532 );
    retrieveWf( wfs, w_cx, nevt, 535 );
    retrieveWf( wfs, w_cx, nevt, 556 );
#endif
#endif

    // *** DIAGRAM 7081 OF 15495 ***
    // Wavefunction(s) for diagram number 7081
    // (none)
    // Amplitude(s) for diagram number 7081
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[209], w_fp[450], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[427] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[429] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[469] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[471] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7082 OF 15495 ***
    // Wavefunction(s) for diagram number 7082
    // (none)
    // Amplitude(s) for diagram number 7082
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[208], w_fp[169], w_fp[450], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[372] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[432] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[433] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7083 OF 15495 ***
    // Wavefunction(s) for diagram number 7083
    // (none)
    // Amplitude(s) for diagram number 7083
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[535], w_fp[363], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[372] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[432] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[433] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[437] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7084 OF 15495 ***
    // Wavefunction(s) for diagram number 7084
    // (none)
    // Amplitude(s) for diagram number 7084
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[208], w_fp[535], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[372] += amp_sv[0];
    jamp_sv[373] -= amp_sv[0];
    jamp_sv[432] -= amp_sv[0];
    jamp_sv[433] += amp_sv[0];

    // *** DIAGRAM 7085 OF 15495 ***
    // Wavefunction(s) for diagram number 7085
    // (none)
    // Amplitude(s) for diagram number 7085
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[532], w_fp[477], w_fp[102], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[394] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[404] -= amp_sv[0];
    jamp_sv[405] += amp_sv[0];

    // *** DIAGRAM 7086 OF 15495 ***
    // Wavefunction(s) for diagram number 7086
    // (none)
    // Amplitude(s) for diagram number 7086
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[532], w_fp[169], w_fp[363], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[394] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[429] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[471] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7087 OF 15495 ***
    // Wavefunction(s) for diagram number 7087
    // (none)
    // Amplitude(s) for diagram number 7087
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[532], w_fp[209], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[428] += amp_sv[0];
    jamp_sv[429] -= amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[471] += amp_sv[0];

    // *** DIAGRAM 7088 OF 15495 ***
    // Wavefunction(s) for diagram number 7088
    // (none)
    // Amplitude(s) for diagram number 7088
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[477], w_fp[556], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[386] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7089 OF 15495 ***
    // Wavefunction(s) for diagram number 7089
    // (none)
    // Amplitude(s) for diagram number 7089
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[556], w_fp[1], w_fp[194], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[377] += amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[394] -= amp_sv[0];
    jamp_sv[395] += amp_sv[0];
    jamp_sv[396] -= amp_sv[0];
    jamp_sv[397] += amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[427] += amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[435] += amp_sv[0];
    jamp_sv[437] -= amp_sv[0];
    jamp_sv[469] -= amp_sv[0];
    jamp_sv[472] += amp_sv[0];

    // *** DIAGRAM 7090 OF 15495 ***
    // Wavefunction(s) for diagram number 7090
    // (none)
    // Amplitude(s) for diagram number 7090
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[208], w_fp[477], w_fp[530], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[386] += amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[396] -= amp_sv[0];
    jamp_sv[397] += amp_sv[0];

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
  diagramgroup710( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 169 );
    retrieveWf( wfs, w_cx, nevt, 194 );
    retrieveWf( wfs, w_cx, nevt, 225 );
    retrieveWf( wfs, w_cx, nevt, 226 );
    retrieveWf( wfs, w_cx, nevt, 228 );
    retrieveWf( wfs, w_cx, nevt, 363 );
    retrieveWf( wfs, w_cx, nevt, 437 );
    retrieveWf( wfs, w_cx, nevt, 447 );
    retrieveWf( wfs, w_cx, nevt, 450 );
    retrieveWf( wfs, w_cx, nevt, 473 );
    retrieveWf( wfs, w_cx, nevt, 524 );
    retrieveWf( wfs, w_cx, nevt, 530 );
    retrieveWf( wfs, w_cx, nevt, 531 );
    retrieveWf( wfs, w_cx, nevt, 541 );
    retrieveWf( wfs, w_cx, nevt, 546 );
    retrieveWf( wfs, w_cx, nevt, 585 );
#endif
#endif

    // *** DIAGRAM 7091 OF 15495 ***
    // Wavefunction(s) for diagram number 7091
    // (none)
    // Amplitude(s) for diagram number 7091
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[530], w_fp[363], w_fp[194], COUPs[0], 1.0, &amp_fp[0] );
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

    // *** DIAGRAM 7092 OF 15495 ***
    // Wavefunction(s) for diagram number 7092
    // (none)
    // Amplitude(s) for diagram number 7092
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[169], w_fp[546], COUPs[1], 1.0, &amp_fp[0] );
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
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[169], w_fp[541], COUPs[1], 1.0, &amp_fp[0] );
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
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[169], w_fp[585], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[372] += amp_sv[0];
    jamp_sv[373] -= amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[396] += amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[427] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[429] -= amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[432] -= amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[469] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[471] += amp_sv[0];
    jamp_sv[472] -= amp_sv[0];

    // *** DIAGRAM 7093 OF 15495 ***
    // Wavefunction(s) for diagram number 7093
    // (none)
    // Amplitude(s) for diagram number 7093
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[450], w_fp[228], w_fp[4], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[618] += amp_sv[0];
    jamp_sv[619] -= amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[662] += amp_sv[0];
    jamp_sv[663] -= amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[685] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];
    jamp_sv[687] += amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[696] -= amp_sv[0];
    jamp_sv[697] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[450], w_fp[228], w_fp[4], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[618] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    jamp_sv[652] -= amp_sv[0];
    jamp_sv[655] += amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[676] -= amp_sv[0];
    jamp_sv[679] += amp_sv[0];
    jamp_sv[682] -= amp_sv[0];
    jamp_sv[685] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];
    jamp_sv[687] += amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[696] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[450], w_fp[228], w_fp[4], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[619] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    jamp_sv[652] -= amp_sv[0];
    jamp_sv[655] += amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    jamp_sv[661] += amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    jamp_sv[663] += amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[676] -= amp_sv[0];
    jamp_sv[679] += amp_sv[0];
    jamp_sv[682] -= amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[697] -= amp_sv[0];

    // *** DIAGRAM 7094 OF 15495 ***
    // Wavefunction(s) for diagram number 7094
    // (none)
    // Amplitude(s) for diagram number 7094
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[228], w_fp[5], w_fp[531], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[618] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    jamp_sv[652] -= amp_sv[0];
    jamp_sv[655] += amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[676] -= amp_sv[0];
    jamp_sv[679] += amp_sv[0];
    jamp_sv[682] -= amp_sv[0];
    jamp_sv[685] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];
    jamp_sv[687] += amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[696] -= amp_sv[0];

    // *** DIAGRAM 7095 OF 15495 ***
    // Wavefunction(s) for diagram number 7095
    // (none)
    // Amplitude(s) for diagram number 7095
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[228], w_fp[4], w_fp[447], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[619] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    jamp_sv[652] -= amp_sv[0];
    jamp_sv[655] += amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    jamp_sv[661] += amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    jamp_sv[663] += amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[676] -= amp_sv[0];
    jamp_sv[679] += amp_sv[0];
    jamp_sv[682] -= amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[697] -= amp_sv[0];

    // *** DIAGRAM 7096 OF 15495 ***
    // Wavefunction(s) for diagram number 7096
    // (none)
    // Amplitude(s) for diagram number 7096
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[437], w_fp[225], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[661] += amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    jamp_sv[663] += amp_sv[0];
    jamp_sv[664] -= amp_sv[0];

    // *** DIAGRAM 7097 OF 15495 ***
    // Wavefunction(s) for diagram number 7097
    // (none)
    // Amplitude(s) for diagram number 7097
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[225], w_fp[447], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[652] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[655] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[662] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[663] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[664] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7098 OF 15495 ***
    // Wavefunction(s) for diagram number 7098
    // (none)
    // Amplitude(s) for diagram number 7098
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[437], w_fp[226], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[685] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];
    jamp_sv[687] += amp_sv[0];
    jamp_sv[688] -= amp_sv[0];

    // *** DIAGRAM 7099 OF 15495 ***
    // Wavefunction(s) for diagram number 7099
    // (none)
    // Amplitude(s) for diagram number 7099
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[226], w_fp[531], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[679] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[687] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[688] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[690] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7100 OF 15495 ***
    // Wavefunction(s) for diagram number 7100
    // (none)
    // Amplitude(s) for diagram number 7100
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[524], w_fp[473], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[638] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] -= cxtype( 0, 1 ) * amp_sv[0];

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
