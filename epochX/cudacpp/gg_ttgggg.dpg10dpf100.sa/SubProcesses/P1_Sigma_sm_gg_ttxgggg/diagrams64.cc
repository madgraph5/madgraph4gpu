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
  diagramgroup631( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 100 );
    retrieveWf( wfs, w_cx, nevt, 102 );
    retrieveWf( wfs, w_cx, nevt, 103 );
    retrieveWf( wfs, w_cx, nevt, 124 );
    retrieveWf( wfs, w_cx, nevt, 186 );
    retrieveWf( wfs, w_cx, nevt, 208 );
    retrieveWf( wfs, w_cx, nevt, 259 );
    retrieveWf( wfs, w_cx, nevt, 436 );
    retrieveWf( wfs, w_cx, nevt, 442 );
    retrieveWf( wfs, w_cx, nevt, 443 );
    retrieveWf( wfs, w_cx, nevt, 456 );
    retrieveWf( wfs, w_cx, nevt, 463 );
    retrieveWf( wfs, w_cx, nevt, 474 );
    retrieveWf( wfs, w_cx, nevt, 530 );
    retrieveWf( wfs, w_cx, nevt, 532 );
    retrieveWf( wfs, w_cx, nevt, 556 );
#endif
#endif

    // *** DIAGRAM 6301 OF 15495 ***
    // Wavefunction(s) for diagram number 6301
    // (none)
    // Amplitude(s) for diagram number 6301
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[436], w_fp[556], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[179] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[180] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[181] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6302 OF 15495 ***
    // Wavefunction(s) for diagram number 6302
    // (none)
    // Amplitude(s) for diagram number 6302
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[474], w_fp[456], w_fp[102], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[134] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[163] -= amp_sv[0];
    jamp_sv[164] += amp_sv[0];
    jamp_sv[165] -= amp_sv[0];
    jamp_sv[166] += amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[180] += amp_sv[0];
    jamp_sv[181] -= amp_sv[0];
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[195] += amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[225] += amp_sv[0];
    jamp_sv[226] -= amp_sv[0];

    // *** DIAGRAM 6303 OF 15495 ***
    // Wavefunction(s) for diagram number 6303
    // (none)
    // Amplitude(s) for diagram number 6303
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[463], w_fp[474], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[163] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[165] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6304 OF 15495 ***
    // Wavefunction(s) for diagram number 6304
    // (none)
    // Amplitude(s) for diagram number 6304
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[208], w_fp[259], w_fp[474], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[134] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[180] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[181] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6305 OF 15495 ***
    // Wavefunction(s) for diagram number 6305
    // (none)
    // Amplitude(s) for diagram number 6305
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[530], w_fp[456], w_fp[103], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[133] += amp_sv[0];
    jamp_sv[134] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[196] += amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[225] -= amp_sv[0];

    // *** DIAGRAM 6306 OF 15495 ***
    // Wavefunction(s) for diagram number 6306
    // (none)
    // Amplitude(s) for diagram number 6306
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[208], w_fp[436], w_fp[530], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[170] += amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[180] -= amp_sv[0];
    jamp_sv[181] += amp_sv[0];

    // *** DIAGRAM 6307 OF 15495 ***
    // Wavefunction(s) for diagram number 6307
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[530], w_fp[102], w_fp[5], COUPs[2], 1.0, 0., 0., w_fp[550] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[530], w_fp[102], w_fp[5], COUPs[2], 1.0, 0., 0., w_fp[480] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[530], w_fp[102], w_fp[5], COUPs[2], 1.0, 0., 0., w_fp[558] );
    // Amplitude(s) for diagram number 6307
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[259], w_fp[550], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[134] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[136] += amp_sv[0];
    jamp_sv[164] += amp_sv[0];
    jamp_sv[165] -= amp_sv[0];
    jamp_sv[178] -= amp_sv[0];
    jamp_sv[179] += amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[193] += amp_sv[0];
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[195] += amp_sv[0];
    jamp_sv[196] -= amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[225] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[259], w_fp[480], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[134] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[163] -= amp_sv[0];
    jamp_sv[164] += amp_sv[0];
    jamp_sv[165] -= amp_sv[0];
    jamp_sv[166] += amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[180] += amp_sv[0];
    jamp_sv[181] -= amp_sv[0];
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[195] += amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[225] += amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[259], w_fp[558], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[133] += amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[163] -= amp_sv[0];
    jamp_sv[166] += amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[180] += amp_sv[0];
    jamp_sv[181] -= amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[196] += amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[226] -= amp_sv[0];

    // *** DIAGRAM 6308 OF 15495 ***
    // Wavefunction(s) for diagram number 6308
    // (none)
    // Amplitude(s) for diagram number 6308
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[186], w_fp[442], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[193] += amp_sv[0];

    // *** DIAGRAM 6309 OF 15495 ***
    // Wavefunction(s) for diagram number 6309
    // (none)
    // Amplitude(s) for diagram number 6309
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[442], w_fp[124], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6310 OF 15495 ***
    // Wavefunction(s) for diagram number 6310
    // (none)
    // Amplitude(s) for diagram number 6310
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[532], w_fp[443], w_fp[100], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] += amp_sv[0];
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[165] += amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 480 );
    storeWf( wfs, w_cx, nevt, 550 );
    storeWf( wfs, w_cx, nevt, 558 );
#endif
#endif
  }

  //--------------------------------------------------------------------------

#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup632( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 100 );
    retrieveWf( wfs, w_cx, nevt, 124 );
    retrieveWf( wfs, w_cx, nevt, 186 );
    retrieveWf( wfs, w_cx, nevt, 259 );
    retrieveWf( wfs, w_cx, nevt, 443 );
    retrieveWf( wfs, w_cx, nevt, 456 );
    retrieveWf( wfs, w_cx, nevt, 466 );
    retrieveWf( wfs, w_cx, nevt, 510 );
    retrieveWf( wfs, w_cx, nevt, 523 );
    retrieveWf( wfs, w_cx, nevt, 530 );
    retrieveWf( wfs, w_cx, nevt, 532 );
#endif
#endif

    // *** DIAGRAM 6311 OF 15495 ***
    // Wavefunction(s) for diagram number 6311
    // (none)
    // Amplitude(s) for diagram number 6311
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[532], w_fp[466], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[188] += amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[231] += amp_sv[0];

    // *** DIAGRAM 6312 OF 15495 ***
    // Wavefunction(s) for diagram number 6312
    // (none)
    // Amplitude(s) for diagram number 6312
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[532], w_fp[259], w_fp[124], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[155] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[165] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6313 OF 15495 ***
    // Wavefunction(s) for diagram number 6313
    // (none)
    // Amplitude(s) for diagram number 6313
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[523], w_fp[456], w_fp[100], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[156] += amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[193] += amp_sv[0];
    jamp_sv[229] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[231] += amp_sv[0];
    jamp_sv[232] -= amp_sv[0];

    // *** DIAGRAM 6314 OF 15495 ***
    // Wavefunction(s) for diagram number 6314
    // (none)
    // Amplitude(s) for diagram number 6314
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[466], w_fp[523], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[187] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6315 OF 15495 ***
    // Wavefunction(s) for diagram number 6315
    // (none)
    // Amplitude(s) for diagram number 6315
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[186], w_fp[259], w_fp[523], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[146] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[147] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6316 OF 15495 ***
    // Wavefunction(s) for diagram number 6316
    // (none)
    // Amplitude(s) for diagram number 6316
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[510], w_fp[456], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[135] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[154] += amp_sv[0];
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[156] += amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[197] += amp_sv[0];
    jamp_sv[229] += amp_sv[0];
    jamp_sv[232] -= amp_sv[0];

    // *** DIAGRAM 6317 OF 15495 ***
    // Wavefunction(s) for diagram number 6317
    // (none)
    // Amplitude(s) for diagram number 6317
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[443], w_fp[510], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[146] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[147] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[154] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[155] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[165] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6318 OF 15495 ***
    // Wavefunction(s) for diagram number 6318
    // (none)
    // Amplitude(s) for diagram number 6318
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[530], w_fp[456], w_fp[124], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[155] += amp_sv[0];
    jamp_sv[164] += amp_sv[0];
    jamp_sv[165] -= amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[193] += amp_sv[0];
    jamp_sv[195] += amp_sv[0];
    jamp_sv[197] -= amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[231] += amp_sv[0];

    // *** DIAGRAM 6319 OF 15495 ***
    // Wavefunction(s) for diagram number 6319
    // (none)
    // Amplitude(s) for diagram number 6319
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[186], w_fp[443], w_fp[530], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[146] += amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[156] -= amp_sv[0];
    jamp_sv[157] += amp_sv[0];

    // *** DIAGRAM 6320 OF 15495 ***
    // Wavefunction(s) for diagram number 6320
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[530], w_fp[4], w_fp[100], COUPs[2], 1.0, 0., 0., w_fp[559] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[530], w_fp[4], w_fp[100], COUPs[2], 1.0, 0., 0., w_fp[560] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[530], w_fp[4], w_fp[100], COUPs[2], 1.0, 0., 0., w_fp[438] );
    // Amplitude(s) for diagram number 6320
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[259], w_fp[559], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[154] += amp_sv[0];
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    jamp_sv[192] += amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[197] += amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[259], w_fp[560], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[135] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[154] += amp_sv[0];
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[156] += amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[197] += amp_sv[0];
    jamp_sv[229] += amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[259], w_fp[438], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[156] += amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[193] += amp_sv[0];
    jamp_sv[229] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[231] += amp_sv[0];
    jamp_sv[232] -= amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 438 );
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
  diagramgroup633( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 29 );
    retrieveWf( wfs, w_cx, nevt, 77 );
    retrieveWf( wfs, w_cx, nevt, 94 );
    retrieveWf( wfs, w_cx, nevt, 139 );
    retrieveWf( wfs, w_cx, nevt, 140 );
    retrieveWf( wfs, w_cx, nevt, 141 );
    retrieveWf( wfs, w_cx, nevt, 259 );
    retrieveWf( wfs, w_cx, nevt, 442 );
    retrieveWf( wfs, w_cx, nevt, 523 );
    retrieveWf( wfs, w_cx, nevt, 530 );
    retrieveWf( wfs, w_cx, nevt, 532 );
    retrieveWf( wfs, w_cx, nevt, 534 );
#endif
#endif

    // *** DIAGRAM 6321 OF 15495 ***
    // Wavefunction(s) for diagram number 6321
    // (none)
    // Amplitude(s) for diagram number 6321
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[442], w_fp[139], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[442], w_fp[140], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[133] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[134] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[196] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[442], w_fp[141], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[134] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[196] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6322 OF 15495 ***
    // Wavefunction(s) for diagram number 6322
    // (none)
    // Amplitude(s) for diagram number 6322
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[532], w_fp[259], w_fp[139], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[155] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[165] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[532], w_fp[259], w_fp[140], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[165] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[179] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[532], w_fp[259], w_fp[141], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[155] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[179] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6323 OF 15495 ***
    // Wavefunction(s) for diagram number 6323
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[530], w_fp[139], COUPs[0], 1.0, 0., 0., w_fp[442] );
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[530], w_fp[140], COUPs[0], 1.0, 0., 0., w_fp[562] );
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[530], w_fp[141], COUPs[0], 1.0, 0., 0., w_fp[437] );
    // Amplitude(s) for diagram number 6323
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[259], w_fp[442], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[154] += amp_sv[0];
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    jamp_sv[192] += amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[197] += amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[259], w_fp[562], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[133] += amp_sv[0];
    jamp_sv[134] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[196] += amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[225] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[259], w_fp[437], COUPs[1], 1.0, &amp_fp[0] );
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

    // *** DIAGRAM 6324 OF 15495 ***
    // Wavefunction(s) for diagram number 6324
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[530], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[449] );
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[449], w_fp[5], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[590] );
    // Amplitude(s) for diagram number 6324
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[94], w_fp[590], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[89] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6325 OF 15495 ***
    // Wavefunction(s) for diagram number 6325
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[449], w_fp[7], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[481] );
    // Amplitude(s) for diagram number 6325
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[94], w_fp[481], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[95] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6326 OF 15495 ***
    // Wavefunction(s) for diagram number 6326
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[449], w_fp[4], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[596] );
    // Amplitude(s) for diagram number 6326
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[29], w_fp[596], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[83] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[491] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6327 OF 15495 ***
    // Wavefunction(s) for diagram number 6327
    // (none)
    // Amplitude(s) for diagram number 6327
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[29], w_fp[481], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[93] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[501] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6328 OF 15495 ***
    // Wavefunction(s) for diagram number 6328
    // (none)
    // Amplitude(s) for diagram number 6328
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[77], w_fp[596], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[81] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6329 OF 15495 ***
    // Wavefunction(s) for diagram number 6329
    // (none)
    // Amplitude(s) for diagram number 6329
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[77], w_fp[590], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[87] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6330 OF 15495 ***
    // Wavefunction(s) for diagram number 6330
    // (none)
    // Amplitude(s) for diagram number 6330
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[523], w_fp[534], w_fp[5], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[81] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    jamp_sv[469] += amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[475] += amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[491] += amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[523], w_fp[534], w_fp[5], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[81] += amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[435] += amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[523], w_fp[534], w_fp[5], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[83] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[435] += amp_sv[0];
    jamp_sv[461] += amp_sv[0];
    jamp_sv[469] -= amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[475] -= amp_sv[0];
    jamp_sv[491] -= amp_sv[0];
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 442 );
    storeWf( wfs, w_cx, nevt, 449 );
    storeWf( wfs, w_cx, nevt, 481 );
    storeWf( wfs, w_cx, nevt, 562 );
    storeWf( wfs, w_cx, nevt, 590 );
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
  diagramgroup634( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 29 );
    retrieveWf( wfs, w_cx, nevt, 77 );
    retrieveWf( wfs, w_cx, nevt, 94 );
    retrieveWf( wfs, w_cx, nevt, 474 );
    retrieveWf( wfs, w_cx, nevt, 479 );
    retrieveWf( wfs, w_cx, nevt, 515 );
    retrieveWf( wfs, w_cx, nevt, 523 );
    retrieveWf( wfs, w_cx, nevt, 534 );
    retrieveWf( wfs, w_cx, nevt, 536 );
    retrieveWf( wfs, w_cx, nevt, 580 );
#endif
#endif

    // *** DIAGRAM 6331 OF 15495 ***
    // Wavefunction(s) for diagram number 6331
    // (none)
    // Amplitude(s) for diagram number 6331
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[534], w_fp[7], w_fp[580], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[81] += amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[435] += amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];

    // *** DIAGRAM 6332 OF 15495 ***
    // Wavefunction(s) for diagram number 6332
    // (none)
    // Amplitude(s) for diagram number 6332
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[534], w_fp[5], w_fp[515], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[83] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[435] += amp_sv[0];
    jamp_sv[461] += amp_sv[0];
    jamp_sv[469] -= amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[475] -= amp_sv[0];
    jamp_sv[491] -= amp_sv[0];
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];

    // *** DIAGRAM 6333 OF 15495 ***
    // Wavefunction(s) for diagram number 6333
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[523], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[595] );
    // Amplitude(s) for diagram number 6333
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[29], w_fp[595], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[83] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[491] -= amp_sv[0];

    // *** DIAGRAM 6334 OF 15495 ***
    // Wavefunction(s) for diagram number 6334
    // (none)
    // Amplitude(s) for diagram number 6334
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[29], w_fp[2], w_fp[515], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[83] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[491] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[653] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6335 OF 15495 ***
    // Wavefunction(s) for diagram number 6335
    // (none)
    // Amplitude(s) for diagram number 6335
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[77], w_fp[595], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[81] += amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[489] -= amp_sv[0];

    // *** DIAGRAM 6336 OF 15495 ***
    // Wavefunction(s) for diagram number 6336
    // (none)
    // Amplitude(s) for diagram number 6336
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[77], w_fp[2], w_fp[580], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[81] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6337 OF 15495 ***
    // Wavefunction(s) for diagram number 6337
    // (none)
    // Amplitude(s) for diagram number 6337
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[474], w_fp[534], w_fp[4], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[87] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[341] -= amp_sv[0];
    jamp_sv[349] += amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[355] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[377] += amp_sv[0];
    jamp_sv[435] += amp_sv[0];
    jamp_sv[437] -= amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[497] += amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[474], w_fp[534], w_fp[4], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[87] += amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[435] += amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[701] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[474], w_fp[534], w_fp[4], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[89] += amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[341] += amp_sv[0];
    jamp_sv[349] -= amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[355] -= amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[701] += amp_sv[0];

    // *** DIAGRAM 6338 OF 15495 ***
    // Wavefunction(s) for diagram number 6338
    // (none)
    // Amplitude(s) for diagram number 6338
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[534], w_fp[7], w_fp[479], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[87] += amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[435] += amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[701] += amp_sv[0];

    // *** DIAGRAM 6339 OF 15495 ***
    // Wavefunction(s) for diagram number 6339
    // (none)
    // Amplitude(s) for diagram number 6339
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[534], w_fp[4], w_fp[536], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[89] += amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[341] += amp_sv[0];
    jamp_sv[349] -= amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[355] -= amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[701] += amp_sv[0];

    // *** DIAGRAM 6340 OF 15495 ***
    // Wavefunction(s) for diagram number 6340
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[474], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[594] );
    // Amplitude(s) for diagram number 6340
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[94], w_fp[594], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[89] += amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[497] -= amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 594 );
    storeWf( wfs, w_cx, nevt, 595 );
#endif
#endif
  }

  //--------------------------------------------------------------------------


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup635( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 29 );
    retrieveWf( wfs, w_cx, nevt, 77 );
    retrieveWf( wfs, w_cx, nevt, 94 );
    retrieveWf( wfs, w_cx, nevt, 479 );
    retrieveWf( wfs, w_cx, nevt, 519 );
    retrieveWf( wfs, w_cx, nevt, 522 );
    retrieveWf( wfs, w_cx, nevt, 534 );
    retrieveWf( wfs, w_cx, nevt, 536 );
    retrieveWf( wfs, w_cx, nevt, 547 );
    retrieveWf( wfs, w_cx, nevt, 594 );
#endif
#endif

    // *** DIAGRAM 6341 OF 15495 ***
    // Wavefunction(s) for diagram number 6341
    // (none)
    // Amplitude(s) for diagram number 6341
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[94], w_fp[2], w_fp[536], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[89] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[437] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[623] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[691] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6342 OF 15495 ***
    // Wavefunction(s) for diagram number 6342
    // (none)
    // Amplitude(s) for diagram number 6342
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[77], w_fp[594], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[87] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[435] += amp_sv[0];
    jamp_sv[495] -= amp_sv[0];

    // *** DIAGRAM 6343 OF 15495 ***
    // Wavefunction(s) for diagram number 6343
    // (none)
    // Amplitude(s) for diagram number 6343
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[77], w_fp[2], w_fp[479], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[87] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6344 OF 15495 ***
    // Wavefunction(s) for diagram number 6344
    // (none)
    // Amplitude(s) for diagram number 6344
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[547], w_fp[534], w_fp[4], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[93] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[301] += amp_sv[0];
    jamp_sv[307] += amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[427] -= amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[501] -= amp_sv[0];
    jamp_sv[503] += amp_sv[0];
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[623] += amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    jamp_sv[701] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[547], w_fp[534], w_fp[4], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[93] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[341] += amp_sv[0];
    jamp_sv[355] -= amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[427] -= amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[461] += amp_sv[0];
    jamp_sv[475] -= amp_sv[0];
    jamp_sv[501] -= amp_sv[0];
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[547], w_fp[534], w_fp[4], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[95] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[307] -= amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[341] += amp_sv[0];
    jamp_sv[355] -= amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[461] += amp_sv[0];
    jamp_sv[475] -= amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[701] += amp_sv[0];

    // *** DIAGRAM 6345 OF 15495 ***
    // Wavefunction(s) for diagram number 6345
    // (none)
    // Amplitude(s) for diagram number 6345
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[534], w_fp[5], w_fp[522], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[93] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[341] += amp_sv[0];
    jamp_sv[355] -= amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[427] -= amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[461] += amp_sv[0];
    jamp_sv[475] -= amp_sv[0];
    jamp_sv[501] -= amp_sv[0];
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];

    // *** DIAGRAM 6346 OF 15495 ***
    // Wavefunction(s) for diagram number 6346
    // (none)
    // Amplitude(s) for diagram number 6346
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[534], w_fp[4], w_fp[519], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[95] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[307] -= amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[341] += amp_sv[0];
    jamp_sv[355] -= amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[461] += amp_sv[0];
    jamp_sv[475] -= amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[701] += amp_sv[0];

    // *** DIAGRAM 6347 OF 15495 ***
    // Wavefunction(s) for diagram number 6347
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[547], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[571] );
    // Amplitude(s) for diagram number 6347
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[94], w_fp[571], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[95] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[701] += amp_sv[0];

    // *** DIAGRAM 6348 OF 15495 ***
    // Wavefunction(s) for diagram number 6348
    // (none)
    // Amplitude(s) for diagram number 6348
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[94], w_fp[2], w_fp[519], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[95] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[437] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[475] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[623] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6349 OF 15495 ***
    // Wavefunction(s) for diagram number 6349
    // (none)
    // Amplitude(s) for diagram number 6349
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[29], w_fp[571], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[93] += amp_sv[0];
    jamp_sv[501] -= amp_sv[0];
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];

    // *** DIAGRAM 6350 OF 15495 ***
    // Wavefunction(s) for diagram number 6350
    // (none)
    // Amplitude(s) for diagram number 6350
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[29], w_fp[2], w_fp[522], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[93] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[501] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] -= cxtype( 0, 1 ) * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 571 );
#endif
#endif
  }

  //--------------------------------------------------------------------------


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup636( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 29 );
    retrieveWf( wfs, w_cx, nevt, 77 );
    retrieveWf( wfs, w_cx, nevt, 94 );
    retrieveWf( wfs, w_cx, nevt, 156 );
    retrieveWf( wfs, w_cx, nevt, 157 );
    retrieveWf( wfs, w_cx, nevt, 158 );
    retrieveWf( wfs, w_cx, nevt, 163 );
    retrieveWf( wfs, w_cx, nevt, 451 );
    retrieveWf( wfs, w_cx, nevt, 454 );
    retrieveWf( wfs, w_cx, nevt, 486 );
    retrieveWf( wfs, w_cx, nevt, 520 );
    retrieveWf( wfs, w_cx, nevt, 530 );
    retrieveWf( wfs, w_cx, nevt, 534 );
    retrieveWf( wfs, w_cx, nevt, 575 );
    retrieveWf( wfs, w_cx, nevt, 576 );
    retrieveWf( wfs, w_cx, nevt, 586 );
    retrieveWf( wfs, w_cx, nevt, 588 );
    retrieveWf( wfs, w_cx, nevt, 589 );
#endif
#endif

    // *** DIAGRAM 6351 OF 15495 ***
    // Wavefunction(s) for diagram number 6351
    // (none)
    // Amplitude(s) for diagram number 6351
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[588], w_fp[534], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[81] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[301] += amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[623] += amp_sv[0];
    jamp_sv[661] += amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    jamp_sv[701] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[586], w_fp[534], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[301] += amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[623] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];
    jamp_sv[661] += amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[701] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[534], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[489] += amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[685] += amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[699] -= amp_sv[0];

    // *** DIAGRAM 6352 OF 15495 ***
    // Wavefunction(s) for diagram number 6352
    // (none)
    // Amplitude(s) for diagram number 6352
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[77], w_fp[2], w_fp[588], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[81] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[77], w_fp[2], w_fp[586], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[87] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[77], w_fp[2], w_fp[451], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[81] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6353 OF 15495 ***
    // Wavefunction(s) for diagram number 6353
    // (none)
    // Amplitude(s) for diagram number 6353
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[520], w_fp[534], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[83] += amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[341] -= amp_sv[0];
    jamp_sv[355] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[377] += amp_sv[0];
    jamp_sv[427] += amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[435] += amp_sv[0];
    jamp_sv[437] -= amp_sv[0];
    jamp_sv[469] -= amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[491] -= amp_sv[0];
    jamp_sv[501] += amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[486], w_fp[534], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[341] -= amp_sv[0];
    jamp_sv[355] += amp_sv[0];
    jamp_sv[377] += amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[427] += amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[437] -= amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    jamp_sv[475] += amp_sv[0];
    jamp_sv[501] += amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[699] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[576], w_fp[534], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    jamp_sv[469] += amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[475] += amp_sv[0];
    jamp_sv[491] += amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[699] -= amp_sv[0];

    // *** DIAGRAM 6354 OF 15495 ***
    // Wavefunction(s) for diagram number 6354
    // (none)
    // Amplitude(s) for diagram number 6354
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[29], w_fp[2], w_fp[520], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[83] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[93] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[491] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[501] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[653] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[29], w_fp[2], w_fp[486], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[93] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[501] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[29], w_fp[2], w_fp[576], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[83] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[491] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[653] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6355 OF 15495 ***
    // Wavefunction(s) for diagram number 6355
    // (none)
    // Amplitude(s) for diagram number 6355
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[454], w_fp[534], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[89] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[307] += amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[349] -= amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    jamp_sv[475] += amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[503] += amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[575], w_fp[534], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[301] += amp_sv[0];
    jamp_sv[307] += amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[341] -= amp_sv[0];
    jamp_sv[355] += amp_sv[0];
    jamp_sv[377] += amp_sv[0];
    jamp_sv[437] -= amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    jamp_sv[475] += amp_sv[0];
    jamp_sv[503] += amp_sv[0];
    jamp_sv[623] += amp_sv[0];
    jamp_sv[701] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[589], w_fp[534], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[301] += amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[341] -= amp_sv[0];
    jamp_sv[349] += amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[355] += amp_sv[0];
    jamp_sv[377] += amp_sv[0];
    jamp_sv[437] -= amp_sv[0];
    jamp_sv[497] += amp_sv[0];
    jamp_sv[623] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[701] -= amp_sv[0];

    // *** DIAGRAM 6356 OF 15495 ***
    // Wavefunction(s) for diagram number 6356
    // (none)
    // Amplitude(s) for diagram number 6356
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[94], w_fp[2], w_fp[454], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[89] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[95] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[475] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[691] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[94], w_fp[2], w_fp[575], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[95] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[437] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[475] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[623] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[94], w_fp[2], w_fp[589], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[89] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[437] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[623] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[691] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6357 OF 15495 ***
    // Wavefunction(s) for diagram number 6357
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[157], w_fp[530], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[495] );
    // Amplitude(s) for diagram number 6357
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[495], w_fp[158], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[307] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[310] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6358 OF 15495 ***
    // Wavefunction(s) for diagram number 6358
    // (none)
    // Amplitude(s) for diagram number 6358
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[495], w_fp[163], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6359 OF 15495 ***
    // Wavefunction(s) for diagram number 6359
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[156], w_fp[530], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[553] );
    // Amplitude(s) for diagram number 6359
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[29], w_fp[553], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[257] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6360 OF 15495 ***
    // Wavefunction(s) for diagram number 6360
    // (none)
    // Amplitude(s) for diagram number 6360
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[77], w_fp[553], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[255] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 495 );
    storeWf( wfs, w_cx, nevt, 553 );
#endif
#endif
  }

  //--------------------------------------------------------------------------


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup637( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 29 );
    retrieveWf( wfs, w_cx, nevt, 77 );
    retrieveWf( wfs, w_cx, nevt, 100 );
    retrieveWf( wfs, w_cx, nevt, 156 );
    retrieveWf( wfs, w_cx, nevt, 157 );
    retrieveWf( wfs, w_cx, nevt, 158 );
    retrieveWf( wfs, w_cx, nevt, 163 );
    retrieveWf( wfs, w_cx, nevt, 454 );
    retrieveWf( wfs, w_cx, nevt, 474 );
    retrieveWf( wfs, w_cx, nevt, 495 );
    retrieveWf( wfs, w_cx, nevt, 497 );
    retrieveWf( wfs, w_cx, nevt, 530 );
    retrieveWf( wfs, w_cx, nevt, 547 );
    retrieveWf( wfs, w_cx, nevt, 575 );
    retrieveWf( wfs, w_cx, nevt, 589 );
#endif
#endif

    // *** DIAGRAM 6361 OF 15495 ***
    // Wavefunction(s) for diagram number 6361
    // (none)
    // Amplitude(s) for diagram number 6361
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[474], w_fp[497], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[255] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6362 OF 15495 ***
    // Wavefunction(s) for diagram number 6362
    // (none)
    // Amplitude(s) for diagram number 6362
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[77], w_fp[156], w_fp[474], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[255] += amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[301] += amp_sv[0];
    jamp_sv[315] -= amp_sv[0];

    // *** DIAGRAM 6363 OF 15495 ***
    // Wavefunction(s) for diagram number 6363
    // (none)
    // Amplitude(s) for diagram number 6363
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[157], w_fp[163], w_fp[474], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[341] += amp_sv[0];
    jamp_sv[349] -= amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[355] -= amp_sv[0];

    // *** DIAGRAM 6364 OF 15495 ***
    // Wavefunction(s) for diagram number 6364
    // (none)
    // Amplitude(s) for diagram number 6364
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[547], w_fp[497], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[257] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[307] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[310] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6365 OF 15495 ***
    // Wavefunction(s) for diagram number 6365
    // (none)
    // Amplitude(s) for diagram number 6365
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[29], w_fp[156], w_fp[547], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[257] += amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[341] -= amp_sv[0];
    jamp_sv[355] += amp_sv[0];

    // *** DIAGRAM 6366 OF 15495 ***
    // Wavefunction(s) for diagram number 6366
    // (none)
    // Amplitude(s) for diagram number 6366
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[157], w_fp[158], w_fp[547], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[291] += amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[307] -= amp_sv[0];
    jamp_sv[310] += amp_sv[0];

    // *** DIAGRAM 6367 OF 15495 ***
    // Wavefunction(s) for diagram number 6367
    // (none)
    // Amplitude(s) for diagram number 6367
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[29], w_fp[163], w_fp[530], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[341] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6368 OF 15495 ***
    // Wavefunction(s) for diagram number 6368
    // (none)
    // Amplitude(s) for diagram number 6368
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[77], w_fp[158], w_fp[530], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6369 OF 15495 ***
    // Wavefunction(s) for diagram number 6369
    // (none)
    // Amplitude(s) for diagram number 6369
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[157], w_fp[156], w_fp[454], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[255] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[307] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[310] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[157], w_fp[156], w_fp[575], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[257] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[307] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[310] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[157], w_fp[156], w_fp[589], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[255] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6370 OF 15495 ***
    // Wavefunction(s) for diagram number 6370
    // (none)
    // Amplitude(s) for diagram number 6370
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[495], w_fp[156], w_fp[100], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[307] += amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[349] -= amp_sv[0];
    jamp_sv[352] += amp_sv[0];

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
  diagramgroup638( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 77 );
    retrieveWf( wfs, w_cx, nevt, 94 );
    retrieveWf( wfs, w_cx, nevt, 100 );
    retrieveWf( wfs, w_cx, nevt, 156 );
    retrieveWf( wfs, w_cx, nevt, 157 );
    retrieveWf( wfs, w_cx, nevt, 169 );
    retrieveWf( wfs, w_cx, nevt, 190 );
    retrieveWf( wfs, w_cx, nevt, 193 );
    retrieveWf( wfs, w_cx, nevt, 495 );
    retrieveWf( wfs, w_cx, nevt, 510 );
    retrieveWf( wfs, w_cx, nevt, 523 );
    retrieveWf( wfs, w_cx, nevt, 530 );
    retrieveWf( wfs, w_cx, nevt, 540 );
    retrieveWf( wfs, w_cx, nevt, 547 );
    retrieveWf( wfs, w_cx, nevt, 553 );
#endif
#endif

    // *** DIAGRAM 6371 OF 15495 ***
    // Wavefunction(s) for diagram number 6371
    // (none)
    // Amplitude(s) for diagram number 6371
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[157], w_fp[553], w_fp[100], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[255] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];

    // *** DIAGRAM 6372 OF 15495 ***
    // Wavefunction(s) for diagram number 6372
    // (none)
    // Amplitude(s) for diagram number 6372
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[157], w_fp[156], w_fp[510], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[255] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[307] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[310] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6373 OF 15495 ***
    // Wavefunction(s) for diagram number 6373
    // (none)
    // Amplitude(s) for diagram number 6373
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[495], w_fp[190], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[427] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6374 OF 15495 ***
    // Wavefunction(s) for diagram number 6374
    // (none)
    // Amplitude(s) for diagram number 6374
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[495], w_fp[193], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[469] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6375 OF 15495 ***
    // Wavefunction(s) for diagram number 6375
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[169], w_fp[530], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[535] );
    // Amplitude(s) for diagram number 6375
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[94], w_fp[535], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[437] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6376 OF 15495 ***
    // Wavefunction(s) for diagram number 6376
    // (none)
    // Amplitude(s) for diagram number 6376
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[77], w_fp[535], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6377 OF 15495 ***
    // Wavefunction(s) for diagram number 6377
    // (none)
    // Amplitude(s) for diagram number 6377
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[523], w_fp[540], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[469] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[475] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6378 OF 15495 ***
    // Wavefunction(s) for diagram number 6378
    // (none)
    // Amplitude(s) for diagram number 6378
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[77], w_fp[169], w_fp[523], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[375] += amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[435] -= amp_sv[0];

    // *** DIAGRAM 6379 OF 15495 ***
    // Wavefunction(s) for diagram number 6379
    // (none)
    // Amplitude(s) for diagram number 6379
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[157], w_fp[193], w_fp[523], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[461] += amp_sv[0];
    jamp_sv[469] -= amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[475] -= amp_sv[0];

    // *** DIAGRAM 6380 OF 15495 ***
    // Wavefunction(s) for diagram number 6380
    // (none)
    // Amplitude(s) for diagram number 6380
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[547], w_fp[540], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[427] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[437] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[475] += cxtype( 0, 1 ) * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 535 );
#endif
#endif
  }

  //--------------------------------------------------------------------------


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup639( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 77 );
    retrieveWf( wfs, w_cx, nevt, 94 );
    retrieveWf( wfs, w_cx, nevt, 102 );
    retrieveWf( wfs, w_cx, nevt, 157 );
    retrieveWf( wfs, w_cx, nevt, 169 );
    retrieveWf( wfs, w_cx, nevt, 190 );
    retrieveWf( wfs, w_cx, nevt, 193 );
    retrieveWf( wfs, w_cx, nevt, 225 );
    retrieveWf( wfs, w_cx, nevt, 226 );
    retrieveWf( wfs, w_cx, nevt, 486 );
    retrieveWf( wfs, w_cx, nevt, 495 );
    retrieveWf( wfs, w_cx, nevt, 520 );
    retrieveWf( wfs, w_cx, nevt, 530 );
    retrieveWf( wfs, w_cx, nevt, 535 );
    retrieveWf( wfs, w_cx, nevt, 547 );
    retrieveWf( wfs, w_cx, nevt, 556 );
    retrieveWf( wfs, w_cx, nevt, 576 );
#endif
#endif

    // *** DIAGRAM 6381 OF 15495 ***
    // Wavefunction(s) for diagram number 6381
    // (none)
    // Amplitude(s) for diagram number 6381
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[94], w_fp[169], w_fp[547], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[377] += amp_sv[0];
    jamp_sv[437] -= amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    jamp_sv[475] += amp_sv[0];

    // *** DIAGRAM 6382 OF 15495 ***
    // Wavefunction(s) for diagram number 6382
    // (none)
    // Amplitude(s) for diagram number 6382
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[157], w_fp[190], w_fp[547], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[411] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[427] -= amp_sv[0];
    jamp_sv[430] += amp_sv[0];

    // *** DIAGRAM 6383 OF 15495 ***
    // Wavefunction(s) for diagram number 6383
    // (none)
    // Amplitude(s) for diagram number 6383
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[94], w_fp[193], w_fp[530], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[461] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[475] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6384 OF 15495 ***
    // Wavefunction(s) for diagram number 6384
    // (none)
    // Amplitude(s) for diagram number 6384
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[77], w_fp[190], w_fp[530], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[411] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6385 OF 15495 ***
    // Wavefunction(s) for diagram number 6385
    // (none)
    // Amplitude(s) for diagram number 6385
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[157], w_fp[169], w_fp[520], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[427] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[437] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[469] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[157], w_fp[169], w_fp[486], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[427] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[437] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[475] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[157], w_fp[169], w_fp[576], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[469] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[475] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6386 OF 15495 ***
    // Wavefunction(s) for diagram number 6386
    // (none)
    // Amplitude(s) for diagram number 6386
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[495], w_fp[169], w_fp[102], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[427] += amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[469] -= amp_sv[0];
    jamp_sv[472] += amp_sv[0];

    // *** DIAGRAM 6387 OF 15495 ***
    // Wavefunction(s) for diagram number 6387
    // (none)
    // Amplitude(s) for diagram number 6387
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[157], w_fp[535], w_fp[102], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[375] += amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[437] += amp_sv[0];

    // *** DIAGRAM 6388 OF 15495 ***
    // Wavefunction(s) for diagram number 6388
    // (none)
    // Amplitude(s) for diagram number 6388
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[157], w_fp[169], w_fp[556], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[427] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[437] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[469] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6389 OF 15495 ***
    // Wavefunction(s) for diagram number 6389
    // (none)
    // Amplitude(s) for diagram number 6389
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[495], w_fp[225], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[664] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6390 OF 15495 ***
    // Wavefunction(s) for diagram number 6390
    // (none)
    // Amplitude(s) for diagram number 6390
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[495], w_fp[226], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[685] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[688] -= cxtype( 0, 1 ) * amp_sv[0];

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
  diagramgroup640( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 29 );
    retrieveWf( wfs, w_cx, nevt, 94 );
    retrieveWf( wfs, w_cx, nevt, 157 );
    retrieveWf( wfs, w_cx, nevt, 215 );
    retrieveWf( wfs, w_cx, nevt, 225 );
    retrieveWf( wfs, w_cx, nevt, 226 );
    retrieveWf( wfs, w_cx, nevt, 474 );
    retrieveWf( wfs, w_cx, nevt, 523 );
    retrieveWf( wfs, w_cx, nevt, 530 );
    retrieveWf( wfs, w_cx, nevt, 544 );
#endif
#endif

    // *** DIAGRAM 6391 OF 15495 ***
    // Wavefunction(s) for diagram number 6391
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[215], w_fp[530], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[435] );
    // Amplitude(s) for diagram number 6391
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[94], w_fp[435], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[623] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6392 OF 15495 ***
    // Wavefunction(s) for diagram number 6392
    // (none)
    // Amplitude(s) for diagram number 6392
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[29], w_fp[435], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[621] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6393 OF 15495 ***
    // Wavefunction(s) for diagram number 6393
    // (none)
    // Amplitude(s) for diagram number 6393
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[523], w_fp[544], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[621] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[653] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[688] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[691] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6394 OF 15495 ***
    // Wavefunction(s) for diagram number 6394
    // (none)
    // Amplitude(s) for diagram number 6394
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[29], w_fp[215], w_fp[523], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[621] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[699] -= amp_sv[0];

    // *** DIAGRAM 6395 OF 15495 ***
    // Wavefunction(s) for diagram number 6395
    // (none)
    // Amplitude(s) for diagram number 6395
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[157], w_fp[226], w_fp[523], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[677] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[691] -= amp_sv[0];

    // *** DIAGRAM 6396 OF 15495 ***
    // Wavefunction(s) for diagram number 6396
    // (none)
    // Amplitude(s) for diagram number 6396
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[474], w_fp[544], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[623] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[653] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[664] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[691] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6397 OF 15495 ***
    // Wavefunction(s) for diagram number 6397
    // (none)
    // Amplitude(s) for diagram number 6397
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[94], w_fp[215], w_fp[474], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[623] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[701] -= amp_sv[0];

    // *** DIAGRAM 6398 OF 15495 ***
    // Wavefunction(s) for diagram number 6398
    // (none)
    // Amplitude(s) for diagram number 6398
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[157], w_fp[225], w_fp[474], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[653] += amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];

    // *** DIAGRAM 6399 OF 15495 ***
    // Wavefunction(s) for diagram number 6399
    // (none)
    // Amplitude(s) for diagram number 6399
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[94], w_fp[226], w_fp[530], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[691] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6400 OF 15495 ***
    // Wavefunction(s) for diagram number 6400
    // (none)
    // Amplitude(s) for diagram number 6400
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[29], w_fp[225], w_fp[530], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[653] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 435 );
#endif
#endif
  }

  //--------------------------------------------------------------------------

}
