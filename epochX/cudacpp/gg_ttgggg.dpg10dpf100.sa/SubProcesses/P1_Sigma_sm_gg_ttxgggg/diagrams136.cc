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
  diagramgroup1351( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 1 );
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 128 );
    retrieveWf( wfs, w_cx, nevt, 131 );
    retrieveWf( wfs, w_cx, nevt, 216 );
    retrieveWf( wfs, w_cx, nevt, 238 );
    retrieveWf( wfs, w_cx, nevt, 239 );
    retrieveWf( wfs, w_cx, nevt, 272 );
    retrieveWf( wfs, w_cx, nevt, 355 );
    retrieveWf( wfs, w_cx, nevt, 361 );
    retrieveWf( wfs, w_cx, nevt, 497 );
    retrieveWf( wfs, w_cx, nevt, 588 );
    retrieveWf( wfs, w_cx, nevt, 600 );
    retrieveWf( wfs, w_cx, nevt, 674 );
    retrieveWf( wfs, w_cx, nevt, 697 );
#endif
#endif

    // *** DIAGRAM 13501 OF 15495 ***
    // Wavefunction(s) for diagram number 13501
    // (none)
    // Amplitude(s) for diagram number 13501
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[216], w_fp[497], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[91] += amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[697] += amp_sv[0];

    // *** DIAGRAM 13502 OF 15495 ***
    // Wavefunction(s) for diagram number 13502
    // (none)
    // Amplitude(s) for diagram number 13502
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[674], w_fp[128], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13503 OF 15495 ***
    // Wavefunction(s) for diagram number 13503
    // (none)
    // Amplitude(s) for diagram number 13503
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[355], w_fp[697], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[580] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13504 OF 15495 ***
    // Wavefunction(s) for diagram number 13504
    // (none)
    // Amplitude(s) for diagram number 13504
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[674], w_fp[2], w_fp[131], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[450] += amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];

    // *** DIAGRAM 13505 OF 15495 ***
    // Wavefunction(s) for diagram number 13505
    // (none)
    // Amplitude(s) for diagram number 13505
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[355], w_fp[2], w_fp[600], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13506 OF 15495 ***
    // Wavefunction(s) for diagram number 13506
    // (none)
    // Amplitude(s) for diagram number 13506
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[361], w_fp[239], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[235] += amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[382] -= amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[361], w_fp[239], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[53] += amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[586] += amp_sv[0];
    jamp_sv[706] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[361], w_fp[239], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
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

    // *** DIAGRAM 13507 OF 15495 ***
    // Wavefunction(s) for diagram number 13507
    // (none)
    // Amplitude(s) for diagram number 13507
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[239], w_fp[5], w_fp[238], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[235] += amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[382] -= amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];

    // *** DIAGRAM 13508 OF 15495 ***
    // Wavefunction(s) for diagram number 13508
    // (none)
    // Amplitude(s) for diagram number 13508
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[361], w_fp[5], w_fp[272], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[53] += amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[586] += amp_sv[0];
    jamp_sv[706] -= amp_sv[0];

    // *** DIAGRAM 13509 OF 15495 ***
    // Wavefunction(s) for diagram number 13509
    // (none)
    // Amplitude(s) for diagram number 13509
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[216], w_fp[2], w_fp[238], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13510 OF 15495 ***
    // Wavefunction(s) for diagram number 13510
    // (none)
    // Amplitude(s) for diagram number 13510
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[588], w_fp[2], w_fp[361], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[211] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[703] += amp_sv[0];

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
  diagramgroup1352( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 1 );
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 65 );
    retrieveWf( wfs, w_cx, nevt, 71 );
    retrieveWf( wfs, w_cx, nevt, 128 );
    retrieveWf( wfs, w_cx, nevt, 131 );
    retrieveWf( wfs, w_cx, nevt, 132 );
    retrieveWf( wfs, w_cx, nevt, 216 );
    retrieveWf( wfs, w_cx, nevt, 239 );
    retrieveWf( wfs, w_cx, nevt, 249 );
    retrieveWf( wfs, w_cx, nevt, 257 );
    retrieveWf( wfs, w_cx, nevt, 272 );
    retrieveWf( wfs, w_cx, nevt, 354 );
    retrieveWf( wfs, w_cx, nevt, 355 );
    retrieveWf( wfs, w_cx, nevt, 451 );
    retrieveWf( wfs, w_cx, nevt, 549 );
    retrieveWf( wfs, w_cx, nevt, 571 );
    retrieveWf( wfs, w_cx, nevt, 588 );
    retrieveWf( wfs, w_cx, nevt, 600 );
    retrieveWf( wfs, w_cx, nevt, 697 );
#endif
#endif

    // *** DIAGRAM 13511 OF 15495 ***
    // Wavefunction(s) for diagram number 13511
    // (none)
    // Amplitude(s) for diagram number 13511
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[239], w_fp[131], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[714] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[239], w_fp[131], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[239], w_fp[131], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];

    // *** DIAGRAM 13512 OF 15495 ***
    // Wavefunction(s) for diagram number 13512
    // (none)
    // Amplitude(s) for diagram number 13512
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[131], w_fp[272], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];

    // *** DIAGRAM 13513 OF 15495 ***
    // Wavefunction(s) for diagram number 13513
    // (none)
    // Amplitude(s) for diagram number 13513
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[239], w_fp[600], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];

    // *** DIAGRAM 13514 OF 15495 ***
    // Wavefunction(s) for diagram number 13514
    // (none)
    // Amplitude(s) for diagram number 13514
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[216], w_fp[697], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[577] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[697] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13515 OF 15495 ***
    // Wavefunction(s) for diagram number 13515
    // (none)
    // Amplitude(s) for diagram number 13515
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[588], w_fp[128], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13516 OF 15495 ***
    // Wavefunction(s) for diagram number 13516
    // (none)
    // Amplitude(s) for diagram number 13516
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[71], w_fp[239], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[235] += amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[382] -= amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[132], w_fp[239], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[235] += amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[382] -= amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[577] += amp_sv[0];
    jamp_sv[697] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[65], w_fp[239], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[363] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[577] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[697] -= amp_sv[0];
    jamp_sv[703] += amp_sv[0];

    // *** DIAGRAM 13517 OF 15495 ***
    // Wavefunction(s) for diagram number 13517
    // (none)
    // Amplitude(s) for diagram number 13517
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[216], w_fp[2], w_fp[71], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[216], w_fp[2], w_fp[132], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[91] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[697] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[216], w_fp[2], w_fp[65], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[697] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13518 OF 15495 ***
    // Wavefunction(s) for diagram number 13518
    // (none)
    // Amplitude(s) for diagram number 13518
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[355], w_fp[2], w_fp[549], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[355], w_fp[2], w_fp[451], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[355], w_fp[2], w_fp[571], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13519 OF 15495 ***
    // Wavefunction(s) for diagram number 13519
    // (none)
    // Amplitude(s) for diagram number 13519
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[549], w_fp[1], w_fp[239], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[1], w_fp[239], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[94] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[173] -= amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[211] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[376] -= amp_sv[0];
    jamp_sv[382] += amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[580] -= amp_sv[0];
    jamp_sv[700] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[571], w_fp[1], w_fp[239], COUPs[0], 1.0, &amp_fp[0] );
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

    // *** DIAGRAM 13520 OF 15495 ***
    // Wavefunction(s) for diagram number 13520
    // (none)
    // Amplitude(s) for diagram number 13520
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[257], w_fp[239], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[184] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[716] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[249], w_fp[239], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[53] += amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[586] += amp_sv[0];
    jamp_sv[706] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[354], w_fp[239], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[53] += amp_sv[0];
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    jamp_sv[586] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[706] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];

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
  diagramgroup1353( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 1 );
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 23 );
    retrieveWf( wfs, w_cx, nevt, 24 );
    retrieveWf( wfs, w_cx, nevt, 25 );
    retrieveWf( wfs, w_cx, nevt, 28 );
    retrieveWf( wfs, w_cx, nevt, 38 );
    retrieveWf( wfs, w_cx, nevt, 57 );
    retrieveWf( wfs, w_cx, nevt, 58 );
    retrieveWf( wfs, w_cx, nevt, 60 );
    retrieveWf( wfs, w_cx, nevt, 95 );
    retrieveWf( wfs, w_cx, nevt, 118 );
    retrieveWf( wfs, w_cx, nevt, 196 );
    retrieveWf( wfs, w_cx, nevt, 229 );
    retrieveWf( wfs, w_cx, nevt, 239 );
    retrieveWf( wfs, w_cx, nevt, 345 );
    retrieveWf( wfs, w_cx, nevt, 355 );
    retrieveWf( wfs, w_cx, nevt, 362 );
    retrieveWf( wfs, w_cx, nevt, 434 );
    retrieveWf( wfs, w_cx, nevt, 437 );
    retrieveWf( wfs, w_cx, nevt, 452 );
    retrieveWf( wfs, w_cx, nevt, 488 );
    retrieveWf( wfs, w_cx, nevt, 536 );
    retrieveWf( wfs, w_cx, nevt, 543 );
    retrieveWf( wfs, w_cx, nevt, 651 );
    retrieveWf( wfs, w_cx, nevt, 652 );
    retrieveWf( wfs, w_cx, nevt, 653 );
    retrieveWf( wfs, w_cx, nevt, 660 );
    retrieveWf( wfs, w_cx, nevt, 747 );
#endif
#endif

    // *** DIAGRAM 13521 OF 15495 ***
    // Wavefunction(s) for diagram number 13521
    // (none)
    // Amplitude(s) for diagram number 13521
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[536], w_fp[2], w_fp[345], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[184] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[596] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[716] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[536], w_fp[2], w_fp[95], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[190] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[596] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[536], w_fp[2], w_fp[434], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[184] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[716] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13522 OF 15495 ***
    // Wavefunction(s) for diagram number 13522
    // (none)
    // Amplitude(s) for diagram number 13522
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[536], w_fp[229], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[452] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[536], w_fp[28], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[536], w_fp[60], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    jamp_sv[716] -= amp_sv[0];

    // *** DIAGRAM 13523 OF 15495 ***
    // Wavefunction(s) for diagram number 13523
    // (none)
    // Amplitude(s) for diagram number 13523
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[355], w_fp[2], w_fp[452], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[355], w_fp[2], w_fp[488], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[690] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[355], w_fp[2], w_fp[437], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[690] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13524 OF 15495 ***
    // Wavefunction(s) for diagram number 13524
    // (none)
    // Amplitude(s) for diagram number 13524
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[452], w_fp[1], w_fp[239], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[1], w_fp[239], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[70] += amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[437], w_fp[1], w_fp[239], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] += amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[184] += amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[714] -= amp_sv[0];

    // *** DIAGRAM 13525 OF 15495 ***
    // Wavefunction(s) for diagram number 13525
    // (none)
    // Amplitude(s) for diagram number 13525
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[355], w_fp[229], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[450] += amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[355], w_fp[28], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[355], w_fp[60], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[714] -= amp_sv[0];

    // *** DIAGRAM 13526 OF 15495 ***
    // Wavefunction(s) for diagram number 13526
    // (none)
    // Amplitude(s) for diagram number 13526
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[345], w_fp[239], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[184] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[716] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[95], w_fp[239], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[434], w_fp[239], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];

    // *** DIAGRAM 13527 OF 15495 ***
    // Wavefunction(s) for diagram number 13527
    // (none)
    // Amplitude(s) for diagram number 13527
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[2], w_fp[23], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[184] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[716] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[2], w_fp[24], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[184] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[714] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[2], w_fp[25], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[714] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[2], w_fp[38], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[2], w_fp[57], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[2], w_fp[58], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[15] += amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[2], w_fp[653], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[2], w_fp[652], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[2], w_fp[651], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[15] += amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    jamp_sv[714] += amp_sv[0];
    jamp_sv[716] -= amp_sv[0];

    // *** DIAGRAM 13528 OF 15495 ***
    // Wavefunction(s) for diagram number 13528
    // (none)
    // Amplitude(s) for diagram number 13528
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[362], w_fp[7], w_fp[747], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13529 OF 15495 ***
    // Wavefunction(s) for diagram number 13529
    // (none)
    // Amplitude(s) for diagram number 13529
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[543], w_fp[2], w_fp[362], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] += amp_sv[0];
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[538] += amp_sv[0];

    // *** DIAGRAM 13530 OF 15495 ***
    // Wavefunction(s) for diagram number 13530
    // (none)
    // Amplitude(s) for diagram number 13530
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[660], w_fp[118], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];

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
  diagramgroup1354( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 88 );
    retrieveWf( wfs, w_cx, nevt, 118 );
    retrieveWf( wfs, w_cx, nevt, 154 );
    retrieveWf( wfs, w_cx, nevt, 171 );
    retrieveWf( wfs, w_cx, nevt, 246 );
    retrieveWf( wfs, w_cx, nevt, 256 );
    retrieveWf( wfs, w_cx, nevt, 367 );
    retrieveWf( wfs, w_cx, nevt, 368 );
    retrieveWf( wfs, w_cx, nevt, 369 );
    retrieveWf( wfs, w_cx, nevt, 527 );
    retrieveWf( wfs, w_cx, nevt, 543 );
    retrieveWf( wfs, w_cx, nevt, 576 );
    retrieveWf( wfs, w_cx, nevt, 604 );
    retrieveWf( wfs, w_cx, nevt, 606 );
    retrieveWf( wfs, w_cx, nevt, 660 );
    retrieveWf( wfs, w_cx, nevt, 747 );
#endif
#endif

    // *** DIAGRAM 13531 OF 15495 ***
    // Wavefunction(s) for diagram number 13531
    // (none)
    // Amplitude(s) for diagram number 13531
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[543], w_fp[118], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13532 OF 15495 ***
    // Wavefunction(s) for diagram number 13532
    // (none)
    // Amplitude(s) for diagram number 13532
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[660], w_fp[2], w_fp[88], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[332] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];

    // *** DIAGRAM 13533 OF 15495 ***
    // Wavefunction(s) for diagram number 13533
    // (none)
    // Amplitude(s) for diagram number 13533
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[88], w_fp[747], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13534 OF 15495 ***
    // Wavefunction(s) for diagram number 13534
    // (none)
    // Amplitude(s) for diagram number 13534
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[527], w_fp[2], w_fp[367], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[527], w_fp[2], w_fp[368], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[527], w_fp[2], w_fp[369], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13535 OF 15495 ***
    // Wavefunction(s) for diagram number 13535
    // (none)
    // Amplitude(s) for diagram number 13535
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[256], w_fp[246], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] += amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[316] -= amp_sv[0];
    jamp_sv[532] += amp_sv[0];

    // *** DIAGRAM 13536 OF 15495 ***
    // Wavefunction(s) for diagram number 13536
    // (none)
    // Amplitude(s) for diagram number 13536
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[256], w_fp[2], w_fp[606], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13537 OF 15495 ***
    // Wavefunction(s) for diagram number 13537
    // (none)
    // Amplitude(s) for diagram number 13537
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[576], w_fp[1], w_fp[154], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[157] += amp_sv[0];
    jamp_sv[199] -= amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[634] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[708] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[576], w_fp[1], w_fp[154], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[316] += amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[634] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[576], w_fp[1], w_fp[154], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[316] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[610] += amp_sv[0];
    jamp_sv[620] -= amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];

    // *** DIAGRAM 13538 OF 15495 ***
    // Wavefunction(s) for diagram number 13538
    // (none)
    // Amplitude(s) for diagram number 13538
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[154], w_fp[7], w_fp[604], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[157] += amp_sv[0];
    jamp_sv[199] -= amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[634] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[708] -= amp_sv[0];

    // *** DIAGRAM 13539 OF 15495 ***
    // Wavefunction(s) for diagram number 13539
    // (none)
    // Amplitude(s) for diagram number 13539
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[154], w_fp[606], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[316] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[610] += amp_sv[0];
    jamp_sv[620] -= amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];

    // *** DIAGRAM 13540 OF 15495 ***
    // Wavefunction(s) for diagram number 13540
    // (none)
    // Amplitude(s) for diagram number 13540
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[171], w_fp[2], w_fp[604], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[123] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[313] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] -= cxtype( 0, 1 ) * amp_sv[0];

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
  diagramgroup1355( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 1 );
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 64 );
    retrieveWf( wfs, w_cx, nevt, 88 );
    retrieveWf( wfs, w_cx, nevt, 92 );
    retrieveWf( wfs, w_cx, nevt, 118 );
    retrieveWf( wfs, w_cx, nevt, 154 );
    retrieveWf( wfs, w_cx, nevt, 171 );
    retrieveWf( wfs, w_cx, nevt, 246 );
    retrieveWf( wfs, w_cx, nevt, 256 );
    retrieveWf( wfs, w_cx, nevt, 362 );
    retrieveWf( wfs, w_cx, nevt, 450 );
    retrieveWf( wfs, w_cx, nevt, 609 );
    retrieveWf( wfs, w_cx, nevt, 658 );
    retrieveWf( wfs, w_cx, nevt, 749 );
#endif
#endif

    // *** DIAGRAM 13541 OF 15495 ***
    // Wavefunction(s) for diagram number 13541
    // (none)
    // Amplitude(s) for diagram number 13541
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[171], w_fp[246], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];

    // *** DIAGRAM 13542 OF 15495 ***
    // Wavefunction(s) for diagram number 13542
    // (none)
    // Amplitude(s) for diagram number 13542
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[658], w_fp[118], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[330] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13543 OF 15495 ***
    // Wavefunction(s) for diagram number 13543
    // (none)
    // Amplitude(s) for diagram number 13543
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[256], w_fp[92], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[316] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13544 OF 15495 ***
    // Wavefunction(s) for diagram number 13544
    // (none)
    // Amplitude(s) for diagram number 13544
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[658], w_fp[2], w_fp[88], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[330] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];

    // *** DIAGRAM 13545 OF 15495 ***
    // Wavefunction(s) for diagram number 13545
    // (none)
    // Amplitude(s) for diagram number 13545
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[256], w_fp[2], w_fp[609], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[330] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13546 OF 15495 ***
    // Wavefunction(s) for diagram number 13546
    // (none)
    // Amplitude(s) for diagram number 13546
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[362], w_fp[154], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[319] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[604] += amp_sv[0];
    jamp_sv[610] += amp_sv[0];
    jamp_sv[620] -= amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[362], w_fp[154], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[362], w_fp[154], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[157] += amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[199] -= amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[535] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[601] += amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];

    // *** DIAGRAM 13547 OF 15495 ***
    // Wavefunction(s) for diagram number 13547
    // (none)
    // Amplitude(s) for diagram number 13547
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[154], w_fp[7], w_fp[64], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[319] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[604] += amp_sv[0];
    jamp_sv[610] += amp_sv[0];
    jamp_sv[620] -= amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];

    // *** DIAGRAM 13548 OF 15495 ***
    // Wavefunction(s) for diagram number 13548
    // (none)
    // Amplitude(s) for diagram number 13548
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[362], w_fp[7], w_fp[749], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];

    // *** DIAGRAM 13549 OF 15495 ***
    // Wavefunction(s) for diagram number 13549
    // (none)
    // Amplitude(s) for diagram number 13549
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[171], w_fp[2], w_fp[64], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[37] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13550 OF 15495 ***
    // Wavefunction(s) for diagram number 13550
    // (none)
    // Amplitude(s) for diagram number 13550
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[450], w_fp[2], w_fp[362], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[157] += amp_sv[0];
    jamp_sv[199] -= amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[535] += amp_sv[0];

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
  diagramgroup1356( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 1 );
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 20 );
    retrieveWf( wfs, w_cx, nevt, 21 );
    retrieveWf( wfs, w_cx, nevt, 22 );
    retrieveWf( wfs, w_cx, nevt, 88 );
    retrieveWf( wfs, w_cx, nevt, 92 );
    retrieveWf( wfs, w_cx, nevt, 118 );
    retrieveWf( wfs, w_cx, nevt, 154 );
    retrieveWf( wfs, w_cx, nevt, 171 );
    retrieveWf( wfs, w_cx, nevt, 256 );
    retrieveWf( wfs, w_cx, nevt, 367 );
    retrieveWf( wfs, w_cx, nevt, 368 );
    retrieveWf( wfs, w_cx, nevt, 369 );
    retrieveWf( wfs, w_cx, nevt, 450 );
    retrieveWf( wfs, w_cx, nevt, 480 );
    retrieveWf( wfs, w_cx, nevt, 550 );
    retrieveWf( wfs, w_cx, nevt, 558 );
    retrieveWf( wfs, w_cx, nevt, 609 );
    retrieveWf( wfs, w_cx, nevt, 749 );
#endif
#endif

    // *** DIAGRAM 13551 OF 15495 ***
    // Wavefunction(s) for diagram number 13551
    // (none)
    // Amplitude(s) for diagram number 13551
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[154], w_fp[88], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    jamp_sv[708] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[154], w_fp[88], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[154], w_fp[88], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];

    // *** DIAGRAM 13552 OF 15495 ***
    // Wavefunction(s) for diagram number 13552
    // (none)
    // Amplitude(s) for diagram number 13552
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[88], w_fp[749], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];

    // *** DIAGRAM 13553 OF 15495 ***
    // Wavefunction(s) for diagram number 13553
    // (none)
    // Amplitude(s) for diagram number 13553
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[154], w_fp[609], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];

    // *** DIAGRAM 13554 OF 15495 ***
    // Wavefunction(s) for diagram number 13554
    // (none)
    // Amplitude(s) for diagram number 13554
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[171], w_fp[92], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[313] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13555 OF 15495 ***
    // Wavefunction(s) for diagram number 13555
    // (none)
    // Amplitude(s) for diagram number 13555
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[450], w_fp[118], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[319] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13556 OF 15495 ***
    // Wavefunction(s) for diagram number 13556
    // (none)
    // Amplitude(s) for diagram number 13556
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[22], w_fp[154], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[319] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[604] += amp_sv[0];
    jamp_sv[610] += amp_sv[0];
    jamp_sv[620] -= amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[154], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[313] += amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[610] += amp_sv[0];
    jamp_sv[620] -= amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[628] += amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[20], w_fp[154], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[13] += amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[313] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[535] += amp_sv[0];
    jamp_sv[601] += amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[628] += amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[668] += amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    jamp_sv[710] -= amp_sv[0];

    // *** DIAGRAM 13557 OF 15495 ***
    // Wavefunction(s) for diagram number 13557
    // (none)
    // Amplitude(s) for diagram number 13557
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[171], w_fp[2], w_fp[22], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[37] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[171], w_fp[2], w_fp[21], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[123] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[313] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[171], w_fp[2], w_fp[20], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[123] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[313] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13558 OF 15495 ***
    // Wavefunction(s) for diagram number 13558
    // (none)
    // Amplitude(s) for diagram number 13558
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[256], w_fp[2], w_fp[550], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[330] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[256], w_fp[2], w_fp[480], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[330] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[256], w_fp[2], w_fp[558], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13559 OF 15495 ***
    // Wavefunction(s) for diagram number 13559
    // (none)
    // Amplitude(s) for diagram number 13559
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[550], w_fp[1], w_fp[154], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[480], w_fp[1], w_fp[154], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[106] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[157] += amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[199] -= amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[217] += amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[316] -= amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[558], w_fp[1], w_fp[154], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] += amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[157] += amp_sv[0];
    jamp_sv[199] -= amp_sv[0];
    jamp_sv[217] += amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[316] -= amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[708] -= amp_sv[0];

    // *** DIAGRAM 13560 OF 15495 ***
    // Wavefunction(s) for diagram number 13560
    // (none)
    // Amplitude(s) for diagram number 13560
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[367], w_fp[154], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[13] += amp_sv[0];
    jamp_sv[19] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[40] += amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[160] += amp_sv[0];
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[668] += amp_sv[0];
    jamp_sv[710] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[368], w_fp[154], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[40] += amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[369], w_fp[154], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];

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
  diagramgroup1357( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 98 );
    retrieveWf( wfs, w_cx, nevt, 104 );
    retrieveWf( wfs, w_cx, nevt, 154 );
    retrieveWf( wfs, w_cx, nevt, 245 );
    retrieveWf( wfs, w_cx, nevt, 256 );
    retrieveWf( wfs, w_cx, nevt, 363 );
    retrieveWf( wfs, w_cx, nevt, 371 );
    retrieveWf( wfs, w_cx, nevt, 372 );
    retrieveWf( wfs, w_cx, nevt, 373 );
    retrieveWf( wfs, w_cx, nevt, 518 );
    retrieveWf( wfs, w_cx, nevt, 527 );
    retrieveWf( wfs, w_cx, nevt, 538 );
    retrieveWf( wfs, w_cx, nevt, 613 );
    retrieveWf( wfs, w_cx, nevt, 660 );
    retrieveWf( wfs, w_cx, nevt, 747 );
#endif
#endif

    // *** DIAGRAM 13561 OF 15495 ***
    // Wavefunction(s) for diagram number 13561
    // (none)
    // Amplitude(s) for diagram number 13561
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[363], w_fp[6], w_fp[747], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13562 OF 15495 ***
    // Wavefunction(s) for diagram number 13562
    // (none)
    // Amplitude(s) for diagram number 13562
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[538], w_fp[2], w_fp[363], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[166] += amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    jamp_sv[346] -= amp_sv[0];
    jamp_sv[658] += amp_sv[0];

    // *** DIAGRAM 13563 OF 15495 ***
    // Wavefunction(s) for diagram number 13563
    // (none)
    // Amplitude(s) for diagram number 13563
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[660], w_fp[98], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13564 OF 15495 ***
    // Wavefunction(s) for diagram number 13564
    // (none)
    // Amplitude(s) for diagram number 13564
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[538], w_fp[98], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[346] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13565 OF 15495 ***
    // Wavefunction(s) for diagram number 13565
    // (none)
    // Amplitude(s) for diagram number 13565
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[660], w_fp[2], w_fp[104], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[356] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[668] -= amp_sv[0];

    // *** DIAGRAM 13566 OF 15495 ***
    // Wavefunction(s) for diagram number 13566
    // (none)
    // Amplitude(s) for diagram number 13566
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[104], w_fp[747], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13567 OF 15495 ***
    // Wavefunction(s) for diagram number 13567
    // (none)
    // Amplitude(s) for diagram number 13567
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[527], w_fp[2], w_fp[371], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[527], w_fp[2], w_fp[372], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[527], w_fp[2], w_fp[373], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13568 OF 15495 ***
    // Wavefunction(s) for diagram number 13568
    // (none)
    // Amplitude(s) for diagram number 13568
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[256], w_fp[245], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] += amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[340] -= amp_sv[0];
    jamp_sv[652] += amp_sv[0];

    // *** DIAGRAM 13569 OF 15495 ***
    // Wavefunction(s) for diagram number 13569
    // (none)
    // Amplitude(s) for diagram number 13569
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[256], w_fp[2], w_fp[613], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[340] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[500] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13570 OF 15495 ***
    // Wavefunction(s) for diagram number 13570
    // (none)
    // Amplitude(s) for diagram number 13570
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[518], w_fp[1], w_fp[154], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[43] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[163] += amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    jamp_sv[337] -= amp_sv[0];
    jamp_sv[490] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[514] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[518], w_fp[1], w_fp[154], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[43] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[196] += amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[337] -= amp_sv[0];
    jamp_sv[340] += amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[514] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    jamp_sv[652] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[518], w_fp[1], w_fp[154], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[163] -= amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[196] += amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[340] += amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[500] -= amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[652] -= amp_sv[0];

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
  diagramgroup1358( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 1 );
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 61 );
    retrieveWf( wfs, w_cx, nevt, 98 );
    retrieveWf( wfs, w_cx, nevt, 104 );
    retrieveWf( wfs, w_cx, nevt, 110 );
    retrieveWf( wfs, w_cx, nevt, 154 );
    retrieveWf( wfs, w_cx, nevt, 170 );
    retrieveWf( wfs, w_cx, nevt, 245 );
    retrieveWf( wfs, w_cx, nevt, 256 );
    retrieveWf( wfs, w_cx, nevt, 363 );
    retrieveWf( wfs, w_cx, nevt, 611 );
    retrieveWf( wfs, w_cx, nevt, 613 );
    retrieveWf( wfs, w_cx, nevt, 616 );
    retrieveWf( wfs, w_cx, nevt, 658 );
#endif
#endif

    // *** DIAGRAM 13571 OF 15495 ***
    // Wavefunction(s) for diagram number 13571
    // (none)
    // Amplitude(s) for diagram number 13571
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[154], w_fp[6], w_fp[611], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[43] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[163] += amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    jamp_sv[337] -= amp_sv[0];
    jamp_sv[490] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[514] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];

    // *** DIAGRAM 13572 OF 15495 ***
    // Wavefunction(s) for diagram number 13572
    // (none)
    // Amplitude(s) for diagram number 13572
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[154], w_fp[613], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[163] -= amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[196] += amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[340] += amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[500] -= amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[652] -= amp_sv[0];

    // *** DIAGRAM 13573 OF 15495 ***
    // Wavefunction(s) for diagram number 13573
    // (none)
    // Amplitude(s) for diagram number 13573
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[170], w_fp[2], w_fp[611], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[43] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[163] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13574 OF 15495 ***
    // Wavefunction(s) for diagram number 13574
    // (none)
    // Amplitude(s) for diagram number 13574
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[170], w_fp[245], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[43] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[337] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];

    // *** DIAGRAM 13575 OF 15495 ***
    // Wavefunction(s) for diagram number 13575
    // (none)
    // Amplitude(s) for diagram number 13575
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[658], w_fp[98], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[354] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13576 OF 15495 ***
    // Wavefunction(s) for diagram number 13576
    // (none)
    // Amplitude(s) for diagram number 13576
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[256], w_fp[110], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[340] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13577 OF 15495 ***
    // Wavefunction(s) for diagram number 13577
    // (none)
    // Amplitude(s) for diagram number 13577
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[658], w_fp[2], w_fp[104], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[354] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[666] -= amp_sv[0];

    // *** DIAGRAM 13578 OF 15495 ***
    // Wavefunction(s) for diagram number 13578
    // (none)
    // Amplitude(s) for diagram number 13578
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[256], w_fp[2], w_fp[616], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[354] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13579 OF 15495 ***
    // Wavefunction(s) for diagram number 13579
    // (none)
    // Amplitude(s) for diagram number 13579
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[363], w_fp[154], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[163] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[343] += amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[484] += amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[500] -= amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[363], w_fp[154], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[346] += amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[363], w_fp[154], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[163] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[343] -= amp_sv[0];
    jamp_sv[346] += amp_sv[0];
    jamp_sv[481] += amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[490] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[655] += amp_sv[0];
    jamp_sv[658] -= amp_sv[0];

    // *** DIAGRAM 13580 OF 15495 ***
    // Wavefunction(s) for diagram number 13580
    // (none)
    // Amplitude(s) for diagram number 13580
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[154], w_fp[6], w_fp[61], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[163] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[343] += amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[484] += amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[500] -= amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];

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
  diagramgroup1359( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 1 );
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 17 );
    retrieveWf( wfs, w_cx, nevt, 18 );
    retrieveWf( wfs, w_cx, nevt, 19 );
    retrieveWf( wfs, w_cx, nevt, 61 );
    retrieveWf( wfs, w_cx, nevt, 98 );
    retrieveWf( wfs, w_cx, nevt, 104 );
    retrieveWf( wfs, w_cx, nevt, 110 );
    retrieveWf( wfs, w_cx, nevt, 154 );
    retrieveWf( wfs, w_cx, nevt, 170 );
    retrieveWf( wfs, w_cx, nevt, 363 );
    retrieveWf( wfs, w_cx, nevt, 535 );
    retrieveWf( wfs, w_cx, nevt, 616 );
    retrieveWf( wfs, w_cx, nevt, 749 );
#endif
#endif

    // *** DIAGRAM 13581 OF 15495 ***
    // Wavefunction(s) for diagram number 13581
    // (none)
    // Amplitude(s) for diagram number 13581
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[363], w_fp[6], w_fp[749], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[346] += amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[658] -= amp_sv[0];

    // *** DIAGRAM 13582 OF 15495 ***
    // Wavefunction(s) for diagram number 13582
    // (none)
    // Amplitude(s) for diagram number 13582
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[170], w_fp[2], w_fp[61], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[43] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[163] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[343] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[655] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13583 OF 15495 ***
    // Wavefunction(s) for diagram number 13583
    // (none)
    // Amplitude(s) for diagram number 13583
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[535], w_fp[2], w_fp[363], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[163] += amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    jamp_sv[343] -= amp_sv[0];
    jamp_sv[655] += amp_sv[0];

    // *** DIAGRAM 13584 OF 15495 ***
    // Wavefunction(s) for diagram number 13584
    // (none)
    // Amplitude(s) for diagram number 13584
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[154], w_fp[104], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[354] -= amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[154], w_fp[104], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[154], w_fp[104], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[136] += amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[354] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[666] -= amp_sv[0];

    // *** DIAGRAM 13585 OF 15495 ***
    // Wavefunction(s) for diagram number 13585
    // (none)
    // Amplitude(s) for diagram number 13585
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[104], w_fp[749], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[668] -= amp_sv[0];

    // *** DIAGRAM 13586 OF 15495 ***
    // Wavefunction(s) for diagram number 13586
    // (none)
    // Amplitude(s) for diagram number 13586
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[154], w_fp[616], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[136] += amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[354] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[666] -= amp_sv[0];

    // *** DIAGRAM 13587 OF 15495 ***
    // Wavefunction(s) for diagram number 13587
    // (none)
    // Amplitude(s) for diagram number 13587
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[170], w_fp[110], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[337] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13588 OF 15495 ***
    // Wavefunction(s) for diagram number 13588
    // (none)
    // Amplitude(s) for diagram number 13588
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[535], w_fp[98], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[343] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[655] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13589 OF 15495 ***
    // Wavefunction(s) for diagram number 13589
    // (none)
    // Amplitude(s) for diagram number 13589
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[19], w_fp[154], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[163] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[343] += amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[484] += amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[500] -= amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[18], w_fp[154], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[163] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[337] += amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[500] -= amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[649] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[17], w_fp[154], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[19] += amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[337] += amp_sv[0];
    jamp_sv[343] -= amp_sv[0];
    jamp_sv[481] += amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    jamp_sv[649] -= amp_sv[0];
    jamp_sv[655] += amp_sv[0];

    // *** DIAGRAM 13590 OF 15495 ***
    // Wavefunction(s) for diagram number 13590
    // (none)
    // Amplitude(s) for diagram number 13590
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[170], w_fp[2], w_fp[19], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[43] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[163] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[343] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[655] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[170], w_fp[2], w_fp[18], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[43] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[163] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[170], w_fp[2], w_fp[17], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[343] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[655] -= cxtype( 0, 1 ) * amp_sv[0];

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
  diagramgroup1360( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 1 );
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 128 );
    retrieveWf( wfs, w_cx, nevt, 130 );
    retrieveWf( wfs, w_cx, nevt, 154 );
    retrieveWf( wfs, w_cx, nevt, 248 );
    retrieveWf( wfs, w_cx, nevt, 256 );
    retrieveWf( wfs, w_cx, nevt, 264 );
    retrieveWf( wfs, w_cx, nevt, 291 );
    retrieveWf( wfs, w_cx, nevt, 361 );
    retrieveWf( wfs, w_cx, nevt, 371 );
    retrieveWf( wfs, w_cx, nevt, 372 );
    retrieveWf( wfs, w_cx, nevt, 373 );
    retrieveWf( wfs, w_cx, nevt, 486 );
    retrieveWf( wfs, w_cx, nevt, 509 );
    retrieveWf( wfs, w_cx, nevt, 527 );
    retrieveWf( wfs, w_cx, nevt, 529 );
    retrieveWf( wfs, w_cx, nevt, 547 );
    retrieveWf( wfs, w_cx, nevt, 660 );
    retrieveWf( wfs, w_cx, nevt, 747 );
#endif
#endif

    // *** DIAGRAM 13591 OF 15495 ***
    // Wavefunction(s) for diagram number 13591
    // (none)
    // Amplitude(s) for diagram number 13591
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[256], w_fp[2], w_fp[486], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[354] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[256], w_fp[2], w_fp[547], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[82] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[340] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[354] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[500] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[256], w_fp[2], w_fp[529], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[340] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[500] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13592 OF 15495 ***
    // Wavefunction(s) for diagram number 13592
    // (none)
    // Amplitude(s) for diagram number 13592
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[486], w_fp[1], w_fp[154], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[136] += amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[354] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[547], w_fp[1], w_fp[154], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[82] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[136] += amp_sv[0];
    jamp_sv[163] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[193] += amp_sv[0];
    jamp_sv[196] -= amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[340] -= amp_sv[0];
    jamp_sv[354] += amp_sv[0];
    jamp_sv[490] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[652] += amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[529], w_fp[1], w_fp[154], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] += amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[163] += amp_sv[0];
    jamp_sv[193] += amp_sv[0];
    jamp_sv[196] -= amp_sv[0];
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    jamp_sv[340] -= amp_sv[0];
    jamp_sv[490] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[652] += amp_sv[0];

    // *** DIAGRAM 13593 OF 15495 ***
    // Wavefunction(s) for diagram number 13593
    // (none)
    // Amplitude(s) for diagram number 13593
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[371], w_fp[154], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[13] += amp_sv[0];
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[19] += amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[166] += amp_sv[0];
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    jamp_sv[356] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    jamp_sv[668] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[372], w_fp[154], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] += amp_sv[0];
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[346] += amp_sv[0];
    jamp_sv[356] -= amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    jamp_sv[668] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[373], w_fp[154], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[346] += amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[658] -= amp_sv[0];

    // *** DIAGRAM 13594 OF 15495 ***
    // Wavefunction(s) for diagram number 13594
    // (none)
    // Amplitude(s) for diagram number 13594
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[361], w_fp[4], w_fp[747], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13595 OF 15495 ***
    // Wavefunction(s) for diagram number 13595
    // (none)
    // Amplitude(s) for diagram number 13595
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[509], w_fp[2], w_fp[361], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[212] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[704] += amp_sv[0];

    // *** DIAGRAM 13596 OF 15495 ***
    // Wavefunction(s) for diagram number 13596
    // (none)
    // Amplitude(s) for diagram number 13596
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[660], w_fp[128], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[590] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13597 OF 15495 ***
    // Wavefunction(s) for diagram number 13597
    // (none)
    // Amplitude(s) for diagram number 13597
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[509], w_fp[128], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[584] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13598 OF 15495 ***
    // Wavefunction(s) for diagram number 13598
    // (none)
    // Amplitude(s) for diagram number 13598
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[660], w_fp[2], w_fp[130], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[332] += amp_sv[0];
    jamp_sv[356] -= amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];

    // *** DIAGRAM 13599 OF 15495 ***
    // Wavefunction(s) for diagram number 13599
    // (none)
    // Amplitude(s) for diagram number 13599
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[130], w_fp[747], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13600 OF 15495 ***
    // Wavefunction(s) for diagram number 13600
    // (none)
    // Amplitude(s) for diagram number 13600
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[527], w_fp[2], w_fp[291], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[527], w_fp[2], w_fp[248], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[527], w_fp[2], w_fp[264], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] -= cxtype( 0, 1 ) * amp_sv[0];

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
