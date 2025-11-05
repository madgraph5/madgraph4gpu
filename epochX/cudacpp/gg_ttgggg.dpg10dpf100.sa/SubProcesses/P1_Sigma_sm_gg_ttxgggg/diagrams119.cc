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
  diagramgroup1181( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 66 );
    retrieveWf( wfs, w_cx, nevt, 168 );
    retrieveWf( wfs, w_cx, nevt, 197 );
    retrieveWf( wfs, w_cx, nevt, 211 );
    retrieveWf( wfs, w_cx, nevt, 214 );
    retrieveWf( wfs, w_cx, nevt, 218 );
    retrieveWf( wfs, w_cx, nevt, 219 );
    retrieveWf( wfs, w_cx, nevt, 221 );
    retrieveWf( wfs, w_cx, nevt, 222 );
    retrieveWf( wfs, w_cx, nevt, 325 );
    retrieveWf( wfs, w_cx, nevt, 332 );
    retrieveWf( wfs, w_cx, nevt, 344 );
    retrieveWf( wfs, w_cx, nevt, 347 );
    retrieveWf( wfs, w_cx, nevt, 349 );
    retrieveWf( wfs, w_cx, nevt, 504 );
    retrieveWf( wfs, w_cx, nevt, 508 );
    retrieveWf( wfs, w_cx, nevt, 527 );
    retrieveWf( wfs, w_cx, nevt, 540 );
    retrieveWf( wfs, w_cx, nevt, 544 );
    retrieveWf( wfs, w_cx, nevt, 551 );
    retrieveWf( wfs, w_cx, nevt, 662 );
#endif
#endif

    // *** DIAGRAM 11801 OF 15495 ***
    // Wavefunction(s) for diagram number 11801
    // (none)
    // Amplitude(s) for diagram number 11801
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[527], w_fp[211], w_fp[325], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[538] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11802 OF 15495 ***
    // Wavefunction(s) for diagram number 11802
    // (none)
    // Amplitude(s) for diagram number 11802
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[218], w_fp[508], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[522] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11803 OF 15495 ***
    // Wavefunction(s) for diagram number 11803
    // (none)
    // Amplitude(s) for diagram number 11803
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[347], w_fp[211], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[529] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11804 OF 15495 ***
    // Wavefunction(s) for diagram number 11804
    // (none)
    // Amplitude(s) for diagram number 11804
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[332], w_fp[219], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[498] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11805 OF 15495 ***
    // Wavefunction(s) for diagram number 11805
    // (none)
    // Amplitude(s) for diagram number 11805
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[168], w_fp[197], w_fp[551], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[498] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[168], w_fp[197], w_fp[504], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[487] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[168], w_fp[197], w_fp[544], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[498] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11806 OF 15495 ***
    // Wavefunction(s) for diagram number 11806
    // (none)
    // Amplitude(s) for diagram number 11806
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[214], w_fp[66], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] += amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[499] += amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[565] -= amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];

    // *** DIAGRAM 11807 OF 15495 ***
    // Wavefunction(s) for diagram number 11807
    // (none)
    // Amplitude(s) for diagram number 11807
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[222], w_fp[662], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11808 OF 15495 ***
    // Wavefunction(s) for diagram number 11808
    // (none)
    // Amplitude(s) for diagram number 11808
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[221], w_fp[197], w_fp[662], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[485] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[498] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[499] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11809 OF 15495 ***
    // Wavefunction(s) for diagram number 11809
    // (none)
    // Amplitude(s) for diagram number 11809
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[344], w_fp[540], w_fp[66], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[488] += amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[494] -= amp_sv[0];
    jamp_sv[495] += amp_sv[0];

    // *** DIAGRAM 11810 OF 15495 ***
    // Wavefunction(s) for diagram number 11810
    // (none)
    // Amplitude(s) for diagram number 11810
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[540], w_fp[349], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[485] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[498] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[499] -= cxtype( 0, 1 ) * amp_sv[0];

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

#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup1182( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 197 );
    retrieveWf( wfs, w_cx, nevt, 214 );
    retrieveWf( wfs, w_cx, nevt, 221 );
    retrieveWf( wfs, w_cx, nevt, 222 );
    retrieveWf( wfs, w_cx, nevt, 239 );
    retrieveWf( wfs, w_cx, nevt, 267 );
    retrieveWf( wfs, w_cx, nevt, 268 );
    retrieveWf( wfs, w_cx, nevt, 325 );
    retrieveWf( wfs, w_cx, nevt, 344 );
    retrieveWf( wfs, w_cx, nevt, 349 );
    retrieveWf( wfs, w_cx, nevt, 508 );
    retrieveWf( wfs, w_cx, nevt, 514 );
    retrieveWf( wfs, w_cx, nevt, 540 );
    retrieveWf( wfs, w_cx, nevt, 662 );
    retrieveWf( wfs, w_cx, nevt, 706 );
    retrieveWf( wfs, w_cx, nevt, 711 );
#endif
#endif

    // *** DIAGRAM 11811 OF 15495 ***
    // Wavefunction(s) for diagram number 11811
    // (none)
    // Amplitude(s) for diagram number 11811
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[221], w_fp[540], w_fp[325], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] += amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[499] += amp_sv[0];

    // *** DIAGRAM 11812 OF 15495 ***
    // Wavefunction(s) for diagram number 11812
    // (none)
    // Amplitude(s) for diagram number 11812
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[508], w_fp[514], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[522] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11813 OF 15495 ***
    // Wavefunction(s) for diagram number 11813
    // (none)
    // Amplitude(s) for diagram number 11813
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[344], w_fp[197], w_fp[514], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[488] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11814 OF 15495 ***
    // Wavefunction(s) for diagram number 11814
    // (none)
    // Amplitude(s) for diagram number 11814
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[514], w_fp[325], w_fp[214], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[489] += amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[587] -= amp_sv[0];

    // *** DIAGRAM 11815 OF 15495 ***
    // Wavefunction(s) for diagram number 11815
    // (none)
    // Amplitude(s) for diagram number 11815
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[221], w_fp[508], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[522] += amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];

    // *** DIAGRAM 11816 OF 15495 ***
    // Wavefunction(s) for diagram number 11816
    // (none)
    // Amplitude(s) for diagram number 11816
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[344], w_fp[222], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[540] += amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    jamp_sv[565] += amp_sv[0];

    // *** DIAGRAM 11817 OF 15495 ***
    // Wavefunction(s) for diagram number 11817
    // (none)
    // Amplitude(s) for diagram number 11817
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[349], w_fp[214], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[494] -= amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[498] += amp_sv[0];
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];

    // *** DIAGRAM 11818 OF 15495 ***
    // Wavefunction(s) for diagram number 11818
    // (none)
    // Amplitude(s) for diagram number 11818
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[197], w_fp[711], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[494] -= amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[498] += amp_sv[0];
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[197], w_fp[268], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[488] += amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[494] -= amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[565] -= amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[197], w_fp[267], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] += amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[499] += amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[565] -= amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];

    // *** DIAGRAM 11819 OF 15495 ***
    // Wavefunction(s) for diagram number 11819
    // (none)
    // Amplitude(s) for diagram number 11819
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[239], w_fp[5], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[101] += amp_sv[0];
    jamp_sv[219] -= amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[436] += amp_sv[0];
    jamp_sv[442] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[239], w_fp[5], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[219] -= amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[379] += amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[499] += amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[239], w_fp[5], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[23] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[379] += amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[499] += amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];

    // *** DIAGRAM 11820 OF 15495 ***
    // Wavefunction(s) for diagram number 11820
    // (none)
    // Amplitude(s) for diagram number 11820
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[239], w_fp[6], w_fp[706], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[219] -= amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[379] += amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[499] += amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[627] += amp_sv[0];

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


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup1183( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 198 );
    retrieveWf( wfs, w_cx, nevt, 216 );
    retrieveWf( wfs, w_cx, nevt, 239 );
    retrieveWf( wfs, w_cx, nevt, 333 );
    retrieveWf( wfs, w_cx, nevt, 334 );
    retrieveWf( wfs, w_cx, nevt, 499 );
    retrieveWf( wfs, w_cx, nevt, 500 );
    retrieveWf( wfs, w_cx, nevt, 501 );
    retrieveWf( wfs, w_cx, nevt, 560 );
    retrieveWf( wfs, w_cx, nevt, 586 );
    retrieveWf( wfs, w_cx, nevt, 662 );
    retrieveWf( wfs, w_cx, nevt, 706 );
#endif
#endif

    // *** DIAGRAM 11821 OF 15495 ***
    // Wavefunction(s) for diagram number 11821
    // (none)
    // Amplitude(s) for diagram number 11821
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[239], w_fp[5], w_fp[500], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[23] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[379] += amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[499] += amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];

    // *** DIAGRAM 11822 OF 15495 ***
    // Wavefunction(s) for diagram number 11822
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[662], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[508] );
    // Amplitude(s) for diagram number 11822
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[216], w_fp[508], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[23] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];

    // *** DIAGRAM 11823 OF 15495 ***
    // Wavefunction(s) for diagram number 11823
    // (none)
    // Amplitude(s) for diagram number 11823
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[216], w_fp[2], w_fp[500], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[221] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[485] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[499] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11824 OF 15495 ***
    // Wavefunction(s) for diagram number 11824
    // (none)
    // Amplitude(s) for diagram number 11824
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[198], w_fp[508], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[219] -= amp_sv[0];
    jamp_sv[627] += amp_sv[0];

    // *** DIAGRAM 11825 OF 15495 ***
    // Wavefunction(s) for diagram number 11825
    // (none)
    // Amplitude(s) for diagram number 11825
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[198], w_fp[2], w_fp[706], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11826 OF 15495 ***
    // Wavefunction(s) for diagram number 11826
    // (none)
    // Amplitude(s) for diagram number 11826
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[586], w_fp[501], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11827 OF 15495 ***
    // Wavefunction(s) for diagram number 11827
    // (none)
    // Amplitude(s) for diagram number 11827
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[560], w_fp[501], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[232] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11828 OF 15495 ***
    // Wavefunction(s) for diagram number 11828
    // (none)
    // Amplitude(s) for diagram number 11828
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[333], w_fp[6], w_fp[499], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[232] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11829 OF 15495 ***
    // Wavefunction(s) for diagram number 11829
    // (none)
    // Amplitude(s) for diagram number 11829
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[560], w_fp[2], w_fp[333], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[232] += amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    jamp_sv[466] += amp_sv[0];
    jamp_sv[640] -= amp_sv[0];

    // *** DIAGRAM 11830 OF 15495 ***
    // Wavefunction(s) for diagram number 11830
    // (none)
    // Amplitude(s) for diagram number 11830
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[334], w_fp[5], w_fp[499], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] += cxtype( 0, 1 ) * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 508 );
#endif
#endif
  }

  //--------------------------------------------------------------------------


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup1184( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 198 );
    retrieveWf( wfs, w_cx, nevt, 216 );
    retrieveWf( wfs, w_cx, nevt, 239 );
    retrieveWf( wfs, w_cx, nevt, 270 );
    retrieveWf( wfs, w_cx, nevt, 272 );
    retrieveWf( wfs, w_cx, nevt, 333 );
    retrieveWf( wfs, w_cx, nevt, 334 );
    retrieveWf( wfs, w_cx, nevt, 341 );
    retrieveWf( wfs, w_cx, nevt, 342 );
    retrieveWf( wfs, w_cx, nevt, 343 );
    retrieveWf( wfs, w_cx, nevt, 501 );
    retrieveWf( wfs, w_cx, nevt, 536 );
    retrieveWf( wfs, w_cx, nevt, 541 );
    retrieveWf( wfs, w_cx, nevt, 586 );
    retrieveWf( wfs, w_cx, nevt, 588 );
#endif
#endif

    // *** DIAGRAM 11831 OF 15495 ***
    // Wavefunction(s) for diagram number 11831
    // (none)
    // Amplitude(s) for diagram number 11831
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[586], w_fp[2], w_fp[334], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[238] += amp_sv[0];
    jamp_sv[526] -= amp_sv[0];
    jamp_sv[586] += amp_sv[0];
    jamp_sv[646] -= amp_sv[0];

    // *** DIAGRAM 11832 OF 15495 ***
    // Wavefunction(s) for diagram number 11832
    // (none)
    // Amplitude(s) for diagram number 11832
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[536], w_fp[2], w_fp[341], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[232] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[536], w_fp[2], w_fp[342], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[536], w_fp[2], w_fp[343], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[232] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11833 OF 15495 ***
    // Wavefunction(s) for diagram number 11833
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[501], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[507] );
    // Amplitude(s) for diagram number 11833
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[216], w_fp[507], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[221] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11834 OF 15495 ***
    // Wavefunction(s) for diagram number 11834
    // (none)
    // Amplitude(s) for diagram number 11834
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[588], w_fp[501], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11835 OF 15495 ***
    // Wavefunction(s) for diagram number 11835
    // (none)
    // Amplitude(s) for diagram number 11835
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[198], w_fp[507], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11836 OF 15495 ***
    // Wavefunction(s) for diagram number 11836
    // (none)
    // Amplitude(s) for diagram number 11836
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[541], w_fp[501], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11837 OF 15495 ***
    // Wavefunction(s) for diagram number 11837
    // (none)
    // Amplitude(s) for diagram number 11837
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[333], w_fp[239], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[229] -= amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[493] += amp_sv[0];
    jamp_sv[496] -= amp_sv[0];
    jamp_sv[499] += amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[586] -= amp_sv[0];
    jamp_sv[637] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[333], w_fp[239], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[586] -= amp_sv[0];
    jamp_sv[640] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[333], w_fp[239], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[229] += amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    jamp_sv[493] -= amp_sv[0];
    jamp_sv[496] += amp_sv[0];
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[637] -= amp_sv[0];
    jamp_sv[640] += amp_sv[0];

    // *** DIAGRAM 11838 OF 15495 ***
    // Wavefunction(s) for diagram number 11838
    // (none)
    // Amplitude(s) for diagram number 11838
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[239], w_fp[6], w_fp[270], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[229] -= amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[493] += amp_sv[0];
    jamp_sv[496] -= amp_sv[0];
    jamp_sv[499] += amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[586] -= amp_sv[0];
    jamp_sv[637] += amp_sv[0];

    // *** DIAGRAM 11839 OF 15495 ***
    // Wavefunction(s) for diagram number 11839
    // (none)
    // Amplitude(s) for diagram number 11839
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[333], w_fp[6], w_fp[272], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[586] -= amp_sv[0];
    jamp_sv[640] += amp_sv[0];

    // *** DIAGRAM 11840 OF 15495 ***
    // Wavefunction(s) for diagram number 11840
    // (none)
    // Amplitude(s) for diagram number 11840
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[198], w_fp[2], w_fp[270], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] -= cxtype( 0, 1 ) * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 507 );
#endif
#endif
  }

  //--------------------------------------------------------------------------


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup1185( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 198 );
    retrieveWf( wfs, w_cx, nevt, 216 );
    retrieveWf( wfs, w_cx, nevt, 239 );
    retrieveWf( wfs, w_cx, nevt, 272 );
    retrieveWf( wfs, w_cx, nevt, 274 );
    retrieveWf( wfs, w_cx, nevt, 333 );
    retrieveWf( wfs, w_cx, nevt, 334 );
    retrieveWf( wfs, w_cx, nevt, 503 );
    retrieveWf( wfs, w_cx, nevt, 541 );
    retrieveWf( wfs, w_cx, nevt, 588 );
    retrieveWf( wfs, w_cx, nevt, 670 );
    retrieveWf( wfs, w_cx, nevt, 703 );
    retrieveWf( wfs, w_cx, nevt, 709 );
    retrieveWf( wfs, w_cx, nevt, 710 );
    retrieveWf( wfs, w_cx, nevt, 712 );
#endif
#endif

    // *** DIAGRAM 11841 OF 15495 ***
    // Wavefunction(s) for diagram number 11841
    // (none)
    // Amplitude(s) for diagram number 11841
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[541], w_fp[2], w_fp[333], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[229] += amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[637] -= amp_sv[0];

    // *** DIAGRAM 11842 OF 15495 ***
    // Wavefunction(s) for diagram number 11842
    // (none)
    // Amplitude(s) for diagram number 11842
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[334], w_fp[239], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[23] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[373] += amp_sv[0];
    jamp_sv[376] -= amp_sv[0];
    jamp_sv[379] += amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[334], w_fp[239], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[23] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[586] -= amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[334], w_fp[239], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[235] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[373] -= amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[586] -= amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[646] += amp_sv[0];

    // *** DIAGRAM 11843 OF 15495 ***
    // Wavefunction(s) for diagram number 11843
    // (none)
    // Amplitude(s) for diagram number 11843
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[239], w_fp[5], w_fp[503], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[23] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[373] += amp_sv[0];
    jamp_sv[376] -= amp_sv[0];
    jamp_sv[379] += amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[643] += amp_sv[0];

    // *** DIAGRAM 11844 OF 15495 ***
    // Wavefunction(s) for diagram number 11844
    // (none)
    // Amplitude(s) for diagram number 11844
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[334], w_fp[5], w_fp[272], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[23] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[586] -= amp_sv[0];
    jamp_sv[646] += amp_sv[0];

    // *** DIAGRAM 11845 OF 15495 ***
    // Wavefunction(s) for diagram number 11845
    // (none)
    // Amplitude(s) for diagram number 11845
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[216], w_fp[2], w_fp[503], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[77] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11846 OF 15495 ***
    // Wavefunction(s) for diagram number 11846
    // (none)
    // Amplitude(s) for diagram number 11846
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[588], w_fp[2], w_fp[334], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[235] += amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[643] -= amp_sv[0];

    // *** DIAGRAM 11847 OF 15495 ***
    // Wavefunction(s) for diagram number 11847
    // (none)
    // Amplitude(s) for diagram number 11847
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[709], w_fp[239], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[229] -= amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[493] += amp_sv[0];
    jamp_sv[496] -= amp_sv[0];
    jamp_sv[499] += amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[586] -= amp_sv[0];
    jamp_sv[637] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[703], w_fp[239], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[219] += amp_sv[0];
    jamp_sv[229] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[493] += amp_sv[0];
    jamp_sv[496] -= amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[553] -= amp_sv[0];
    jamp_sv[556] += amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[586] -= amp_sv[0];
    jamp_sv[627] -= amp_sv[0];
    jamp_sv[637] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[670], w_fp[239], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[99] += amp_sv[0];
    jamp_sv[219] += amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[553] -= amp_sv[0];
    jamp_sv[556] += amp_sv[0];
    jamp_sv[562] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[627] -= amp_sv[0];

    // *** DIAGRAM 11848 OF 15495 ***
    // Wavefunction(s) for diagram number 11848
    // (none)
    // Amplitude(s) for diagram number 11848
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[198], w_fp[2], w_fp[709], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[198], w_fp[2], w_fp[703], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[198], w_fp[2], w_fp[670], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11849 OF 15495 ***
    // Wavefunction(s) for diagram number 11849
    // (none)
    // Amplitude(s) for diagram number 11849
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[274], w_fp[239], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[23] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[373] += amp_sv[0];
    jamp_sv[376] -= amp_sv[0];
    jamp_sv[379] += amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[710], w_fp[239], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[373] += amp_sv[0];
    jamp_sv[376] -= amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[436] += amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[712], w_fp[239], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[101] += amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[436] += amp_sv[0];
    jamp_sv[442] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];

    // *** DIAGRAM 11850 OF 15495 ***
    // Wavefunction(s) for diagram number 11850
    // (none)
    // Amplitude(s) for diagram number 11850
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[216], w_fp[2], w_fp[274], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[77] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[216], w_fp[2], w_fp[710], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[77] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[221] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[485] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[499] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[216], w_fp[2], w_fp[712], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[221] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[485] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[499] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] += cxtype( 0, 1 ) * amp_sv[0];

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


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup1186( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 113 );
    retrieveWf( wfs, w_cx, nevt, 196 );
    retrieveWf( wfs, w_cx, nevt, 239 );
    retrieveWf( wfs, w_cx, nevt, 243 );
    retrieveWf( wfs, w_cx, nevt, 244 );
    retrieveWf( wfs, w_cx, nevt, 325 );
    retrieveWf( wfs, w_cx, nevt, 341 );
    retrieveWf( wfs, w_cx, nevt, 342 );
    retrieveWf( wfs, w_cx, nevt, 343 );
    retrieveWf( wfs, w_cx, nevt, 346 );
    retrieveWf( wfs, w_cx, nevt, 351 );
    retrieveWf( wfs, w_cx, nevt, 501 );
    retrieveWf( wfs, w_cx, nevt, 536 );
    retrieveWf( wfs, w_cx, nevt, 559 );
    retrieveWf( wfs, w_cx, nevt, 662 );
#endif
#endif

    // *** DIAGRAM 11851 OF 15495 ***
    // Wavefunction(s) for diagram number 11851
    // (none)
    // Amplitude(s) for diagram number 11851
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[341], w_fp[239], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[85] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[99] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[562] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[640] -= amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[342], w_fp[239], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[23] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[586] -= amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[343], w_fp[239], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[586] -= amp_sv[0];
    jamp_sv[640] += amp_sv[0];

    // *** DIAGRAM 11852 OF 15495 ***
    // Wavefunction(s) for diagram number 11852
    // (none)
    // Amplitude(s) for diagram number 11852
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[239], w_fp[113], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[101] += amp_sv[0];
    jamp_sv[219] -= amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[436] += amp_sv[0];
    jamp_sv[442] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];

    // *** DIAGRAM 11853 OF 15495 ***
    // Wavefunction(s) for diagram number 11853
    // (none)
    // Amplitude(s) for diagram number 11853
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[244], w_fp[662], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[433] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[436] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11854 OF 15495 ***
    // Wavefunction(s) for diagram number 11854
    // (none)
    // Amplitude(s) for diagram number 11854
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[243], w_fp[2], w_fp[662], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[221] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11855 OF 15495 ***
    // Wavefunction(s) for diagram number 11855
    // (none)
    // Amplitude(s) for diagram number 11855
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[536], w_fp[501], w_fp[113], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[232] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[640] -= amp_sv[0];
    jamp_sv[646] += amp_sv[0];

    // *** DIAGRAM 11856 OF 15495 ***
    // Wavefunction(s) for diagram number 11856
    // (none)
    // Amplitude(s) for diagram number 11856
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[536], w_fp[2], w_fp[351], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[232] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11857 OF 15495 ***
    // Wavefunction(s) for diagram number 11857
    // (none)
    // Amplitude(s) for diagram number 11857
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[536], w_fp[244], w_fp[325], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[442] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];

    // *** DIAGRAM 11858 OF 15495 ***
    // Wavefunction(s) for diagram number 11858
    // (none)
    // Amplitude(s) for diagram number 11858
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[501], w_fp[559], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[221] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11859 OF 15495 ***
    // Wavefunction(s) for diagram number 11859
    // (none)
    // Amplitude(s) for diagram number 11859
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[346], w_fp[2], w_fp[559], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[433] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[436] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11860 OF 15495 ***
    // Wavefunction(s) for diagram number 11860
    // (none)
    // Amplitude(s) for diagram number 11860
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[559], w_fp[325], w_fp[239], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[61] -= amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[219] += amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[553] -= amp_sv[0];
    jamp_sv[556] += amp_sv[0];
    jamp_sv[627] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[640] += amp_sv[0];
    jamp_sv[646] -= amp_sv[0];

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


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup1187( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 154 );
    retrieveWf( wfs, w_cx, nevt, 170 );
    retrieveWf( wfs, w_cx, nevt, 196 );
    retrieveWf( wfs, w_cx, nevt, 218 );
    retrieveWf( wfs, w_cx, nevt, 239 );
    retrieveWf( wfs, w_cx, nevt, 243 );
    retrieveWf( wfs, w_cx, nevt, 244 );
    retrieveWf( wfs, w_cx, nevt, 346 );
    retrieveWf( wfs, w_cx, nevt, 351 );
    retrieveWf( wfs, w_cx, nevt, 498 );
    retrieveWf( wfs, w_cx, nevt, 500 );
    retrieveWf( wfs, w_cx, nevt, 501 );
    retrieveWf( wfs, w_cx, nevt, 508 );
    retrieveWf( wfs, w_cx, nevt, 548 );
    retrieveWf( wfs, w_cx, nevt, 662 );
    retrieveWf( wfs, w_cx, nevt, 708 );
    retrieveWf( wfs, w_cx, nevt, 748 );
#endif
#endif

    // *** DIAGRAM 11861 OF 15495 ***
    // Wavefunction(s) for diagram number 11861
    // (none)
    // Amplitude(s) for diagram number 11861
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[243], w_fp[501], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[219] += amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[627] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];

    // *** DIAGRAM 11862 OF 15495 ***
    // Wavefunction(s) for diagram number 11862
    // (none)
    // Amplitude(s) for diagram number 11862
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[346], w_fp[244], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[433] += amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[553] -= amp_sv[0];
    jamp_sv[556] += amp_sv[0];

    // *** DIAGRAM 11863 OF 15495 ***
    // Wavefunction(s) for diagram number 11863
    // (none)
    // Amplitude(s) for diagram number 11863
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[351], w_fp[239], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[85] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[99] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[562] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[640] -= amp_sv[0];
    jamp_sv[646] += amp_sv[0];

    // *** DIAGRAM 11864 OF 15495 ***
    // Wavefunction(s) for diagram number 11864
    // (none)
    // Amplitude(s) for diagram number 11864
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[2], w_fp[498], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[85] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[99] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[562] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[640] -= amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[2], w_fp[748], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[61] += amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[85] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[219] -= amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[436] += amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    jamp_sv[640] -= amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[2], w_fp[708], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[101] += amp_sv[0];
    jamp_sv[219] -= amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[436] += amp_sv[0];
    jamp_sv[442] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];

    // *** DIAGRAM 11865 OF 15495 ***
    // Wavefunction(s) for diagram number 11865
    // (none)
    // Amplitude(s) for diagram number 11865
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[154], w_fp[4], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[316] += amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[154], w_fp[4], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[259] += amp_sv[0];
    jamp_sv[283] += amp_sv[0];
    jamp_sv[343] -= amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[498] += amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[154], w_fp[4], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] += amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[259] += amp_sv[0];
    jamp_sv[283] += amp_sv[0];
    jamp_sv[313] += amp_sv[0];
    jamp_sv[316] -= amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[343] -= amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[498] += amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[628] += amp_sv[0];

    // *** DIAGRAM 11866 OF 15495 ***
    // Wavefunction(s) for diagram number 11866
    // (none)
    // Amplitude(s) for diagram number 11866
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[154], w_fp[6], w_fp[548], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[259] += amp_sv[0];
    jamp_sv[283] += amp_sv[0];
    jamp_sv[343] -= amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[498] += amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];

    // *** DIAGRAM 11867 OF 15495 ***
    // Wavefunction(s) for diagram number 11867
    // (none)
    // Amplitude(s) for diagram number 11867
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[154], w_fp[4], w_fp[500], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] += amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[259] += amp_sv[0];
    jamp_sv[283] += amp_sv[0];
    jamp_sv[313] += amp_sv[0];
    jamp_sv[316] -= amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[343] -= amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[498] += amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[628] += amp_sv[0];

    // *** DIAGRAM 11868 OF 15495 ***
    // Wavefunction(s) for diagram number 11868
    // (none)
    // Amplitude(s) for diagram number 11868
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[218], w_fp[508], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] += amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[628] += amp_sv[0];

    // *** DIAGRAM 11869 OF 15495 ***
    // Wavefunction(s) for diagram number 11869
    // (none)
    // Amplitude(s) for diagram number 11869
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[218], w_fp[2], w_fp[500], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[484] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[498] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11870 OF 15495 ***
    // Wavefunction(s) for diagram number 11870
    // (none)
    // Amplitude(s) for diagram number 11870
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[170], w_fp[508], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];

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


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup1188( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 170 );
    retrieveWf( wfs, w_cx, nevt, 218 );
    retrieveWf( wfs, w_cx, nevt, 332 );
    retrieveWf( wfs, w_cx, nevt, 334 );
    retrieveWf( wfs, w_cx, nevt, 338 );
    retrieveWf( wfs, w_cx, nevt, 339 );
    retrieveWf( wfs, w_cx, nevt, 340 );
    retrieveWf( wfs, w_cx, nevt, 481 );
    retrieveWf( wfs, w_cx, nevt, 501 );
    retrieveWf( wfs, w_cx, nevt, 507 );
    retrieveWf( wfs, w_cx, nevt, 509 );
    retrieveWf( wfs, w_cx, nevt, 527 );
    retrieveWf( wfs, w_cx, nevt, 538 );
    retrieveWf( wfs, w_cx, nevt, 548 );
    retrieveWf( wfs, w_cx, nevt, 747 );
#endif
#endif

    // *** DIAGRAM 11871 OF 15495 ***
    // Wavefunction(s) for diagram number 11871
    // (none)
    // Amplitude(s) for diagram number 11871
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[170], w_fp[2], w_fp[548], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[217] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[283] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[343] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11872 OF 15495 ***
    // Wavefunction(s) for diagram number 11872
    // (none)
    // Amplitude(s) for diagram number 11872
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[509], w_fp[501], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11873 OF 15495 ***
    // Wavefunction(s) for diagram number 11873
    // (none)
    // Amplitude(s) for diagram number 11873
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[538], w_fp[501], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11874 OF 15495 ***
    // Wavefunction(s) for diagram number 11874
    // (none)
    // Amplitude(s) for diagram number 11874
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[332], w_fp[6], w_fp[747], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11875 OF 15495 ***
    // Wavefunction(s) for diagram number 11875
    // (none)
    // Amplitude(s) for diagram number 11875
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[538], w_fp[2], w_fp[332], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[226] += amp_sv[0];
    jamp_sv[286] -= amp_sv[0];
    jamp_sv[346] += amp_sv[0];
    jamp_sv[634] -= amp_sv[0];

    // *** DIAGRAM 11876 OF 15495 ***
    // Wavefunction(s) for diagram number 11876
    // (none)
    // Amplitude(s) for diagram number 11876
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[334], w_fp[4], w_fp[747], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11877 OF 15495 ***
    // Wavefunction(s) for diagram number 11877
    // (none)
    // Amplitude(s) for diagram number 11877
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[509], w_fp[2], w_fp[334], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[236] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];

    // *** DIAGRAM 11878 OF 15495 ***
    // Wavefunction(s) for diagram number 11878
    // (none)
    // Amplitude(s) for diagram number 11878
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[527], w_fp[2], w_fp[338], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[527], w_fp[2], w_fp[339], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[527], w_fp[2], w_fp[340], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11879 OF 15495 ***
    // Wavefunction(s) for diagram number 11879
    // (none)
    // Amplitude(s) for diagram number 11879
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[218], w_fp[507], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[220] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11880 OF 15495 ***
    // Wavefunction(s) for diagram number 11880
    // (none)
    // Amplitude(s) for diagram number 11880
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[481], w_fp[501], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];

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


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup1189( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 154 );
    retrieveWf( wfs, w_cx, nevt, 170 );
    retrieveWf( wfs, w_cx, nevt, 332 );
    retrieveWf( wfs, w_cx, nevt, 334 );
    retrieveWf( wfs, w_cx, nevt, 501 );
    retrieveWf( wfs, w_cx, nevt, 503 );
    retrieveWf( wfs, w_cx, nevt, 507 );
    retrieveWf( wfs, w_cx, nevt, 535 );
    retrieveWf( wfs, w_cx, nevt, 743 );
    retrieveWf( wfs, w_cx, nevt, 749 );
#endif
#endif

    // *** DIAGRAM 11881 OF 15495 ***
    // Wavefunction(s) for diagram number 11881
    // (none)
    // Amplitude(s) for diagram number 11881
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[170], w_fp[507], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[217] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11882 OF 15495 ***
    // Wavefunction(s) for diagram number 11882
    // (none)
    // Amplitude(s) for diagram number 11882
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[535], w_fp[501], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[223] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11883 OF 15495 ***
    // Wavefunction(s) for diagram number 11883
    // (none)
    // Amplitude(s) for diagram number 11883
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[332], w_fp[154], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    jamp_sv[283] += amp_sv[0];
    jamp_sv[343] -= amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[487] += amp_sv[0];
    jamp_sv[490] -= amp_sv[0];
    jamp_sv[498] += amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[332], w_fp[154], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[90] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    jamp_sv[286] += amp_sv[0];
    jamp_sv[346] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[332], w_fp[154], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[90] += amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    jamp_sv[283] -= amp_sv[0];
    jamp_sv[286] += amp_sv[0];
    jamp_sv[343] += amp_sv[0];
    jamp_sv[346] -= amp_sv[0];
    jamp_sv[484] += amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[634] += amp_sv[0];

    // *** DIAGRAM 11884 OF 15495 ***
    // Wavefunction(s) for diagram number 11884
    // (none)
    // Amplitude(s) for diagram number 11884
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[154], w_fp[6], w_fp[743], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    jamp_sv[283] += amp_sv[0];
    jamp_sv[343] -= amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[487] += amp_sv[0];
    jamp_sv[490] -= amp_sv[0];
    jamp_sv[498] += amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[631] += amp_sv[0];

    // *** DIAGRAM 11885 OF 15495 ***
    // Wavefunction(s) for diagram number 11885
    // (none)
    // Amplitude(s) for diagram number 11885
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[332], w_fp[6], w_fp[749], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[90] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    jamp_sv[286] += amp_sv[0];
    jamp_sv[346] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[634] += amp_sv[0];

    // *** DIAGRAM 11886 OF 15495 ***
    // Wavefunction(s) for diagram number 11886
    // (none)
    // Amplitude(s) for diagram number 11886
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[170], w_fp[2], w_fp[743], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[43] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[283] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[343] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11887 OF 15495 ***
    // Wavefunction(s) for diagram number 11887
    // (none)
    // Amplitude(s) for diagram number 11887
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[535], w_fp[2], w_fp[332], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[223] += amp_sv[0];
    jamp_sv[283] -= amp_sv[0];
    jamp_sv[343] += amp_sv[0];
    jamp_sv[631] -= amp_sv[0];

    // *** DIAGRAM 11888 OF 15495 ***
    // Wavefunction(s) for diagram number 11888
    // (none)
    // Amplitude(s) for diagram number 11888
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[334], w_fp[154], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] += amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[90] += amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[234] -= amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[253] += amp_sv[0];
    jamp_sv[256] -= amp_sv[0];
    jamp_sv[259] += amp_sv[0];
    jamp_sv[286] += amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[346] -= amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[334], w_fp[154], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[90] += amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[286] += amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[346] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[334], w_fp[154], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[253] -= amp_sv[0];
    jamp_sv[256] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];

    // *** DIAGRAM 11889 OF 15495 ***
    // Wavefunction(s) for diagram number 11889
    // (none)
    // Amplitude(s) for diagram number 11889
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[154], w_fp[4], w_fp[503], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] += amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[90] += amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[234] -= amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[253] += amp_sv[0];
    jamp_sv[256] -= amp_sv[0];
    jamp_sv[259] += amp_sv[0];
    jamp_sv[286] += amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[346] -= amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[642] += amp_sv[0];

    // *** DIAGRAM 11890 OF 15495 ***
    // Wavefunction(s) for diagram number 11890
    // (none)
    // Amplitude(s) for diagram number 11890
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[334], w_fp[4], w_fp[749], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[90] += amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[286] += amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[346] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];

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


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup1190( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 86 );
    retrieveWf( wfs, w_cx, nevt, 118 );
    retrieveWf( wfs, w_cx, nevt, 154 );
    retrieveWf( wfs, w_cx, nevt, 168 );
    retrieveWf( wfs, w_cx, nevt, 170 );
    retrieveWf( wfs, w_cx, nevt, 200 );
    retrieveWf( wfs, w_cx, nevt, 218 );
    retrieveWf( wfs, w_cx, nevt, 274 );
    retrieveWf( wfs, w_cx, nevt, 334 );
    retrieveWf( wfs, w_cx, nevt, 338 );
    retrieveWf( wfs, w_cx, nevt, 339 );
    retrieveWf( wfs, w_cx, nevt, 340 );
    retrieveWf( wfs, w_cx, nevt, 481 );
    retrieveWf( wfs, w_cx, nevt, 503 );
    retrieveWf( wfs, w_cx, nevt, 504 );
    retrieveWf( wfs, w_cx, nevt, 544 );
    retrieveWf( wfs, w_cx, nevt, 551 );
    retrieveWf( wfs, w_cx, nevt, 662 );
    retrieveWf( wfs, w_cx, nevt, 710 );
    retrieveWf( wfs, w_cx, nevt, 712 );
#endif
#endif

    // *** DIAGRAM 11891 OF 15495 ***
    // Wavefunction(s) for diagram number 11891
    // (none)
    // Amplitude(s) for diagram number 11891
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[218], w_fp[2], w_fp[503], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[76] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[90] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11892 OF 15495 ***
    // Wavefunction(s) for diagram number 11892
    // (none)
    // Amplitude(s) for diagram number 11892
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[481], w_fp[2], w_fp[334], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[234] += amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[642] -= amp_sv[0];

    // *** DIAGRAM 11893 OF 15495 ***
    // Wavefunction(s) for diagram number 11893
    // (none)
    // Amplitude(s) for diagram number 11893
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[551], w_fp[154], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    jamp_sv[283] += amp_sv[0];
    jamp_sv[343] -= amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[487] += amp_sv[0];
    jamp_sv[490] -= amp_sv[0];
    jamp_sv[498] += amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[504], w_fp[154], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[217] += amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[487] += amp_sv[0];
    jamp_sv[490] -= amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[544], w_fp[154], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[217] += amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[283] -= amp_sv[0];
    jamp_sv[343] += amp_sv[0];
    jamp_sv[484] += amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[538] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[625] -= amp_sv[0];

    // *** DIAGRAM 11894 OF 15495 ***
    // Wavefunction(s) for diagram number 11894
    // (none)
    // Amplitude(s) for diagram number 11894
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[170], w_fp[2], w_fp[551], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[43] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[283] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[343] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[170], w_fp[2], w_fp[504], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[29] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[43] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[217] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[170], w_fp[2], w_fp[544], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[217] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[283] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[343] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11895 OF 15495 ***
    // Wavefunction(s) for diagram number 11895
    // (none)
    // Amplitude(s) for diagram number 11895
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[274], w_fp[154], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] += amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[90] += amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[234] -= amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[253] += amp_sv[0];
    jamp_sv[256] -= amp_sv[0];
    jamp_sv[259] += amp_sv[0];
    jamp_sv[286] += amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[346] -= amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[710], w_fp[154], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[90] += amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[234] -= amp_sv[0];
    jamp_sv[253] += amp_sv[0];
    jamp_sv[256] -= amp_sv[0];
    jamp_sv[283] -= amp_sv[0];
    jamp_sv[286] += amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[316] += amp_sv[0];
    jamp_sv[343] += amp_sv[0];
    jamp_sv[346] -= amp_sv[0];
    jamp_sv[484] += amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[712], w_fp[154], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[283] -= amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[316] += amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[343] += amp_sv[0];
    jamp_sv[484] += amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];

    // *** DIAGRAM 11896 OF 15495 ***
    // Wavefunction(s) for diagram number 11896
    // (none)
    // Amplitude(s) for diagram number 11896
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[218], w_fp[2], w_fp[274], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[76] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[90] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[218], w_fp[2], w_fp[710], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[76] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[90] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[484] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[498] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[218], w_fp[2], w_fp[712], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[484] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[498] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11897 OF 15495 ***
    // Wavefunction(s) for diagram number 11897
    // (none)
    // Amplitude(s) for diagram number 11897
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[338], w_fp[154], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[538] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[634] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[339], w_fp[154], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[90] += amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[286] += amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[346] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[340], w_fp[154], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[90] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    jamp_sv[286] += amp_sv[0];
    jamp_sv[346] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[634] += amp_sv[0];

    // *** DIAGRAM 11898 OF 15495 ***
    // Wavefunction(s) for diagram number 11898
    // (none)
    // Amplitude(s) for diagram number 11898
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[154], w_fp[86], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[316] += amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];

    // *** DIAGRAM 11899 OF 15495 ***
    // Wavefunction(s) for diagram number 11899
    // (none)
    // Amplitude(s) for diagram number 11899
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[168], w_fp[118], w_fp[662], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[313] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11900 OF 15495 ***
    // Wavefunction(s) for diagram number 11900
    // (none)
    // Amplitude(s) for diagram number 11900
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[200], w_fp[2], w_fp[662], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[217] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] += cxtype( 0, 1 ) * amp_sv[0];

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
