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
  diagramgroup351( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 8 );
    retrieveWf( wfs, w_cx, nevt, 114 );
    retrieveWf( wfs, w_cx, nevt, 115 );
    retrieveWf( wfs, w_cx, nevt, 116 );
    retrieveWf( wfs, w_cx, nevt, 237 );
    retrieveWf( wfs, w_cx, nevt, 359 );
    retrieveWf( wfs, w_cx, nevt, 469 );
    retrieveWf( wfs, w_cx, nevt, 470 );
    retrieveWf( wfs, w_cx, nevt, 478 );
    retrieveWf( wfs, w_cx, nevt, 479 );
    retrieveWf( wfs, w_cx, nevt, 481 );
    retrieveWf( wfs, w_cx, nevt, 487 );
    retrieveWf( wfs, w_cx, nevt, 505 );
    retrieveWf( wfs, w_cx, nevt, 510 );
#endif
#endif

    // *** DIAGRAM 3501 OF 15495 ***
    // Wavefunction(s) for diagram number 3501
    // (none)
    // Amplitude(s) for diagram number 3501
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[510], w_fp[2], w_fp[115], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[305] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];

    // *** DIAGRAM 3502 OF 15495 ***
    // Wavefunction(s) for diagram number 3502
    // (none)
    // Amplitude(s) for diagram number 3502
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[4], w_fp[116], w_fp[470], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[305] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3503 OF 15495 ***
    // Wavefunction(s) for diagram number 3503
    // (none)
    // Amplitude(s) for diagram number 3503
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[478], w_fp[2], w_fp[116], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[455] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[695] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];

    // *** DIAGRAM 3504 OF 15495 ***
    // Wavefunction(s) for diagram number 3504
    // (none)
    // Amplitude(s) for diagram number 3504
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[2], w_fp[237], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[305] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[2], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[353] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[665] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[671] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[2], w_fp[114], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[305] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[665] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[671] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3505 OF 15495 ***
    // Wavefunction(s) for diagram number 3505
    // (none)
    // Amplitude(s) for diagram number 3505
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[359], w_fp[4], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[183] -= amp_sv[0];
    jamp_sv[185] += amp_sv[0];
    jamp_sv[207] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[347] -= amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[443] -= amp_sv[0];
    jamp_sv[561] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[657] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[671] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[359], w_fp[4], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[303] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[347] -= amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[443] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[359], w_fp[4], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[183] += amp_sv[0];
    jamp_sv[207] -= amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[303] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[441] -= amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    jamp_sv[659] -= amp_sv[0];
    jamp_sv[665] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];

    // *** DIAGRAM 3506 OF 15495 ***
    // Wavefunction(s) for diagram number 3506
    // (none)
    // Amplitude(s) for diagram number 3506
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[359], w_fp[7], w_fp[487], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[303] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[347] -= amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[443] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];

    // *** DIAGRAM 3507 OF 15495 ***
    // Wavefunction(s) for diagram number 3507
    // (none)
    // Amplitude(s) for diagram number 3507
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[359], w_fp[4], w_fp[481], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[183] += amp_sv[0];
    jamp_sv[207] -= amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[303] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[441] -= amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    jamp_sv[659] -= amp_sv[0];
    jamp_sv[665] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];

    // *** DIAGRAM 3508 OF 15495 ***
    // Wavefunction(s) for diagram number 3508
    // (none)
    // Amplitude(s) for diagram number 3508
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[115], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[159] += amp_sv[0];
    jamp_sv[183] += amp_sv[0];
    jamp_sv[207] -= amp_sv[0];
    jamp_sv[225] += amp_sv[0];
    jamp_sv[227] -= amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[665] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[115], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[225] += amp_sv[0];
    jamp_sv[227] -= amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[327] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[447] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[635] += amp_sv[0];
    jamp_sv[641] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[115], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += amp_sv[0];
    jamp_sv[159] -= amp_sv[0];
    jamp_sv[183] -= amp_sv[0];
    jamp_sv[207] += amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[327] += amp_sv[0];
    jamp_sv[447] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[635] += amp_sv[0];
    jamp_sv[641] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[671] -= amp_sv[0];
    jamp_sv[695] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];

    // *** DIAGRAM 3509 OF 15495 ***
    // Wavefunction(s) for diagram number 3509
    // (none)
    // Amplitude(s) for diagram number 3509
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[115], w_fp[7], w_fp[469], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[159] += amp_sv[0];
    jamp_sv[183] += amp_sv[0];
    jamp_sv[207] -= amp_sv[0];
    jamp_sv[225] += amp_sv[0];
    jamp_sv[227] -= amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[665] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];

    // *** DIAGRAM 3510 OF 15495 ***
    // Wavefunction(s) for diagram number 3510
    // (none)
    // Amplitude(s) for diagram number 3510
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[115], w_fp[481], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += amp_sv[0];
    jamp_sv[159] -= amp_sv[0];
    jamp_sv[183] -= amp_sv[0];
    jamp_sv[207] += amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[327] += amp_sv[0];
    jamp_sv[447] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[635] += amp_sv[0];
    jamp_sv[641] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[671] -= amp_sv[0];
    jamp_sv[695] -= amp_sv[0];
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
  diagramgroup352( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 8 );
    retrieveWf( wfs, w_cx, nevt, 72 );
    retrieveWf( wfs, w_cx, nevt, 79 );
    retrieveWf( wfs, w_cx, nevt, 114 );
    retrieveWf( wfs, w_cx, nevt, 116 );
    retrieveWf( wfs, w_cx, nevt, 237 );
    retrieveWf( wfs, w_cx, nevt, 244 );
    retrieveWf( wfs, w_cx, nevt, 262 );
    retrieveWf( wfs, w_cx, nevt, 353 );
    retrieveWf( wfs, w_cx, nevt, 359 );
    retrieveWf( wfs, w_cx, nevt, 374 );
    retrieveWf( wfs, w_cx, nevt, 375 );
    retrieveWf( wfs, w_cx, nevt, 438 );
    retrieveWf( wfs, w_cx, nevt, 455 );
    retrieveWf( wfs, w_cx, nevt, 469 );
    retrieveWf( wfs, w_cx, nevt, 479 );
    retrieveWf( wfs, w_cx, nevt, 487 );
    retrieveWf( wfs, w_cx, nevt, 495 );
#endif
#endif

    // *** DIAGRAM 3511 OF 15495 ***
    // Wavefunction(s) for diagram number 3511
    // (none)
    // Amplitude(s) for diagram number 3511
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[4], w_fp[116], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[159] += amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[185] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[4], w_fp[116], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[4], w_fp[116], COUPs[2], 1.0, &amp_fp[0] );
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

    // *** DIAGRAM 3512 OF 15495 ***
    // Wavefunction(s) for diagram number 3512
    // (none)
    // Amplitude(s) for diagram number 3512
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[4], w_fp[116], w_fp[469], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[159] += amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[185] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];

    // *** DIAGRAM 3513 OF 15495 ***
    // Wavefunction(s) for diagram number 3513
    // (none)
    // Amplitude(s) for diagram number 3513
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[116], w_fp[487], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];

    // *** DIAGRAM 3514 OF 15495 ***
    // Wavefunction(s) for diagram number 3514
    // (none)
    // Amplitude(s) for diagram number 3514
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[262], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[159] += amp_sv[0];
    jamp_sv[183] += amp_sv[0];
    jamp_sv[207] -= amp_sv[0];
    jamp_sv[303] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[633] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];
    jamp_sv[665] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[353], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[183] += amp_sv[0];
    jamp_sv[207] -= amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[303] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[441] -= amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    jamp_sv[659] -= amp_sv[0];
    jamp_sv[665] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[72], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += amp_sv[0];
    jamp_sv[159] -= amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[441] -= amp_sv[0];
    jamp_sv[447] += amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[635] += amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    jamp_sv[659] -= amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    jamp_sv[695] -= amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];

    // *** DIAGRAM 3515 OF 15495 ***
    // Wavefunction(s) for diagram number 3515
    // (none)
    // Amplitude(s) for diagram number 3515
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[79], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[209] += amp_sv[0];
    jamp_sv[233] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[353] += amp_sv[0];
    jamp_sv[359] -= amp_sv[0];
    jamp_sv[453] += amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[717] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[374], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[233] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];
    jamp_sv[303] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[347] -= amp_sv[0];
    jamp_sv[443] -= amp_sv[0];
    jamp_sv[453] += amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[375], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[303] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[347] -= amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[443] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];

    // *** DIAGRAM 3516 OF 15495 ***
    // Wavefunction(s) for diagram number 3516
    // (none)
    // Amplitude(s) for diagram number 3516
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[237], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[159] += amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[185] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[8], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] += amp_sv[0];
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[183] -= amp_sv[0];
    jamp_sv[185] += amp_sv[0];
    jamp_sv[207] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[225] -= amp_sv[0];
    jamp_sv[227] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[671] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[114], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += amp_sv[0];
    jamp_sv[159] -= amp_sv[0];
    jamp_sv[183] -= amp_sv[0];
    jamp_sv[207] += amp_sv[0];
    jamp_sv[225] -= amp_sv[0];
    jamp_sv[227] += amp_sv[0];
    jamp_sv[233] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[671] -= amp_sv[0];
    jamp_sv[695] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];

    // *** DIAGRAM 3517 OF 15495 ***
    // Wavefunction(s) for diagram number 3517
    // (none)
    // Amplitude(s) for diagram number 3517
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[359], w_fp[7], w_fp[495], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[209] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[443] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[563] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3518 OF 15495 ***
    // Wavefunction(s) for diagram number 3518
    // (none)
    // Amplitude(s) for diagram number 3518
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[438], w_fp[2], w_fp[359], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[443] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];

    // *** DIAGRAM 3519 OF 15495 ***
    // Wavefunction(s) for diagram number 3519
    // (none)
    // Amplitude(s) for diagram number 3519
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[455], w_fp[244], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[453] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3520 OF 15495 ***
    // Wavefunction(s) for diagram number 3520
    // (none)
    // Amplitude(s) for diagram number 3520
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[438], w_fp[244], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[443] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[563] -= cxtype( 0, 1 ) * amp_sv[0];

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
  diagramgroup353( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 72 );
    retrieveWf( wfs, w_cx, nevt, 79 );
    retrieveWf( wfs, w_cx, nevt, 115 );
    retrieveWf( wfs, w_cx, nevt, 116 );
    retrieveWf( wfs, w_cx, nevt, 244 );
    retrieveWf( wfs, w_cx, nevt, 262 );
    retrieveWf( wfs, w_cx, nevt, 353 );
    retrieveWf( wfs, w_cx, nevt, 359 );
    retrieveWf( wfs, w_cx, nevt, 374 );
    retrieveWf( wfs, w_cx, nevt, 375 );
    retrieveWf( wfs, w_cx, nevt, 435 );
    retrieveWf( wfs, w_cx, nevt, 451 );
    retrieveWf( wfs, w_cx, nevt, 454 );
    retrieveWf( wfs, w_cx, nevt, 455 );
    retrieveWf( wfs, w_cx, nevt, 486 );
    retrieveWf( wfs, w_cx, nevt, 495 );
    retrieveWf( wfs, w_cx, nevt, 515 );
#endif
#endif

    // *** DIAGRAM 3521 OF 15495 ***
    // Wavefunction(s) for diagram number 3521
    // (none)
    // Amplitude(s) for diagram number 3521
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[455], w_fp[2], w_fp[116], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[453] += amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[717] += amp_sv[0];

    // *** DIAGRAM 3522 OF 15495 ***
    // Wavefunction(s) for diagram number 3522
    // (none)
    // Amplitude(s) for diagram number 3522
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[116], w_fp[495], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[209] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[453] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3523 OF 15495 ***
    // Wavefunction(s) for diagram number 3523
    // (none)
    // Amplitude(s) for diagram number 3523
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[435], w_fp[2], w_fp[79], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[209] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[453] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[435], w_fp[2], w_fp[374], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[233] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[443] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[453] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[563] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[435], w_fp[2], w_fp[375], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[209] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[443] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[563] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3524 OF 15495 ***
    // Wavefunction(s) for diagram number 3524
    // (none)
    // Amplitude(s) for diagram number 3524
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[359], w_fp[4], w_fp[486], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[183] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3525 OF 15495 ***
    // Wavefunction(s) for diagram number 3525
    // (none)
    // Amplitude(s) for diagram number 3525
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[454], w_fp[2], w_fp[359], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[183] += amp_sv[0];
    jamp_sv[207] -= amp_sv[0];
    jamp_sv[441] -= amp_sv[0];
    jamp_sv[561] += amp_sv[0];

    // *** DIAGRAM 3526 OF 15495 ***
    // Wavefunction(s) for diagram number 3526
    // (none)
    // Amplitude(s) for diagram number 3526
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[515], w_fp[244], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[447] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3527 OF 15495 ***
    // Wavefunction(s) for diagram number 3527
    // (none)
    // Amplitude(s) for diagram number 3527
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[454], w_fp[244], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3528 OF 15495 ***
    // Wavefunction(s) for diagram number 3528
    // (none)
    // Amplitude(s) for diagram number 3528
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[515], w_fp[2], w_fp[115], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[303] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];

    // *** DIAGRAM 3529 OF 15495 ***
    // Wavefunction(s) for diagram number 3529
    // (none)
    // Amplitude(s) for diagram number 3529
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[115], w_fp[486], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[159] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[183] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3530 OF 15495 ***
    // Wavefunction(s) for diagram number 3530
    // (none)
    // Amplitude(s) for diagram number 3530
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[2], w_fp[262], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[159] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[183] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[2], w_fp[353], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[183] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[2], w_fp[72], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[159] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] -= cxtype( 0, 1 ) * amp_sv[0];

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
  diagramgroup354( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 55 );
    retrieveWf( wfs, w_cx, nevt, 56 );
    retrieveWf( wfs, w_cx, nevt, 122 );
    retrieveWf( wfs, w_cx, nevt, 124 );
    retrieveWf( wfs, w_cx, nevt, 125 );
    retrieveWf( wfs, w_cx, nevt, 236 );
    retrieveWf( wfs, w_cx, nevt, 360 );
    retrieveWf( wfs, w_cx, nevt, 470 );
    retrieveWf( wfs, w_cx, nevt, 478 );
    retrieveWf( wfs, w_cx, nevt, 479 );
    retrieveWf( wfs, w_cx, nevt, 480 );
    retrieveWf( wfs, w_cx, nevt, 487 );
    retrieveWf( wfs, w_cx, nevt, 505 );
    retrieveWf( wfs, w_cx, nevt, 509 );
#endif
#endif

    // *** DIAGRAM 3531 OF 15495 ***
    // Wavefunction(s) for diagram number 3531
    // (none)
    // Amplitude(s) for diagram number 3531
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[478], w_fp[122], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[479] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3532 OF 15495 ***
    // Wavefunction(s) for diagram number 3532
    // (none)
    // Amplitude(s) for diagram number 3532
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[509], w_fp[122], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[473] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[689] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3533 OF 15495 ***
    // Wavefunction(s) for diagram number 3533
    // (none)
    // Amplitude(s) for diagram number 3533
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[124], w_fp[6], w_fp[470], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[311] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[689] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3534 OF 15495 ***
    // Wavefunction(s) for diagram number 3534
    // (none)
    // Amplitude(s) for diagram number 3534
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[509], w_fp[2], w_fp[124], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[311] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[689] += amp_sv[0];

    // *** DIAGRAM 3535 OF 15495 ***
    // Wavefunction(s) for diagram number 3535
    // (none)
    // Amplitude(s) for diagram number 3535
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[4], w_fp[125], w_fp[470], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[311] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3536 OF 15495 ***
    // Wavefunction(s) for diagram number 3536
    // (none)
    // Amplitude(s) for diagram number 3536
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[478], w_fp[2], w_fp[125], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[479] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    jamp_sv[695] -= amp_sv[0];

    // *** DIAGRAM 3537 OF 15495 ***
    // Wavefunction(s) for diagram number 3537
    // (none)
    // Amplitude(s) for diagram number 3537
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[2], w_fp[236], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[311] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[2], w_fp[56], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[329] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[689] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[2], w_fp[55], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[311] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[689] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3538 OF 15495 ***
    // Wavefunction(s) for diagram number 3538
    // (none)
    // Amplitude(s) for diagram number 3538
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[360], w_fp[4], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[231] += amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[323] -= amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[467] -= amp_sv[0];
    jamp_sv[537] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[681] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[360], w_fp[4], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[191] += amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[309] += amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[323] -= amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[467] -= amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[360], w_fp[4], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[189] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[309] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[465] -= amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[539] -= amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[681] += amp_sv[0];

    // *** DIAGRAM 3539 OF 15495 ***
    // Wavefunction(s) for diagram number 3539
    // (none)
    // Amplitude(s) for diagram number 3539
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[360], w_fp[6], w_fp[487], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[191] += amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[309] += amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[323] -= amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[467] -= amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];

    // *** DIAGRAM 3540 OF 15495 ***
    // Wavefunction(s) for diagram number 3540
    // (none)
    // Amplitude(s) for diagram number 3540
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[360], w_fp[4], w_fp[480], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[189] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[309] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[465] -= amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[539] -= amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[681] += amp_sv[0];

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
  diagramgroup355( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 46 );
    retrieveWf( wfs, w_cx, nevt, 55 );
    retrieveWf( wfs, w_cx, nevt, 56 );
    retrieveWf( wfs, w_cx, nevt, 81 );
    retrieveWf( wfs, w_cx, nevt, 124 );
    retrieveWf( wfs, w_cx, nevt, 125 );
    retrieveWf( wfs, w_cx, nevt, 236 );
    retrieveWf( wfs, w_cx, nevt, 261 );
    retrieveWf( wfs, w_cx, nevt, 278 );
    retrieveWf( wfs, w_cx, nevt, 352 );
    retrieveWf( wfs, w_cx, nevt, 358 );
    retrieveWf( wfs, w_cx, nevt, 360 );
    retrieveWf( wfs, w_cx, nevt, 469 );
    retrieveWf( wfs, w_cx, nevt, 479 );
    retrieveWf( wfs, w_cx, nevt, 480 );
    retrieveWf( wfs, w_cx, nevt, 487 );
    retrieveWf( wfs, w_cx, nevt, 495 );
#endif
#endif

    // *** DIAGRAM 3541 OF 15495 ***
    // Wavefunction(s) for diagram number 3541
    // (none)
    // Amplitude(s) for diagram number 3541
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[124], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    jamp_sv[201] += amp_sv[0];
    jamp_sv[203] -= amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];
    jamp_sv[689] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[124], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[201] += amp_sv[0];
    jamp_sv[203] -= amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[351] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[471] += amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[515] += amp_sv[0];
    jamp_sv[521] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    jamp_sv[689] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[124], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] += amp_sv[0];
    jamp_sv[165] -= amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[231] += amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[351] += amp_sv[0];
    jamp_sv[471] += amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[515] += amp_sv[0];
    jamp_sv[521] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];

    // *** DIAGRAM 3542 OF 15495 ***
    // Wavefunction(s) for diagram number 3542
    // (none)
    // Amplitude(s) for diagram number 3542
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[124], w_fp[6], w_fp[469], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    jamp_sv[201] += amp_sv[0];
    jamp_sv[203] -= amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];
    jamp_sv[689] += amp_sv[0];

    // *** DIAGRAM 3543 OF 15495 ***
    // Wavefunction(s) for diagram number 3543
    // (none)
    // Amplitude(s) for diagram number 3543
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[124], w_fp[480], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] += amp_sv[0];
    jamp_sv[165] -= amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[231] += amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[351] += amp_sv[0];
    jamp_sv[471] += amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[515] += amp_sv[0];
    jamp_sv[521] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];

    // *** DIAGRAM 3544 OF 15495 ***
    // Wavefunction(s) for diagram number 3544
    // (none)
    // Amplitude(s) for diagram number 3544
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[4], w_fp[125], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[159] += amp_sv[0];
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[4], w_fp[125], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[191] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[4], w_fp[125], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] += amp_sv[0];
    jamp_sv[159] -= amp_sv[0];
    jamp_sv[161] += amp_sv[0];
    jamp_sv[165] -= amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    jamp_sv[695] -= amp_sv[0];

    // *** DIAGRAM 3545 OF 15495 ***
    // Wavefunction(s) for diagram number 3545
    // (none)
    // Amplitude(s) for diagram number 3545
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[4], w_fp[125], w_fp[469], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[159] += amp_sv[0];
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];

    // *** DIAGRAM 3546 OF 15495 ***
    // Wavefunction(s) for diagram number 3546
    // (none)
    // Amplitude(s) for diagram number 3546
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[125], w_fp[487], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[191] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[693] += amp_sv[0];

    // *** DIAGRAM 3547 OF 15495 ***
    // Wavefunction(s) for diagram number 3547
    // (none)
    // Amplitude(s) for diagram number 3547
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[261], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[309] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    jamp_sv[513] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];
    jamp_sv[687] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[352], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[189] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[309] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[465] -= amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[539] -= amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[278], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] += amp_sv[0];
    jamp_sv[165] -= amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[465] -= amp_sv[0];
    jamp_sv[471] += amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[515] += amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[539] -= amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];

    // *** DIAGRAM 3548 OF 15495 ***
    // Wavefunction(s) for diagram number 3548
    // (none)
    // Amplitude(s) for diagram number 3548
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[358], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[209] += amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[233] += amp_sv[0];
    jamp_sv[275] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[281] += amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[353] += amp_sv[0];
    jamp_sv[477] += amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[597] += amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[81], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[209] += amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[281] += amp_sv[0];
    jamp_sv[309] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[323] -= amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[353] += amp_sv[0];
    jamp_sv[467] -= amp_sv[0];
    jamp_sv[477] += amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[46], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[191] += amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[309] += amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[323] -= amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[467] -= amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];

    // *** DIAGRAM 3549 OF 15495 ***
    // Wavefunction(s) for diagram number 3549
    // (none)
    // Amplitude(s) for diagram number 3549
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[236], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[159] += amp_sv[0];
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[56], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] += amp_sv[0];
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[201] -= amp_sv[0];
    jamp_sv[203] += amp_sv[0];
    jamp_sv[231] += amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[473] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[689] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[55], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] += amp_sv[0];
    jamp_sv[165] -= amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[201] -= amp_sv[0];
    jamp_sv[203] += amp_sv[0];
    jamp_sv[209] += amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[231] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[353] += amp_sv[0];
    jamp_sv[473] += amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    jamp_sv[689] -= amp_sv[0];

    // *** DIAGRAM 3550 OF 15495 ***
    // Wavefunction(s) for diagram number 3550
    // (none)
    // Amplitude(s) for diagram number 3550
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[360], w_fp[6], w_fp[495], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[191] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[467] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[521] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] -= cxtype( 0, 1 ) * amp_sv[0];

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
  diagramgroup356( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 46 );
    retrieveWf( wfs, w_cx, nevt, 81 );
    retrieveWf( wfs, w_cx, nevt, 122 );
    retrieveWf( wfs, w_cx, nevt, 125 );
    retrieveWf( wfs, w_cx, nevt, 358 );
    retrieveWf( wfs, w_cx, nevt, 360 );
    retrieveWf( wfs, w_cx, nevt, 435 );
    retrieveWf( wfs, w_cx, nevt, 437 );
    retrieveWf( wfs, w_cx, nevt, 450 );
    retrieveWf( wfs, w_cx, nevt, 455 );
    retrieveWf( wfs, w_cx, nevt, 485 );
    retrieveWf( wfs, w_cx, nevt, 495 );
    retrieveWf( wfs, w_cx, nevt, 514 );
#endif
#endif

    // *** DIAGRAM 3551 OF 15495 ***
    // Wavefunction(s) for diagram number 3551
    // (none)
    // Amplitude(s) for diagram number 3551
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[437], w_fp[2], w_fp[360], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[191] += amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[467] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];

    // *** DIAGRAM 3552 OF 15495 ***
    // Wavefunction(s) for diagram number 3552
    // (none)
    // Amplitude(s) for diagram number 3552
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[455], w_fp[122], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[477] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3553 OF 15495 ***
    // Wavefunction(s) for diagram number 3553
    // (none)
    // Amplitude(s) for diagram number 3553
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[437], w_fp[122], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[467] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3554 OF 15495 ***
    // Wavefunction(s) for diagram number 3554
    // (none)
    // Amplitude(s) for diagram number 3554
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[455], w_fp[2], w_fp[125], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[477] += amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[597] += amp_sv[0];
    jamp_sv[693] -= amp_sv[0];

    // *** DIAGRAM 3555 OF 15495 ***
    // Wavefunction(s) for diagram number 3555
    // (none)
    // Amplitude(s) for diagram number 3555
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[125], w_fp[495], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[191] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[209] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[477] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3556 OF 15495 ***
    // Wavefunction(s) for diagram number 3556
    // (none)
    // Amplitude(s) for diagram number 3556
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[435], w_fp[2], w_fp[358], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[191] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[209] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[477] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[435], w_fp[2], w_fp[81], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[209] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[467] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[477] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[521] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[435], w_fp[2], w_fp[46], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[191] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[467] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[521] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3557 OF 15495 ***
    // Wavefunction(s) for diagram number 3557
    // (none)
    // Amplitude(s) for diagram number 3557
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[360], w_fp[4], w_fp[485], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[189] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[309] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[351] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3558 OF 15495 ***
    // Wavefunction(s) for diagram number 3558
    // (none)
    // Amplitude(s) for diagram number 3558
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[450], w_fp[2], w_fp[360], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[189] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[465] -= amp_sv[0];
    jamp_sv[681] += amp_sv[0];

    // *** DIAGRAM 3559 OF 15495 ***
    // Wavefunction(s) for diagram number 3559
    // (none)
    // Amplitude(s) for diagram number 3559
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[514], w_fp[122], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[471] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[687] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3560 OF 15495 ***
    // Wavefunction(s) for diagram number 3560
    // (none)
    // Amplitude(s) for diagram number 3560
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[450], w_fp[122], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[465] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] -= cxtype( 0, 1 ) * amp_sv[0];

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
  diagramgroup357( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 124 );
    retrieveWf( wfs, w_cx, nevt, 128 );
    retrieveWf( wfs, w_cx, nevt, 130 );
    retrieveWf( wfs, w_cx, nevt, 131 );
    retrieveWf( wfs, w_cx, nevt, 138 );
    retrieveWf( wfs, w_cx, nevt, 155 );
    retrieveWf( wfs, w_cx, nevt, 238 );
    retrieveWf( wfs, w_cx, nevt, 261 );
    retrieveWf( wfs, w_cx, nevt, 278 );
    retrieveWf( wfs, w_cx, nevt, 352 );
    retrieveWf( wfs, w_cx, nevt, 447 );
    retrieveWf( wfs, w_cx, nevt, 470 );
    retrieveWf( wfs, w_cx, nevt, 478 );
    retrieveWf( wfs, w_cx, nevt, 485 );
    retrieveWf( wfs, w_cx, nevt, 505 );
    retrieveWf( wfs, w_cx, nevt, 511 );
    retrieveWf( wfs, w_cx, nevt, 514 );
#endif
#endif

    // *** DIAGRAM 3561 OF 15495 ***
    // Wavefunction(s) for diagram number 3561
    // (none)
    // Amplitude(s) for diagram number 3561
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[514], w_fp[2], w_fp[124], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[309] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    jamp_sv[687] += amp_sv[0];

    // *** DIAGRAM 3562 OF 15495 ***
    // Wavefunction(s) for diagram number 3562
    // (none)
    // Amplitude(s) for diagram number 3562
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[124], w_fp[485], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[165] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[309] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[351] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[471] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[687] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3563 OF 15495 ***
    // Wavefunction(s) for diagram number 3563
    // (none)
    // Amplitude(s) for diagram number 3563
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[447], w_fp[2], w_fp[261], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[165] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[309] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[351] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[471] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[687] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[447], w_fp[2], w_fp[352], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[189] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[309] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[351] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[447], w_fp[2], w_fp[278], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[165] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[471] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[687] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3564 OF 15495 ***
    // Wavefunction(s) for diagram number 3564
    // (none)
    // Amplitude(s) for diagram number 3564
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[478], w_fp[128], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[599] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3565 OF 15495 ***
    // Wavefunction(s) for diagram number 3565
    // (none)
    // Amplitude(s) for diagram number 3565
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[511], w_fp[128], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[593] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[713] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3566 OF 15495 ***
    // Wavefunction(s) for diagram number 3566
    // (none)
    // Amplitude(s) for diagram number 3566
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[130], w_fp[5], w_fp[470], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[335] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[713] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3567 OF 15495 ***
    // Wavefunction(s) for diagram number 3567
    // (none)
    // Amplitude(s) for diagram number 3567
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[511], w_fp[2], w_fp[130], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[335] += amp_sv[0];
    jamp_sv[359] -= amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[713] += amp_sv[0];

    // *** DIAGRAM 3568 OF 15495 ***
    // Wavefunction(s) for diagram number 3568
    // (none)
    // Amplitude(s) for diagram number 3568
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[4], w_fp[131], w_fp[470], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[305] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[311] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3569 OF 15495 ***
    // Wavefunction(s) for diagram number 3569
    // (none)
    // Amplitude(s) for diagram number 3569
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[478], w_fp[2], w_fp[131], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[455] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    jamp_sv[599] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];

    // *** DIAGRAM 3570 OF 15495 ***
    // Wavefunction(s) for diagram number 3570
    // (none)
    // Amplitude(s) for diagram number 3570
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[2], w_fp[155], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[305] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[311] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[2], w_fp[138], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[335] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[713] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[2], w_fp[238], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[305] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[311] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[713] += cxtype( 0, 1 ) * amp_sv[0];
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
  diagramgroup358( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 130 );
    retrieveWf( wfs, w_cx, nevt, 131 );
    retrieveWf( wfs, w_cx, nevt, 248 );
    retrieveWf( wfs, w_cx, nevt, 264 );
    retrieveWf( wfs, w_cx, nevt, 291 );
    retrieveWf( wfs, w_cx, nevt, 361 );
    retrieveWf( wfs, w_cx, nevt, 469 );
    retrieveWf( wfs, w_cx, nevt, 479 );
    retrieveWf( wfs, w_cx, nevt, 482 );
    retrieveWf( wfs, w_cx, nevt, 487 );
#endif
#endif

    // *** DIAGRAM 3571 OF 15495 ***
    // Wavefunction(s) for diagram number 3571
    // (none)
    // Amplitude(s) for diagram number 3571
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[361], w_fp[4], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[237] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[299] -= amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[417] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[587] -= amp_sv[0];
    jamp_sv[705] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[361], w_fp[4], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[215] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[299] -= amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[333] += amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];
    jamp_sv[453] += amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[587] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[361], w_fp[4], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[213] += amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];
    jamp_sv[333] += amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[419] -= amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[453] += amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];

    // *** DIAGRAM 3572 OF 15495 ***
    // Wavefunction(s) for diagram number 3572
    // (none)
    // Amplitude(s) for diagram number 3572
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[361], w_fp[5], w_fp[487], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[215] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[299] -= amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[333] += amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];
    jamp_sv[453] += amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[587] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];

    // *** DIAGRAM 3573 OF 15495 ***
    // Wavefunction(s) for diagram number 3573
    // (none)
    // Amplitude(s) for diagram number 3573
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[361], w_fp[4], w_fp[482], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[213] += amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];
    jamp_sv[333] += amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[419] -= amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[453] += amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];

    // *** DIAGRAM 3574 OF 15495 ***
    // Wavefunction(s) for diagram number 3574
    // (none)
    // Amplitude(s) for diagram number 3574
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[130], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[167] += amp_sv[0];
    jamp_sv[177] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[213] += amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[359] -= amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[713] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[130], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[177] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[357] += amp_sv[0];
    jamp_sv[359] -= amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[395] += amp_sv[0];
    jamp_sv[401] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[591] += amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[711] -= amp_sv[0];
    jamp_sv[713] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[130], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] += amp_sv[0];
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[237] += amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[357] += amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[395] += amp_sv[0];
    jamp_sv[401] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    jamp_sv[591] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];

    // *** DIAGRAM 3575 OF 15495 ***
    // Wavefunction(s) for diagram number 3575
    // (none)
    // Amplitude(s) for diagram number 3575
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[130], w_fp[5], w_fp[469], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[167] += amp_sv[0];
    jamp_sv[177] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[213] += amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[359] -= amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[713] += amp_sv[0];

    // *** DIAGRAM 3576 OF 15495 ***
    // Wavefunction(s) for diagram number 3576
    // (none)
    // Amplitude(s) for diagram number 3576
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[130], w_fp[482], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] += amp_sv[0];
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[237] += amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[357] += amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[395] += amp_sv[0];
    jamp_sv[401] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    jamp_sv[591] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];

    // *** DIAGRAM 3577 OF 15495 ***
    // Wavefunction(s) for diagram number 3577
    // (none)
    // Amplitude(s) for diagram number 3577
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[4], w_fp[131], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[155] += amp_sv[0];
    jamp_sv[161] += amp_sv[0];
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[185] += amp_sv[0];
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[4], w_fp[131], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] += amp_sv[0];
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[275] += amp_sv[0];
    jamp_sv[281] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[477] += amp_sv[0];
    jamp_sv[597] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[4], w_fp[131], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += amp_sv[0];
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[167] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[275] += amp_sv[0];
    jamp_sv[281] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[477] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    jamp_sv[597] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];

    // *** DIAGRAM 3578 OF 15495 ***
    // Wavefunction(s) for diagram number 3578
    // (none)
    // Amplitude(s) for diagram number 3578
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[4], w_fp[131], w_fp[469], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[155] += amp_sv[0];
    jamp_sv[161] += amp_sv[0];
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[185] += amp_sv[0];
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];

    // *** DIAGRAM 3579 OF 15495 ***
    // Wavefunction(s) for diagram number 3579
    // (none)
    // Amplitude(s) for diagram number 3579
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[131], w_fp[487], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] += amp_sv[0];
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[275] += amp_sv[0];
    jamp_sv[281] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[477] += amp_sv[0];
    jamp_sv[597] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];

    // *** DIAGRAM 3580 OF 15495 ***
    // Wavefunction(s) for diagram number 3580
    // (none)
    // Amplitude(s) for diagram number 3580
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[291], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[167] += amp_sv[0];
    jamp_sv[213] += amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[333] += amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[393] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[711] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[248], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[213] += amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];
    jamp_sv[333] += amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[419] -= amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[453] += amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[264], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] += amp_sv[0];
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[395] += amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[419] -= amp_sv[0];
    jamp_sv[453] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[591] += amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];

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
  diagramgroup359( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 128 );
    retrieveWf( wfs, w_cx, nevt, 131 );
    retrieveWf( wfs, w_cx, nevt, 138 );
    retrieveWf( wfs, w_cx, nevt, 155 );
    retrieveWf( wfs, w_cx, nevt, 238 );
    retrieveWf( wfs, w_cx, nevt, 249 );
    retrieveWf( wfs, w_cx, nevt, 257 );
    retrieveWf( wfs, w_cx, nevt, 354 );
    retrieveWf( wfs, w_cx, nevt, 361 );
    retrieveWf( wfs, w_cx, nevt, 435 );
    retrieveWf( wfs, w_cx, nevt, 440 );
    retrieveWf( wfs, w_cx, nevt, 455 );
    retrieveWf( wfs, w_cx, nevt, 479 );
    retrieveWf( wfs, w_cx, nevt, 484 );
    retrieveWf( wfs, w_cx, nevt, 495 );
#endif
#endif

    // *** DIAGRAM 3581 OF 15495 ***
    // Wavefunction(s) for diagram number 3581
    // (none)
    // Amplitude(s) for diagram number 3581
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[257], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[359] -= amp_sv[0];
    jamp_sv[453] += amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[717] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[249], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[215] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[299] -= amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[333] += amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];
    jamp_sv[453] += amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[587] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[354], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] += amp_sv[0];
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[275] += amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[299] -= amp_sv[0];
    jamp_sv[333] += amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];
    jamp_sv[587] -= amp_sv[0];
    jamp_sv[597] += amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];

    // *** DIAGRAM 3582 OF 15495 ***
    // Wavefunction(s) for diagram number 3582
    // (none)
    // Amplitude(s) for diagram number 3582
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[155], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[155] += amp_sv[0];
    jamp_sv[161] += amp_sv[0];
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[185] += amp_sv[0];
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[138], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] += amp_sv[0];
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[177] -= amp_sv[0];
    jamp_sv[179] += amp_sv[0];
    jamp_sv[185] += amp_sv[0];
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[237] += amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    jamp_sv[593] += amp_sv[0];
    jamp_sv[713] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[1], w_fp[238], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += amp_sv[0];
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[177] -= amp_sv[0];
    jamp_sv[179] += amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[237] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];
    jamp_sv[593] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];
    jamp_sv[713] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];

    // *** DIAGRAM 3583 OF 15495 ***
    // Wavefunction(s) for diagram number 3583
    // (none)
    // Amplitude(s) for diagram number 3583
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[361], w_fp[5], w_fp[495], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[215] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[401] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[453] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[477] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3584 OF 15495 ***
    // Wavefunction(s) for diagram number 3584
    // (none)
    // Amplitude(s) for diagram number 3584
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[440], w_fp[2], w_fp[361], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[215] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];
    jamp_sv[587] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];

    // *** DIAGRAM 3585 OF 15495 ***
    // Wavefunction(s) for diagram number 3585
    // (none)
    // Amplitude(s) for diagram number 3585
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[455], w_fp[128], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[597] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3586 OF 15495 ***
    // Wavefunction(s) for diagram number 3586
    // (none)
    // Amplitude(s) for diagram number 3586
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[440], w_fp[128], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[587] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3587 OF 15495 ***
    // Wavefunction(s) for diagram number 3587
    // (none)
    // Amplitude(s) for diagram number 3587
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[455], w_fp[2], w_fp[131], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[453] += amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[717] += amp_sv[0];

    // *** DIAGRAM 3588 OF 15495 ***
    // Wavefunction(s) for diagram number 3588
    // (none)
    // Amplitude(s) for diagram number 3588
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[131], w_fp[495], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[453] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[477] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3589 OF 15495 ***
    // Wavefunction(s) for diagram number 3589
    // (none)
    // Amplitude(s) for diagram number 3589
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[435], w_fp[2], w_fp[257], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[453] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[477] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[435], w_fp[2], w_fp[249], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[401] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[453] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[477] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[435], w_fp[2], w_fp[354], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[401] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3590 OF 15495 ***
    // Wavefunction(s) for diagram number 3590
    // (none)
    // Amplitude(s) for diagram number 3590
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[361], w_fp[4], w_fp[484], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];

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
  diagramgroup360( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 74 );
    retrieveWf( wfs, w_cx, nevt, 76 );
    retrieveWf( wfs, w_cx, nevt, 128 );
    retrieveWf( wfs, w_cx, nevt, 130 );
    retrieveWf( wfs, w_cx, nevt, 142 );
    retrieveWf( wfs, w_cx, nevt, 172 );
    retrieveWf( wfs, w_cx, nevt, 248 );
    retrieveWf( wfs, w_cx, nevt, 264 );
    retrieveWf( wfs, w_cx, nevt, 276 );
    retrieveWf( wfs, w_cx, nevt, 277 );
    retrieveWf( wfs, w_cx, nevt, 291 );
    retrieveWf( wfs, w_cx, nevt, 303 );
    retrieveWf( wfs, w_cx, nevt, 356 );
    retrieveWf( wfs, w_cx, nevt, 357 );
    retrieveWf( wfs, w_cx, nevt, 361 );
    retrieveWf( wfs, w_cx, nevt, 442 );
    retrieveWf( wfs, w_cx, nevt, 446 );
    retrieveWf( wfs, w_cx, nevt, 447 );
    retrieveWf( wfs, w_cx, nevt, 479 );
    retrieveWf( wfs, w_cx, nevt, 484 );
    retrieveWf( wfs, w_cx, nevt, 513 );
#endif
#endif

    // *** DIAGRAM 3591 OF 15495 ***
    // Wavefunction(s) for diagram number 3591
    // (none)
    // Amplitude(s) for diagram number 3591
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[446], w_fp[2], w_fp[361], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[213] += amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];

    // *** DIAGRAM 3592 OF 15495 ***
    // Wavefunction(s) for diagram number 3592
    // (none)
    // Amplitude(s) for diagram number 3592
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[513], w_fp[128], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[591] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3593 OF 15495 ***
    // Wavefunction(s) for diagram number 3593
    // (none)
    // Amplitude(s) for diagram number 3593
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[446], w_fp[128], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3594 OF 15495 ***
    // Wavefunction(s) for diagram number 3594
    // (none)
    // Amplitude(s) for diagram number 3594
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[513], w_fp[2], w_fp[130], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[333] += amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[711] += amp_sv[0];

    // *** DIAGRAM 3595 OF 15495 ***
    // Wavefunction(s) for diagram number 3595
    // (none)
    // Amplitude(s) for diagram number 3595
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[130], w_fp[484], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3596 OF 15495 ***
    // Wavefunction(s) for diagram number 3596
    // (none)
    // Amplitude(s) for diagram number 3596
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[442], w_fp[2], w_fp[291], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[442], w_fp[2], w_fp[248], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[442], w_fp[2], w_fp[264], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3597 OF 15495 ***
    // Wavefunction(s) for diagram number 3597
    // (none)
    // Amplitude(s) for diagram number 3597
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[303], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[177] += amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[417] -= amp_sv[0];
    jamp_sv[513] += amp_sv[0];
    jamp_sv[519] -= amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[633] += amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    jamp_sv[663] -= amp_sv[0];
    jamp_sv[687] += amp_sv[0];
    jamp_sv[705] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    jamp_sv[713] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[172], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[177] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[519] -= amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    jamp_sv[663] -= amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];
    jamp_sv[713] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[76], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
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

    // *** DIAGRAM 3598 OF 15495 ***
    // Wavefunction(s) for diagram number 3598
    // (none)
    // Amplitude(s) for diagram number 3598
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[277], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[179] += amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[419] -= amp_sv[0];
    jamp_sv[513] += amp_sv[0];
    jamp_sv[519] -= amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[593] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];
    jamp_sv[633] += amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    jamp_sv[663] -= amp_sv[0];
    jamp_sv[687] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[357], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[179] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[519] -= amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[593] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    jamp_sv[663] -= amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[74], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
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

    // *** DIAGRAM 3599 OF 15495 ***
    // Wavefunction(s) for diagram number 3599
    // (none)
    // Amplitude(s) for diagram number 3599
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[276], w_fp[6], w_fp[7], w_fp[479], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[177] += amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[417] -= amp_sv[0];
    jamp_sv[513] += amp_sv[0];
    jamp_sv[519] -= amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[633] += amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    jamp_sv[663] -= amp_sv[0];
    jamp_sv[687] += amp_sv[0];
    jamp_sv[705] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    jamp_sv[713] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[276], w_fp[6], w_fp[7], w_fp[479], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[179] += amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[419] -= amp_sv[0];
    jamp_sv[513] += amp_sv[0];
    jamp_sv[519] -= amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[593] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];
    jamp_sv[633] += amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    jamp_sv[663] -= amp_sv[0];
    jamp_sv[687] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[276], w_fp[6], w_fp[7], w_fp[479], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += amp_sv[0];
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[177] -= amp_sv[0];
    jamp_sv[179] += amp_sv[0];
    jamp_sv[297] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[419] -= amp_sv[0];
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[593] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[707] -= amp_sv[0];
    jamp_sv[713] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[356], w_fp[6], w_fp[7], w_fp[479], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[177] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[519] -= amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    jamp_sv[663] -= amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];
    jamp_sv[713] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[356], w_fp[6], w_fp[7], w_fp[479], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[179] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[519] -= amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[593] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    jamp_sv[663] -= amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[356], w_fp[6], w_fp[7], w_fp[479], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[177] -= amp_sv[0];
    jamp_sv[179] += amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[297] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[393] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[593] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[707] -= amp_sv[0];
    jamp_sv[711] += amp_sv[0];
    jamp_sv[713] -= amp_sv[0];
    jamp_sv[717] += amp_sv[0];
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[142], w_fp[6], w_fp[7], w_fp[479], COUPs[2], 1.0, &amp_fp[0] );
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
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[142], w_fp[6], w_fp[7], w_fp[479], COUPs[2], 1.0, &amp_fp[0] );
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
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[142], w_fp[6], w_fp[7], w_fp[479], COUPs[2], 1.0, &amp_fp[0] );
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

    // *** DIAGRAM 3600 OF 15495 ***
    // Wavefunction(s) for diagram number 3600
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[276], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[484] );
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[356], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[513] );
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[142], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[446] );
    // Amplitude(s) for diagram number 3600
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[447], w_fp[484], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[299] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[447], w_fp[513], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[275] += amp_sv[0];
    jamp_sv[299] -= amp_sv[0];
    jamp_sv[395] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[447], w_fp[446], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[275] += amp_sv[0];
    jamp_sv[395] += amp_sv[0];
    jamp_sv[419] -= amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 446 );
    storeWf( wfs, w_cx, nevt, 484 );
    storeWf( wfs, w_cx, nevt, 513 );
#endif
#endif
  }

  //--------------------------------------------------------------------------

}
