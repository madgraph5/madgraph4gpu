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
  diagramgroup271( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 45 );
    retrieveWf( wfs, w_cx, nevt, 100 );
    retrieveWf( wfs, w_cx, nevt, 122 );
    retrieveWf( wfs, w_cx, nevt, 125 );
    retrieveWf( wfs, w_cx, nevt, 258 );
    retrieveWf( wfs, w_cx, nevt, 263 );
    retrieveWf( wfs, w_cx, nevt, 287 );
    retrieveWf( wfs, w_cx, nevt, 447 );
    retrieveWf( wfs, w_cx, nevt, 461 );
    retrieveWf( wfs, w_cx, nevt, 479 );
    retrieveWf( wfs, w_cx, nevt, 483 );
#endif
#endif

    // *** DIAGRAM 2701 OF 15495 ***
    // Wavefunction(s) for diagram number 2701
    // (none)
    // Amplitude(s) for diagram number 2701
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[263], w_fp[100], COUPs[0], 1.0, &amp_fp[0] );
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

    // *** DIAGRAM 2702 OF 15495 ***
    // Wavefunction(s) for diagram number 2702
    // (none)
    // Amplitude(s) for diagram number 2702
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[258], w_fp[125], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[159] += amp_sv[0];
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[275] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[281] += amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[477] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[597] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];

    // *** DIAGRAM 2703 OF 15495 ***
    // Wavefunction(s) for diagram number 2703
    // (none)
    // Amplitude(s) for diagram number 2703
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[258], w_fp[100], w_fp[6], w_fp[479], COUPs[2], 1.0, &amp_fp[0] );
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
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[258], w_fp[100], w_fp[6], w_fp[479], COUPs[2], 1.0, &amp_fp[0] );
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
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[258], w_fp[100], w_fp[6], w_fp[479], COUPs[2], 1.0, &amp_fp[0] );
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

    // *** DIAGRAM 2704 OF 15495 ***
    // Wavefunction(s) for diagram number 2704
    // (none)
    // Amplitude(s) for diagram number 2704
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[461], w_fp[483], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] += amp_sv[0];
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[281] += amp_sv[0];

    // *** DIAGRAM 2705 OF 15495 ***
    // Wavefunction(s) for diagram number 2705
    // (none)
    // Amplitude(s) for diagram number 2705
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[461], w_fp[2], w_fp[263], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2706 OF 15495 ***
    // Wavefunction(s) for diagram number 2706
    // (none)
    // Amplitude(s) for diagram number 2706
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[447], w_fp[483], w_fp[100], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] += amp_sv[0];
    jamp_sv[165] -= amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[285] += amp_sv[0];

    // *** DIAGRAM 2707 OF 15495 ***
    // Wavefunction(s) for diagram number 2707
    // (none)
    // Amplitude(s) for diagram number 2707
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[447], w_fp[2], w_fp[287], COUPs[1], 1.0, &amp_fp[0] );
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

    // *** DIAGRAM 2708 OF 15495 ***
    // Wavefunction(s) for diagram number 2708
    // (none)
    // Amplitude(s) for diagram number 2708
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[447], w_fp[122], w_fp[258], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[465] += amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    jamp_sv[681] -= amp_sv[0];
    jamp_sv[687] += amp_sv[0];

    // *** DIAGRAM 2709 OF 15495 ***
    // Wavefunction(s) for diagram number 2709
    // (none)
    // Amplitude(s) for diagram number 2709
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[45], w_fp[483], w_fp[125], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[159] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[165] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2710 OF 15495 ***
    // Wavefunction(s) for diagram number 2710
    // (none)
    // Amplitude(s) for diagram number 2710
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[45], w_fp[122], w_fp[263], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[465] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[471] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[477] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[687] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] += cxtype( 0, 1 ) * amp_sv[0];

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
  diagramgroup272( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 84 );
    retrieveWf( wfs, w_cx, nevt, 128 );
    retrieveWf( wfs, w_cx, nevt, 131 );
    retrieveWf( wfs, w_cx, nevt, 258 );
    retrieveWf( wfs, w_cx, nevt, 260 );
    retrieveWf( wfs, w_cx, nevt, 288 );
    retrieveWf( wfs, w_cx, nevt, 442 );
    retrieveWf( wfs, w_cx, nevt, 458 );
    retrieveWf( wfs, w_cx, nevt, 479 );
    retrieveWf( wfs, w_cx, nevt, 483 );
    retrieveWf( wfs, w_cx, nevt, 487 );
#endif
#endif

    // *** DIAGRAM 2711 OF 15495 ***
    // Wavefunction(s) for diagram number 2711
    // (none)
    // Amplitude(s) for diagram number 2711
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[487], w_fp[128], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[597] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];

    // *** DIAGRAM 2712 OF 15495 ***
    // Wavefunction(s) for diagram number 2712
    // (none)
    // Amplitude(s) for diagram number 2712
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[487], w_fp[2], w_fp[131], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[453] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[477] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2713 OF 15495 ***
    // Wavefunction(s) for diagram number 2713
    // (none)
    // Amplitude(s) for diagram number 2713
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[260], w_fp[84], COUPs[0], 1.0, &amp_fp[0] );
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

    // *** DIAGRAM 2714 OF 15495 ***
    // Wavefunction(s) for diagram number 2714
    // (none)
    // Amplitude(s) for diagram number 2714
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[288], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
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

    // *** DIAGRAM 2715 OF 15495 ***
    // Wavefunction(s) for diagram number 2715
    // (none)
    // Amplitude(s) for diagram number 2715
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[258], w_fp[131], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[155] += amp_sv[0];
    jamp_sv[161] += amp_sv[0];
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];
    jamp_sv[453] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    jamp_sv[717] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];

    // *** DIAGRAM 2716 OF 15495 ***
    // Wavefunction(s) for diagram number 2716
    // (none)
    // Amplitude(s) for diagram number 2716
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[258], w_fp[5], w_fp[84], w_fp[479], COUPs[2], 1.0, &amp_fp[0] );
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
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[258], w_fp[5], w_fp[84], w_fp[479], COUPs[2], 1.0, &amp_fp[0] );
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
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[258], w_fp[5], w_fp[84], w_fp[479], COUPs[2], 1.0, &amp_fp[0] );
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

    // *** DIAGRAM 2717 OF 15495 ***
    // Wavefunction(s) for diagram number 2717
    // (none)
    // Amplitude(s) for diagram number 2717
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[442], w_fp[483], w_fp[84], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] += amp_sv[0];
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];

    // *** DIAGRAM 2718 OF 15495 ***
    // Wavefunction(s) for diagram number 2718
    // (none)
    // Amplitude(s) for diagram number 2718
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[442], w_fp[2], w_fp[288], COUPs[1], 1.0, &amp_fp[0] );
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

    // *** DIAGRAM 2719 OF 15495 ***
    // Wavefunction(s) for diagram number 2719
    // (none)
    // Amplitude(s) for diagram number 2719
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[442], w_fp[128], w_fp[258], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[585] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[705] -= amp_sv[0];
    jamp_sv[711] += amp_sv[0];

    // *** DIAGRAM 2720 OF 15495 ***
    // Wavefunction(s) for diagram number 2720
    // (none)
    // Amplitude(s) for diagram number 2720
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[458], w_fp[483], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += amp_sv[0];
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[275] += amp_sv[0];

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
  diagramgroup273( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 45 );
    retrieveWf( wfs, w_cx, nevt, 128 );
    retrieveWf( wfs, w_cx, nevt, 131 );
    retrieveWf( wfs, w_cx, nevt, 151 );
    retrieveWf( wfs, w_cx, nevt, 152 );
    retrieveWf( wfs, w_cx, nevt, 153 );
    retrieveWf( wfs, w_cx, nevt, 258 );
    retrieveWf( wfs, w_cx, nevt, 260 );
    retrieveWf( wfs, w_cx, nevt, 289 );
    retrieveWf( wfs, w_cx, nevt, 290 );
    retrieveWf( wfs, w_cx, nevt, 455 );
    retrieveWf( wfs, w_cx, nevt, 458 );
    retrieveWf( wfs, w_cx, nevt, 477 );
    retrieveWf( wfs, w_cx, nevt, 478 );
    retrieveWf( wfs, w_cx, nevt, 479 );
    retrieveWf( wfs, w_cx, nevt, 480 );
    retrieveWf( wfs, w_cx, nevt, 481 );
    retrieveWf( wfs, w_cx, nevt, 483 );
    retrieveWf( wfs, w_cx, nevt, 487 );
#endif
#endif

    // *** DIAGRAM 2721 OF 15495 ***
    // Wavefunction(s) for diagram number 2721
    // (none)
    // Amplitude(s) for diagram number 2721
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[458], w_fp[2], w_fp[260], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[155] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[419] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2722 OF 15495 ***
    // Wavefunction(s) for diagram number 2722
    // (none)
    // Amplitude(s) for diagram number 2722
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[45], w_fp[483], w_fp[131], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[155] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2723 OF 15495 ***
    // Wavefunction(s) for diagram number 2723
    // (none)
    // Amplitude(s) for diagram number 2723
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[45], w_fp[128], w_fp[260], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[585] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2724 OF 15495 ***
    // Wavefunction(s) for diagram number 2724
    // (none)
    // Amplitude(s) for diagram number 2724
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[487], w_fp[2], w_fp[151], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[453] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[477] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[487], w_fp[2], w_fp[152], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[477] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[487], w_fp[2], w_fp[153], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[453] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2725 OF 15495 ***
    // Wavefunction(s) for diagram number 2725
    // (none)
    // Amplitude(s) for diagram number 2725
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[258], w_fp[151], w_fp[479], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[155] += amp_sv[0];
    jamp_sv[161] += amp_sv[0];
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];
    jamp_sv[453] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    jamp_sv[717] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[258], w_fp[152], w_fp[479], COUPs[0], 1.0, &amp_fp[0] );
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
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[258], w_fp[153], w_fp[479], COUPs[0], 1.0, &amp_fp[0] );
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

    // *** DIAGRAM 2726 OF 15495 ***
    // Wavefunction(s) for diagram number 2726
    // (none)
    // Amplitude(s) for diagram number 2726
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[455], w_fp[2], w_fp[258], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[155] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[477], w_fp[2], w_fp[258], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[159] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[165] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[478], w_fp[2], w_fp[258], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[159] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[165] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2727 OF 15495 ***
    // Wavefunction(s) for diagram number 2727
    // (none)
    // Amplitude(s) for diagram number 2727
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[289], w_fp[6], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
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
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[289], w_fp[6], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
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
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[289], w_fp[6], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
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

    // *** DIAGRAM 2728 OF 15495 ***
    // Wavefunction(s) for diagram number 2728
    // (none)
    // Amplitude(s) for diagram number 2728
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[289], w_fp[7], w_fp[480], COUPs[0], 1.0, &amp_fp[0] );
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

    // *** DIAGRAM 2729 OF 15495 ***
    // Wavefunction(s) for diagram number 2729
    // (none)
    // Amplitude(s) for diagram number 2729
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[289], w_fp[6], w_fp[481], COUPs[0], 1.0, &amp_fp[0] );
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

    // *** DIAGRAM 2730 OF 15495 ***
    // Wavefunction(s) for diagram number 2730
    // (none)
    // Amplitude(s) for diagram number 2730
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[290], w_fp[4], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[183] -= amp_sv[0];
    jamp_sv[185] += amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[561] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[657] -= amp_sv[0];
    jamp_sv[663] += amp_sv[0];
    jamp_sv[669] += amp_sv[0];
    jamp_sv[671] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[290], w_fp[4], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[290], w_fp[4], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[183] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[519] -= amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    jamp_sv[663] -= amp_sv[0];
    jamp_sv[669] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];

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
  diagramgroup274( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 279 );
    retrieveWf( wfs, w_cx, nevt, 290 );
    retrieveWf( wfs, w_cx, nevt, 292 );
    retrieveWf( wfs, w_cx, nevt, 293 );
    retrieveWf( wfs, w_cx, nevt, 294 );
    retrieveWf( wfs, w_cx, nevt, 295 );
    retrieveWf( wfs, w_cx, nevt, 296 );
    retrieveWf( wfs, w_cx, nevt, 297 );
    retrieveWf( wfs, w_cx, nevt, 298 );
    retrieveWf( wfs, w_cx, nevt, 299 );
    retrieveWf( wfs, w_cx, nevt, 300 );
    retrieveWf( wfs, w_cx, nevt, 301 );
    retrieveWf( wfs, w_cx, nevt, 437 );
    retrieveWf( wfs, w_cx, nevt, 438 );
    retrieveWf( wfs, w_cx, nevt, 479 );
    retrieveWf( wfs, w_cx, nevt, 480 );
    retrieveWf( wfs, w_cx, nevt, 481 );
#endif
#endif

    // *** DIAGRAM 2731 OF 15495 ***
    // Wavefunction(s) for diagram number 2731
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[479], w_fp[4], COUPs[0], 1.0, 0., 0., w_fp[487] );
    // Amplitude(s) for diagram number 2731
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[290], w_fp[7], w_fp[487], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];

    // *** DIAGRAM 2732 OF 15495 ***
    // Wavefunction(s) for diagram number 2732
    // (none)
    // Amplitude(s) for diagram number 2732
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[290], w_fp[4], w_fp[481], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[183] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[519] -= amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    jamp_sv[663] -= amp_sv[0];
    jamp_sv[669] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];

    // *** DIAGRAM 2733 OF 15495 ***
    // Wavefunction(s) for diagram number 2733
    // (none)
    // Amplitude(s) for diagram number 2733
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[292], w_fp[4], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[537] -= amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[549] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[681] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[292], w_fp[4], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[191] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[292], w_fp[4], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[189] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[549] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    jamp_sv[681] += amp_sv[0];

    // *** DIAGRAM 2734 OF 15495 ***
    // Wavefunction(s) for diagram number 2734
    // (none)
    // Amplitude(s) for diagram number 2734
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[292], w_fp[6], w_fp[487], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[191] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];

    // *** DIAGRAM 2735 OF 15495 ***
    // Wavefunction(s) for diagram number 2735
    // (none)
    // Amplitude(s) for diagram number 2735
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[292], w_fp[4], w_fp[480], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[189] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[549] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    jamp_sv[681] += amp_sv[0];

    // *** DIAGRAM 2736 OF 15495 ***
    // Wavefunction(s) for diagram number 2736
    // (none)
    // Amplitude(s) for diagram number 2736
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[293], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[177] -= amp_sv[0];
    jamp_sv[183] += amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[393] += amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[537] -= amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[669] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    jamp_sv[681] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    jamp_sv[711] += amp_sv[0];
    jamp_sv[713] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[294], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[183] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[519] -= amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    jamp_sv[663] -= amp_sv[0];
    jamp_sv[669] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[295], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
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

    // *** DIAGRAM 2737 OF 15495 ***
    // Wavefunction(s) for diagram number 2737
    // (none)
    // Amplitude(s) for diagram number 2737
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[296], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[395] += amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[549] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[561] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[591] += amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[657] -= amp_sv[0];
    jamp_sv[663] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[297], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[189] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[549] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[298], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
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

    // *** DIAGRAM 2738 OF 15495 ***
    // Wavefunction(s) for diagram number 2738
    // (none)
    // Amplitude(s) for diagram number 2738
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[299], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[297] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[357] += amp_sv[0];
    jamp_sv[359] -= amp_sv[0];
    jamp_sv[401] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[707] -= amp_sv[0];
    jamp_sv[717] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[300], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[191] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[301], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];

    // *** DIAGRAM 2739 OF 15495 ***
    // Wavefunction(s) for diagram number 2739
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[279], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[494] );
    // Amplitude(s) for diagram number 2739
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[437], w_fp[494], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[191] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2740 OF 15495 ***
    // Wavefunction(s) for diagram number 2740
    // (none)
    // Amplitude(s) for diagram number 2740
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[438], w_fp[494], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[401] -= cxtype( 0, 1 ) * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 487 );
    storeWf( wfs, w_cx, nevt, 494 );
#endif
#endif
  }

  //--------------------------------------------------------------------------


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup275( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 289 );
    retrieveWf( wfs, w_cx, nevt, 290 );
    retrieveWf( wfs, w_cx, nevt, 292 );
    retrieveWf( wfs, w_cx, nevt, 299 );
    retrieveWf( wfs, w_cx, nevt, 300 );
    retrieveWf( wfs, w_cx, nevt, 301 );
    retrieveWf( wfs, w_cx, nevt, 435 );
    retrieveWf( wfs, w_cx, nevt, 437 );
    retrieveWf( wfs, w_cx, nevt, 438 );
    retrieveWf( wfs, w_cx, nevt, 449 );
    retrieveWf( wfs, w_cx, nevt, 450 );
    retrieveWf( wfs, w_cx, nevt, 485 );
    retrieveWf( wfs, w_cx, nevt, 494 );
#endif
#endif

    // *** DIAGRAM 2741 OF 15495 ***
    // Wavefunction(s) for diagram number 2741
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[435], w_fp[2], COUPs[1], 1.0, 0., 0., w_fp[495] );
    // Amplitude(s) for diagram number 2741
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[290], w_fp[7], w_fp[495], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[401] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[521] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[563] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2742 OF 15495 ***
    // Wavefunction(s) for diagram number 2742
    // (none)
    // Amplitude(s) for diagram number 2742
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[438], w_fp[2], w_fp[290], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] += amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];

    // *** DIAGRAM 2743 OF 15495 ***
    // Wavefunction(s) for diagram number 2743
    // (none)
    // Amplitude(s) for diagram number 2743
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[292], w_fp[6], w_fp[495], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[191] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[521] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[563] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2744 OF 15495 ***
    // Wavefunction(s) for diagram number 2744
    // (none)
    // Amplitude(s) for diagram number 2744
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[437], w_fp[2], w_fp[292], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[191] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];

    // *** DIAGRAM 2745 OF 15495 ***
    // Wavefunction(s) for diagram number 2745
    // (none)
    // Amplitude(s) for diagram number 2745
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[435], w_fp[2], w_fp[299], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[401] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[435], w_fp[2], w_fp[300], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[191] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[521] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[563] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[435], w_fp[2], w_fp[301], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[401] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[521] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[563] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2746 OF 15495 ***
    // Wavefunction(s) for diagram number 2746
    // (none)
    // Amplitude(s) for diagram number 2746
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[450], w_fp[494], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[189] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2747 OF 15495 ***
    // Wavefunction(s) for diagram number 2747
    // (none)
    // Amplitude(s) for diagram number 2747
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[449], w_fp[494], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[179] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2748 OF 15495 ***
    // Wavefunction(s) for diagram number 2748
    // (none)
    // Amplitude(s) for diagram number 2748
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[289], w_fp[7], w_fp[485], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[179] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[299] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[657] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[663] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2749 OF 15495 ***
    // Wavefunction(s) for diagram number 2749
    // (none)
    // Amplitude(s) for diagram number 2749
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[449], w_fp[2], w_fp[289], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[179] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];

    // *** DIAGRAM 2750 OF 15495 ***
    // Wavefunction(s) for diagram number 2750
    // (none)
    // Amplitude(s) for diagram number 2750
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[292], w_fp[4], w_fp[485], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[189] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[299] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[345] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[351] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] += cxtype( 0, 1 ) * amp_sv[0];
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
    storeWf( wfs, w_cx, nevt, 495 );
#endif
#endif
  }

  //--------------------------------------------------------------------------


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup276( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 45 );
    retrieveWf( wfs, w_cx, nevt, 161 );
    retrieveWf( wfs, w_cx, nevt, 279 );
    retrieveWf( wfs, w_cx, nevt, 289 );
    retrieveWf( wfs, w_cx, nevt, 290 );
    retrieveWf( wfs, w_cx, nevt, 292 );
    retrieveWf( wfs, w_cx, nevt, 293 );
    retrieveWf( wfs, w_cx, nevt, 294 );
    retrieveWf( wfs, w_cx, nevt, 295 );
    retrieveWf( wfs, w_cx, nevt, 296 );
    retrieveWf( wfs, w_cx, nevt, 297 );
    retrieveWf( wfs, w_cx, nevt, 298 );
    retrieveWf( wfs, w_cx, nevt, 447 );
    retrieveWf( wfs, w_cx, nevt, 450 );
    retrieveWf( wfs, w_cx, nevt, 451 );
    retrieveWf( wfs, w_cx, nevt, 453 );
    retrieveWf( wfs, w_cx, nevt, 454 );
    retrieveWf( wfs, w_cx, nevt, 486 );
    retrieveWf( wfs, w_cx, nevt, 494 );
#endif
#endif

    // *** DIAGRAM 2751 OF 15495 ***
    // Wavefunction(s) for diagram number 2751
    // (none)
    // Amplitude(s) for diagram number 2751
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[450], w_fp[2], w_fp[292], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[189] += amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    jamp_sv[681] += amp_sv[0];

    // *** DIAGRAM 2752 OF 15495 ***
    // Wavefunction(s) for diagram number 2752
    // (none)
    // Amplitude(s) for diagram number 2752
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[447], w_fp[2], w_fp[296], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[179] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[345] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[351] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[657] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[663] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[447], w_fp[2], w_fp[297], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[189] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[299] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[345] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[351] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[447], w_fp[2], w_fp[298], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[179] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[299] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[657] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[663] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2753 OF 15495 ***
    // Wavefunction(s) for diagram number 2753
    // (none)
    // Amplitude(s) for diagram number 2753
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[454], w_fp[494], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[183] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2754 OF 15495 ***
    // Wavefunction(s) for diagram number 2754
    // (none)
    // Amplitude(s) for diagram number 2754
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[453], w_fp[494], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[177] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2755 OF 15495 ***
    // Wavefunction(s) for diagram number 2755
    // (none)
    // Amplitude(s) for diagram number 2755
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[289], w_fp[6], w_fp[486], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[177] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2756 OF 15495 ***
    // Wavefunction(s) for diagram number 2756
    // (none)
    // Amplitude(s) for diagram number 2756
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[453], w_fp[2], w_fp[289], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[177] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[393] -= amp_sv[0];

    // *** DIAGRAM 2757 OF 15495 ***
    // Wavefunction(s) for diagram number 2757
    // (none)
    // Amplitude(s) for diagram number 2757
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[290], w_fp[4], w_fp[486], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[183] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2758 OF 15495 ***
    // Wavefunction(s) for diagram number 2758
    // (none)
    // Amplitude(s) for diagram number 2758
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[454], w_fp[2], w_fp[290], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[183] += amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[519] -= amp_sv[0];
    jamp_sv[561] += amp_sv[0];

    // *** DIAGRAM 2759 OF 15495 ***
    // Wavefunction(s) for diagram number 2759
    // (none)
    // Amplitude(s) for diagram number 2759
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[2], w_fp[293], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[177] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[183] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[2], w_fp[294], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[183] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[2], w_fp[295], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[177] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2760 OF 15495 ***
    // Wavefunction(s) for diagram number 2760
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[45], w_fp[279], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[496] );
    // Amplitude(s) for diagram number 2760
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[496], w_fp[161], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[333] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] -= cxtype( 0, 1 ) * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 496 );
#endif
#endif
  }

  //--------------------------------------------------------------------------


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup277( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 156 );
    retrieveWf( wfs, w_cx, nevt, 161 );
    retrieveWf( wfs, w_cx, nevt, 163 );
    retrieveWf( wfs, w_cx, nevt, 279 );
    retrieveWf( wfs, w_cx, nevt, 290 );
    retrieveWf( wfs, w_cx, nevt, 292 );
    retrieveWf( wfs, w_cx, nevt, 447 );
    retrieveWf( wfs, w_cx, nevt, 451 );
    retrieveWf( wfs, w_cx, nevt, 496 );
#endif
#endif

    // *** DIAGRAM 2761 OF 15495 ***
    // Wavefunction(s) for diagram number 2761
    // (none)
    // Amplitude(s) for diagram number 2761
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[496], w_fp[163], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[357] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2762 OF 15495 ***
    // Wavefunction(s) for diagram number 2762
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[45], w_fp[156], COUPs[1], 1.0, 0., 0., w_fp[497] );
    // Amplitude(s) for diagram number 2762
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[497], w_fp[290], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[345] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[351] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2763 OF 15495 ***
    // Wavefunction(s) for diagram number 2763
    // (none)
    // Amplitude(s) for diagram number 2763
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[497], w_fp[292], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[275] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[299] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[345] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[351] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2764 OF 15495 ***
    // Wavefunction(s) for diagram number 2764
    // (none)
    // Amplitude(s) for diagram number 2764
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[279], w_fp[6], w_fp[7], w_fp[497], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[345] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[351] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[279], w_fp[6], w_fp[7], w_fp[497], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[275] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[299] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[345] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[351] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[279], w_fp[6], w_fp[7], w_fp[497], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[273] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[299] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2765 OF 15495 ***
    // Wavefunction(s) for diagram number 2765
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[156], w_fp[279], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[498] );
    // Amplitude(s) for diagram number 2765
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[447], w_fp[498], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[275] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[299] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2766 OF 15495 ***
    // Wavefunction(s) for diagram number 2766
    // (none)
    // Amplitude(s) for diagram number 2766
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[447], w_fp[156], w_fp[292], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[275] += amp_sv[0];
    jamp_sv[299] -= amp_sv[0];
    jamp_sv[345] -= amp_sv[0];
    jamp_sv[351] += amp_sv[0];

    // *** DIAGRAM 2767 OF 15495 ***
    // Wavefunction(s) for diagram number 2767
    // (none)
    // Amplitude(s) for diagram number 2767
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[447], w_fp[163], w_fp[279], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[345] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[351] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2768 OF 15495 ***
    // Wavefunction(s) for diagram number 2768
    // (none)
    // Amplitude(s) for diagram number 2768
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[498], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[273] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2769 OF 15495 ***
    // Wavefunction(s) for diagram number 2769
    // (none)
    // Amplitude(s) for diagram number 2769
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[156], w_fp[290], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[273] += amp_sv[0];
    jamp_sv[297] -= amp_sv[0];
    jamp_sv[321] -= amp_sv[0];
    jamp_sv[327] += amp_sv[0];

    // *** DIAGRAM 2770 OF 15495 ***
    // Wavefunction(s) for diagram number 2770
    // (none)
    // Amplitude(s) for diagram number 2770
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[161], w_fp[279], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] -= cxtype( 0, 1 ) * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 497 );
    storeWf( wfs, w_cx, nevt, 498 );
#endif
#endif
  }

  //--------------------------------------------------------------------------


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup278( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 45 );
    retrieveWf( wfs, w_cx, nevt, 84 );
    retrieveWf( wfs, w_cx, nevt, 156 );
    retrieveWf( wfs, w_cx, nevt, 161 );
    retrieveWf( wfs, w_cx, nevt, 163 );
    retrieveWf( wfs, w_cx, nevt, 211 );
    retrieveWf( wfs, w_cx, nevt, 213 );
    retrieveWf( wfs, w_cx, nevt, 279 );
    retrieveWf( wfs, w_cx, nevt, 289 );
    retrieveWf( wfs, w_cx, nevt, 290 );
    retrieveWf( wfs, w_cx, nevt, 292 );
    retrieveWf( wfs, w_cx, nevt, 458 );
    retrieveWf( wfs, w_cx, nevt, 490 );
    retrieveWf( wfs, w_cx, nevt, 496 );
    retrieveWf( wfs, w_cx, nevt, 497 );
#endif
#endif

    // *** DIAGRAM 2771 OF 15495 ***
    // Wavefunction(s) for diagram number 2771
    // (none)
    // Amplitude(s) for diagram number 2771
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[45], w_fp[163], w_fp[290], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[345] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];

    // *** DIAGRAM 2772 OF 15495 ***
    // Wavefunction(s) for diagram number 2772
    // (none)
    // Amplitude(s) for diagram number 2772
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[45], w_fp[161], w_fp[292], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[321] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];

    // *** DIAGRAM 2773 OF 15495 ***
    // Wavefunction(s) for diagram number 2773
    // (none)
    // Amplitude(s) for diagram number 2773
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[496], w_fp[156], w_fp[84], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[333] += amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];

    // *** DIAGRAM 2774 OF 15495 ***
    // Wavefunction(s) for diagram number 2774
    // (none)
    // Amplitude(s) for diagram number 2774
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[279], w_fp[84], w_fp[497], COUPs[0], 1.0, &amp_fp[0] );
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

    // *** DIAGRAM 2775 OF 15495 ***
    // Wavefunction(s) for diagram number 2775
    // (none)
    // Amplitude(s) for diagram number 2775
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[458], w_fp[156], w_fp[279], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[273] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[297] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];

    // *** DIAGRAM 2776 OF 15495 ***
    // Wavefunction(s) for diagram number 2776
    // (none)
    // Amplitude(s) for diagram number 2776
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[496], w_fp[211], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[549] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2777 OF 15495 ***
    // Wavefunction(s) for diagram number 2777
    // (none)
    // Amplitude(s) for diagram number 2777
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[496], w_fp[213], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[591] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2778 OF 15495 ***
    // Wavefunction(s) for diagram number 2778
    // (none)
    // Amplitude(s) for diagram number 2778
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[490], w_fp[289], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2779 OF 15495 ***
    // Wavefunction(s) for diagram number 2779
    // (none)
    // Amplitude(s) for diagram number 2779
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[490], w_fp[292], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[521] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[563] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2780 OF 15495 ***
    // Wavefunction(s) for diagram number 2780
    // (none)
    // Amplitude(s) for diagram number 2780
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[279], w_fp[4], w_fp[7], w_fp[490], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[279], w_fp[4], w_fp[7], w_fp[490], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[521] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[563] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[279], w_fp[4], w_fp[7], w_fp[490], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[519] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[521] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[563] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] -= cxtype( 0, 1 ) * amp_sv[0];

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
  diagramgroup279( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 45 );
    retrieveWf( wfs, w_cx, nevt, 102 );
    retrieveWf( wfs, w_cx, nevt, 197 );
    retrieveWf( wfs, w_cx, nevt, 211 );
    retrieveWf( wfs, w_cx, nevt, 213 );
    retrieveWf( wfs, w_cx, nevt, 279 );
    retrieveWf( wfs, w_cx, nevt, 289 );
    retrieveWf( wfs, w_cx, nevt, 292 );
    retrieveWf( wfs, w_cx, nevt, 435 );
    retrieveWf( wfs, w_cx, nevt, 451 );
    retrieveWf( wfs, w_cx, nevt, 490 );
    retrieveWf( wfs, w_cx, nevt, 496 );
#endif
#endif

    // *** DIAGRAM 2781 OF 15495 ***
    // Wavefunction(s) for diagram number 2781
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[197], w_fp[279], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[499] );
    // Amplitude(s) for diagram number 2781
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[435], w_fp[499], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[521] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[563] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2782 OF 15495 ***
    // Wavefunction(s) for diagram number 2782
    // (none)
    // Amplitude(s) for diagram number 2782
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[435], w_fp[197], w_fp[292], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[521] += amp_sv[0];
    jamp_sv[563] -= amp_sv[0];
    jamp_sv[587] -= amp_sv[0];
    jamp_sv[597] += amp_sv[0];

    // *** DIAGRAM 2783 OF 15495 ***
    // Wavefunction(s) for diagram number 2783
    // (none)
    // Amplitude(s) for diagram number 2783
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[435], w_fp[213], w_fp[279], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[587] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2784 OF 15495 ***
    // Wavefunction(s) for diagram number 2784
    // (none)
    // Amplitude(s) for diagram number 2784
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[499], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[519] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2785 OF 15495 ***
    // Wavefunction(s) for diagram number 2785
    // (none)
    // Amplitude(s) for diagram number 2785
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[197], w_fp[289], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[519] += amp_sv[0];
    jamp_sv[537] -= amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[561] -= amp_sv[0];

    // *** DIAGRAM 2786 OF 15495 ***
    // Wavefunction(s) for diagram number 2786
    // (none)
    // Amplitude(s) for diagram number 2786
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[211], w_fp[279], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[537] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2787 OF 15495 ***
    // Wavefunction(s) for diagram number 2787
    // (none)
    // Amplitude(s) for diagram number 2787
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[45], w_fp[213], w_fp[289], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[587] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[593] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];

    // *** DIAGRAM 2788 OF 15495 ***
    // Wavefunction(s) for diagram number 2788
    // (none)
    // Amplitude(s) for diagram number 2788
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[45], w_fp[211], w_fp[292], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[537] += amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[549] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];

    // *** DIAGRAM 2789 OF 15495 ***
    // Wavefunction(s) for diagram number 2789
    // (none)
    // Amplitude(s) for diagram number 2789
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[496], w_fp[197], w_fp[102], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[549] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[593] += amp_sv[0];

    // *** DIAGRAM 2790 OF 15495 ***
    // Wavefunction(s) for diagram number 2790
    // (none)
    // Amplitude(s) for diagram number 2790
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[279], w_fp[102], w_fp[490], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[521] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[563] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] += cxtype( 0, 1 ) * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 499 );
#endif
#endif
  }

  //--------------------------------------------------------------------------


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup280( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 197 );
    retrieveWf( wfs, w_cx, nevt, 215 );
    retrieveWf( wfs, w_cx, nevt, 225 );
    retrieveWf( wfs, w_cx, nevt, 227 );
    retrieveWf( wfs, w_cx, nevt, 279 );
    retrieveWf( wfs, w_cx, nevt, 289 );
    retrieveWf( wfs, w_cx, nevt, 290 );
    retrieveWf( wfs, w_cx, nevt, 435 );
    retrieveWf( wfs, w_cx, nevt, 447 );
    retrieveWf( wfs, w_cx, nevt, 462 );
    retrieveWf( wfs, w_cx, nevt, 492 );
    retrieveWf( wfs, w_cx, nevt, 496 );
#endif
#endif

    // *** DIAGRAM 2791 OF 15495 ***
    // Wavefunction(s) for diagram number 2791
    // (none)
    // Amplitude(s) for diagram number 2791
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[462], w_fp[197], w_fp[279], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[519] += amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[561] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];

    // *** DIAGRAM 2792 OF 15495 ***
    // Wavefunction(s) for diagram number 2792
    // (none)
    // Amplitude(s) for diagram number 2792
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[496], w_fp[225], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[669] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[671] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2793 OF 15495 ***
    // Wavefunction(s) for diagram number 2793
    // (none)
    // Amplitude(s) for diagram number 2793
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[496], w_fp[227], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[711] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[713] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2794 OF 15495 ***
    // Wavefunction(s) for diagram number 2794
    // (none)
    // Amplitude(s) for diagram number 2794
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[492], w_fp[289], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[639] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[657] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[663] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[713] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2795 OF 15495 ***
    // Wavefunction(s) for diagram number 2795
    // (none)
    // Amplitude(s) for diagram number 2795
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[492], w_fp[290], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[641] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[657] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[663] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[671] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2796 OF 15495 ***
    // Wavefunction(s) for diagram number 2796
    // (none)
    // Amplitude(s) for diagram number 2796
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[279], w_fp[4], w_fp[6], w_fp[492], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[639] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[657] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[663] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[713] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[279], w_fp[4], w_fp[6], w_fp[492], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[641] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[657] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[663] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[671] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[279], w_fp[4], w_fp[6], w_fp[492], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[639] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[671] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[713] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2797 OF 15495 ***
    // Wavefunction(s) for diagram number 2797
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[215], w_fp[279], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[500] );
    // Amplitude(s) for diagram number 2797
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[435], w_fp[500], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[641] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2798 OF 15495 ***
    // Wavefunction(s) for diagram number 2798
    // (none)
    // Amplitude(s) for diagram number 2798
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[435], w_fp[215], w_fp[290], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[641] += amp_sv[0];
    jamp_sv[683] -= amp_sv[0];
    jamp_sv[707] -= amp_sv[0];
    jamp_sv[717] += amp_sv[0];

    // *** DIAGRAM 2799 OF 15495 ***
    // Wavefunction(s) for diagram number 2799
    // (none)
    // Amplitude(s) for diagram number 2799
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[435], w_fp[227], w_fp[279], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[707] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2800 OF 15495 ***
    // Wavefunction(s) for diagram number 2800
    // (none)
    // Amplitude(s) for diagram number 2800
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[447], w_fp[500], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[639] += cxtype( 0, 1 ) * amp_sv[0];
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
    storeWf( wfs, w_cx, nevt, 500 );
#endif
#endif
  }

  //--------------------------------------------------------------------------

}
