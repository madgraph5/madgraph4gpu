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
  diagramgroup951( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 39 );
    retrieveWf( wfs, w_cx, nevt, 65 );
    retrieveWf( wfs, w_cx, nevt, 66 );
    retrieveWf( wfs, w_cx, nevt, 68 );
    retrieveWf( wfs, w_cx, nevt, 77 );
    retrieveWf( wfs, w_cx, nevt, 84 );
    retrieveWf( wfs, w_cx, nevt, 128 );
    retrieveWf( wfs, w_cx, nevt, 148 );
    retrieveWf( wfs, w_cx, nevt, 454 );
    retrieveWf( wfs, w_cx, nevt, 514 );
    retrieveWf( wfs, w_cx, nevt, 519 );
    retrieveWf( wfs, w_cx, nevt, 521 );
    retrieveWf( wfs, w_cx, nevt, 524 );
    retrieveWf( wfs, w_cx, nevt, 526 );
    retrieveWf( wfs, w_cx, nevt, 534 );
    retrieveWf( wfs, w_cx, nevt, 575 );
    retrieveWf( wfs, w_cx, nevt, 601 );
    retrieveWf( wfs, w_cx, nevt, 662 );
    retrieveWf( wfs, w_cx, nevt, 675 );
#endif
#endif

    // *** DIAGRAM 9501 OF 15495 ***
    // Wavefunction(s) for diagram number 9501
    // (none)
    // Amplitude(s) for diagram number 9501
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[675], w_fp[2], w_fp[68], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[301] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[565] += amp_sv[0];

    // *** DIAGRAM 9502 OF 15495 ***
    // Wavefunction(s) for diagram number 9502
    // (none)
    // Amplitude(s) for diagram number 9502
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[77], w_fp[2], w_fp[601], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9503 OF 15495 ***
    // Wavefunction(s) for diagram number 9503
    // (none)
    // Amplitude(s) for diagram number 9503
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[519], w_fp[534], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] += amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[87] += amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[565] -= amp_sv[0];
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[454], w_fp[534], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[87] += amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[489] += amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[685] += amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[699] -= amp_sv[0];
    jamp_sv[701] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[575], w_fp[534], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[489] += amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    jamp_sv[609] += amp_sv[0];
    jamp_sv[615] -= amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[685] += amp_sv[0];
    jamp_sv[699] -= amp_sv[0];
    jamp_sv[701] += amp_sv[0];
    jamp_sv[712] += amp_sv[0];
    jamp_sv[718] -= amp_sv[0];

    // *** DIAGRAM 9504 OF 15495 ***
    // Wavefunction(s) for diagram number 9504
    // (none)
    // Amplitude(s) for diagram number 9504
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[77], w_fp[2], w_fp[519], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[77], w_fp[2], w_fp[454], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[81] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[77], w_fp[2], w_fp[575], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9505 OF 15495 ***
    // Wavefunction(s) for diagram number 9505
    // (none)
    // Amplitude(s) for diagram number 9505
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[521], w_fp[534], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[35] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[111] += amp_sv[0];
    jamp_sv[307] -= amp_sv[0];
    jamp_sv[427] += amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[501] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[598] += amp_sv[0];
    jamp_sv[661] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[526], w_fp[534], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[111] += amp_sv[0];
    jamp_sv[293] += amp_sv[0];
    jamp_sv[307] -= amp_sv[0];
    jamp_sv[413] -= amp_sv[0];
    jamp_sv[427] += amp_sv[0];
    jamp_sv[501] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[579] -= amp_sv[0];
    jamp_sv[581] += amp_sv[0];
    jamp_sv[609] += amp_sv[0];
    jamp_sv[615] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[524], w_fp[534], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[293] += amp_sv[0];
    jamp_sv[413] -= amp_sv[0];
    jamp_sv[489] += amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    jamp_sv[579] -= amp_sv[0];
    jamp_sv[581] += amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[609] += amp_sv[0];
    jamp_sv[615] -= amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[685] += amp_sv[0];

    // *** DIAGRAM 9506 OF 15495 ***
    // Wavefunction(s) for diagram number 9506
    // (none)
    // Amplitude(s) for diagram number 9506
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[39], w_fp[2], w_fp[521], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[35] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[307] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[427] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[39], w_fp[2], w_fp[526], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[105] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[307] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[427] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[39], w_fp[2], w_fp[524], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[35] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9507 OF 15495 ***
    // Wavefunction(s) for diagram number 9507
    // (none)
    // Amplitude(s) for diagram number 9507
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[148], w_fp[84], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[304] += amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[430] += amp_sv[0];

    // *** DIAGRAM 9508 OF 15495 ***
    // Wavefunction(s) for diagram number 9508
    // (none)
    // Amplitude(s) for diagram number 9508
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[128], w_fp[66], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[592] += amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];

    // *** DIAGRAM 9509 OF 15495 ***
    // Wavefunction(s) for diagram number 9509
    // (none)
    // Amplitude(s) for diagram number 9509
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[2], w_fp[65], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[304] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[310] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[424] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9510 OF 15495 ***
    // Wavefunction(s) for diagram number 9510
    // (none)
    // Amplitude(s) for diagram number 9510
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[514], w_fp[534], w_fp[84], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] += amp_sv[0];
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[293] += amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[413] -= amp_sv[0];
    jamp_sv[579] -= amp_sv[0];
    jamp_sv[581] += amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    jamp_sv[701] -= amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];

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
  diagramgroup952( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 65 );
    retrieveWf( wfs, w_cx, nevt, 66 );
    retrieveWf( wfs, w_cx, nevt, 92 );
    retrieveWf( wfs, w_cx, nevt, 118 );
    retrieveWf( wfs, w_cx, nevt, 128 );
    retrieveWf( wfs, w_cx, nevt, 148 );
    retrieveWf( wfs, w_cx, nevt, 157 );
    retrieveWf( wfs, w_cx, nevt, 246 );
    retrieveWf( wfs, w_cx, nevt, 457 );
    retrieveWf( wfs, w_cx, nevt, 467 );
    retrieveWf( wfs, w_cx, nevt, 514 );
    retrieveWf( wfs, w_cx, nevt, 516 );
    retrieveWf( wfs, w_cx, nevt, 534 );
    retrieveWf( wfs, w_cx, nevt, 604 );
    retrieveWf( wfs, w_cx, nevt, 665 );
#endif
#endif

    // *** DIAGRAM 9511 OF 15495 ***
    // Wavefunction(s) for diagram number 9511
    // (none)
    // Amplitude(s) for diagram number 9511
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[246], w_fp[2], w_fp[514], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9512 OF 15495 ***
    // Wavefunction(s) for diagram number 9512
    // (none)
    // Amplitude(s) for diagram number 9512
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[157], w_fp[128], w_fp[514], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[579] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9513 OF 15495 ***
    // Wavefunction(s) for diagram number 9513
    // (none)
    // Amplitude(s) for diagram number 9513
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[516], w_fp[534], w_fp[66], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[93] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[293] += amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[413] -= amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[579] -= amp_sv[0];
    jamp_sv[581] += amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    jamp_sv[701] -= amp_sv[0];

    // *** DIAGRAM 9514 OF 15495 ***
    // Wavefunction(s) for diagram number 9514
    // (none)
    // Amplitude(s) for diagram number 9514
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[92], w_fp[2], w_fp[516], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[93] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[95] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9515 OF 15495 ***
    // Wavefunction(s) for diagram number 9515
    // (none)
    // Amplitude(s) for diagram number 9515
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[157], w_fp[148], w_fp[516], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[304] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[310] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[424] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9516 OF 15495 ***
    // Wavefunction(s) for diagram number 9516
    // (none)
    // Amplitude(s) for diagram number 9516
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[534], w_fp[65], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] += amp_sv[0];
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];

    // *** DIAGRAM 9517 OF 15495 ***
    // Wavefunction(s) for diagram number 9517
    // (none)
    // Amplitude(s) for diagram number 9517
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[92], w_fp[128], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[579] += amp_sv[0];
    jamp_sv[581] -= amp_sv[0];
    jamp_sv[699] -= amp_sv[0];
    jamp_sv[701] += amp_sv[0];

    // *** DIAGRAM 9518 OF 15495 ***
    // Wavefunction(s) for diagram number 9518
    // (none)
    // Amplitude(s) for diagram number 9518
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[246], w_fp[148], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[291] += amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[413] += amp_sv[0];

    // *** DIAGRAM 9519 OF 15495 ***
    // Wavefunction(s) for diagram number 9519
    // (none)
    // Amplitude(s) for diagram number 9519
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[157], w_fp[2], w_fp[467], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[35] += amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[598] += amp_sv[0];
    jamp_sv[712] += amp_sv[0];
    jamp_sv[718] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[157], w_fp[2], w_fp[457], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[93] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[293] += amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[413] -= amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[579] -= amp_sv[0];
    jamp_sv[581] += amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    jamp_sv[701] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[157], w_fp[2], w_fp[604], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] += amp_sv[0];
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[293] += amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[413] -= amp_sv[0];
    jamp_sv[579] -= amp_sv[0];
    jamp_sv[581] += amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    jamp_sv[701] -= amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];

    // *** DIAGRAM 9520 OF 15495 ***
    // Wavefunction(s) for diagram number 9520
    // (none)
    // Amplitude(s) for diagram number 9520
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[665], w_fp[118], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];

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
  diagramgroup953( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 29 );
    retrieveWf( wfs, w_cx, nevt, 64 );
    retrieveWf( wfs, w_cx, nevt, 70 );
    retrieveWf( wfs, w_cx, nevt, 85 );
    retrieveWf( wfs, w_cx, nevt, 87 );
    retrieveWf( wfs, w_cx, nevt, 88 );
    retrieveWf( wfs, w_cx, nevt, 118 );
    retrieveWf( wfs, w_cx, nevt, 534 );
    retrieveWf( wfs, w_cx, nevt, 544 );
    retrieveWf( wfs, w_cx, nevt, 576 );
    retrieveWf( wfs, w_cx, nevt, 605 );
    retrieveWf( wfs, w_cx, nevt, 606 );
    retrieveWf( wfs, w_cx, nevt, 662 );
    retrieveWf( wfs, w_cx, nevt, 664 );
    retrieveWf( wfs, w_cx, nevt, 665 );
#endif
#endif

    // *** DIAGRAM 9521 OF 15495 ***
    // Wavefunction(s) for diagram number 9521
    // (none)
    // Amplitude(s) for diagram number 9521
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[664], w_fp[118], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[328] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[544] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9522 OF 15495 ***
    // Wavefunction(s) for diagram number 9522
    // (none)
    // Amplitude(s) for diagram number 9522
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[87], w_fp[7], w_fp[544], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[328] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[424] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[544] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[688] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[694] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9523 OF 15495 ***
    // Wavefunction(s) for diagram number 9523
    // (none)
    // Amplitude(s) for diagram number 9523
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[664], w_fp[2], w_fp[87], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[328] += amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[544] -= amp_sv[0];

    // *** DIAGRAM 9524 OF 15495 ***
    // Wavefunction(s) for diagram number 9524
    // (none)
    // Amplitude(s) for diagram number 9524
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[88], w_fp[5], w_fp[544], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[424] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9525 OF 15495 ***
    // Wavefunction(s) for diagram number 9525
    // (none)
    // Amplitude(s) for diagram number 9525
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[665], w_fp[2], w_fp[88], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[334] += amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[712] += amp_sv[0];

    // *** DIAGRAM 9526 OF 15495 ***
    // Wavefunction(s) for diagram number 9526
    // (none)
    // Amplitude(s) for diagram number 9526
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[2], w_fp[85], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[328] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[544] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[688] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[694] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[2], w_fp[70], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[424] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[2], w_fp[64], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[328] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[424] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[544] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[688] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[694] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9527 OF 15495 ***
    // Wavefunction(s) for diagram number 9527
    // (none)
    // Amplitude(s) for diagram number 9527
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[576], w_fp[534], w_fp[5], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[39] += amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[459] -= amp_sv[0];
    jamp_sv[461] += amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[533] -= amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[576], w_fp[534], w_fp[5], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[39] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[611] -= amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[576], w_fp[534], w_fp[5], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[41] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[459] += amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    jamp_sv[533] += amp_sv[0];
    jamp_sv[611] -= amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];

    // *** DIAGRAM 9528 OF 15495 ***
    // Wavefunction(s) for diagram number 9528
    // (none)
    // Amplitude(s) for diagram number 9528
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[534], w_fp[7], w_fp[605], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[39] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[611] -= amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];

    // *** DIAGRAM 9529 OF 15495 ***
    // Wavefunction(s) for diagram number 9529
    // (none)
    // Amplitude(s) for diagram number 9529
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[534], w_fp[5], w_fp[606], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[41] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[459] += amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    jamp_sv[533] += amp_sv[0];
    jamp_sv[611] -= amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];

    // *** DIAGRAM 9530 OF 15495 ***
    // Wavefunction(s) for diagram number 9530
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[576], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[246] );
    // Amplitude(s) for diagram number 9530
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[29], w_fp[246], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[41] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[533] += amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 246 );
#endif
#endif
  }

  //--------------------------------------------------------------------------


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup954( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 29 );
    retrieveWf( wfs, w_cx, nevt, 64 );
    retrieveWf( wfs, w_cx, nevt, 70 );
    retrieveWf( wfs, w_cx, nevt, 77 );
    retrieveWf( wfs, w_cx, nevt, 85 );
    retrieveWf( wfs, w_cx, nevt, 87 );
    retrieveWf( wfs, w_cx, nevt, 88 );
    retrieveWf( wfs, w_cx, nevt, 246 );
    retrieveWf( wfs, w_cx, nevt, 534 );
    retrieveWf( wfs, w_cx, nevt, 605 );
    retrieveWf( wfs, w_cx, nevt, 606 );
    retrieveWf( wfs, w_cx, nevt, 608 );
    retrieveWf( wfs, w_cx, nevt, 609 );
    retrieveWf( wfs, w_cx, nevt, 697 );
#endif
#endif

    // *** DIAGRAM 9531 OF 15495 ***
    // Wavefunction(s) for diagram number 9531
    // (none)
    // Amplitude(s) for diagram number 9531
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[29], w_fp[2], w_fp[606], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[41] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[533] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[611] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9532 OF 15495 ***
    // Wavefunction(s) for diagram number 9532
    // (none)
    // Amplitude(s) for diagram number 9532
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[77], w_fp[246], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[39] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[531] += amp_sv[0];

    // *** DIAGRAM 9533 OF 15495 ***
    // Wavefunction(s) for diagram number 9533
    // (none)
    // Amplitude(s) for diagram number 9533
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[77], w_fp[2], w_fp[605], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[39] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[445] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9534 OF 15495 ***
    // Wavefunction(s) for diagram number 9534
    // (none)
    // Amplitude(s) for diagram number 9534
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[534], w_fp[87], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[39] += amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[111] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[534], w_fp[87], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[39] += amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[611] -= amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[534], w_fp[87], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[107] += amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[611] -= amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[621] += amp_sv[0];

    // *** DIAGRAM 9535 OF 15495 ***
    // Wavefunction(s) for diagram number 9535
    // (none)
    // Amplitude(s) for diagram number 9535
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[87], w_fp[7], w_fp[697], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[39] += amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[111] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[712] -= amp_sv[0];

    // *** DIAGRAM 9536 OF 15495 ***
    // Wavefunction(s) for diagram number 9536
    // (none)
    // Amplitude(s) for diagram number 9536
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[534], w_fp[7], w_fp[608], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[39] += amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[611] -= amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[712] -= amp_sv[0];

    // *** DIAGRAM 9537 OF 15495 ***
    // Wavefunction(s) for diagram number 9537
    // (none)
    // Amplitude(s) for diagram number 9537
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[534], w_fp[88], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[41] += amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[69] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[534], w_fp[88], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[41] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    jamp_sv[331] -= amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    jamp_sv[547] += amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[534], w_fp[88], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    jamp_sv[331] -= amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[547] += amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    jamp_sv[712] += amp_sv[0];

    // *** DIAGRAM 9538 OF 15495 ***
    // Wavefunction(s) for diagram number 9538
    // (none)
    // Amplitude(s) for diagram number 9538
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[88], w_fp[5], w_fp[697], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[41] += amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[69] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    jamp_sv[712] -= amp_sv[0];

    // *** DIAGRAM 9539 OF 15495 ***
    // Wavefunction(s) for diagram number 9539
    // (none)
    // Amplitude(s) for diagram number 9539
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[534], w_fp[5], w_fp[609], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[41] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    jamp_sv[331] -= amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    jamp_sv[547] += amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];

    // *** DIAGRAM 9540 OF 15495 ***
    // Wavefunction(s) for diagram number 9540
    // (none)
    // Amplitude(s) for diagram number 9540
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[534], w_fp[85], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[39] += amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[111] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[534], w_fp[70], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[712] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[534], w_fp[64], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[712] += amp_sv[0];

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
  diagramgroup955( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 29 );
    retrieveWf( wfs, w_cx, nevt, 77 );
    retrieveWf( wfs, w_cx, nevt, 87 );
    retrieveWf( wfs, w_cx, nevt, 88 );
    retrieveWf( wfs, w_cx, nevt, 118 );
    retrieveWf( wfs, w_cx, nevt, 447 );
    retrieveWf( wfs, w_cx, nevt, 534 );
    retrieveWf( wfs, w_cx, nevt, 561 );
    retrieveWf( wfs, w_cx, nevt, 590 );
    retrieveWf( wfs, w_cx, nevt, 608 );
    retrieveWf( wfs, w_cx, nevt, 609 );
    retrieveWf( wfs, w_cx, nevt, 670 );
    retrieveWf( wfs, w_cx, nevt, 675 );
#endif
#endif

    // *** DIAGRAM 9541 OF 15495 ***
    // Wavefunction(s) for diagram number 9541
    // (none)
    // Amplitude(s) for diagram number 9541
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[670], w_fp[118], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9542 OF 15495 ***
    // Wavefunction(s) for diagram number 9542
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[118], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[92] );
    // Amplitude(s) for diagram number 9542
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[29], w_fp[92], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[317] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[533] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9543 OF 15495 ***
    // Wavefunction(s) for diagram number 9543
    // (none)
    // Amplitude(s) for diagram number 9543
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[670], w_fp[2], w_fp[88], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[331] += amp_sv[0];
    jamp_sv[547] -= amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[709] += amp_sv[0];

    // *** DIAGRAM 9544 OF 15495 ***
    // Wavefunction(s) for diagram number 9544
    // (none)
    // Amplitude(s) for diagram number 9544
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[29], w_fp[2], w_fp[609], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[41] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9545 OF 15495 ***
    // Wavefunction(s) for diagram number 9545
    // (none)
    // Amplitude(s) for diagram number 9545
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[675], w_fp[118], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9546 OF 15495 ***
    // Wavefunction(s) for diagram number 9546
    // (none)
    // Amplitude(s) for diagram number 9546
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[77], w_fp[92], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9547 OF 15495 ***
    // Wavefunction(s) for diagram number 9547
    // (none)
    // Amplitude(s) for diagram number 9547
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[675], w_fp[2], w_fp[87], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[325] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[445] += amp_sv[0];
    jamp_sv[541] -= amp_sv[0];

    // *** DIAGRAM 9548 OF 15495 ***
    // Wavefunction(s) for diagram number 9548
    // (none)
    // Amplitude(s) for diagram number 9548
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[77], w_fp[2], w_fp[608], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[39] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[445] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9549 OF 15495 ***
    // Wavefunction(s) for diagram number 9549
    // (none)
    // Amplitude(s) for diagram number 9549
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[561], w_fp[534], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[39] += amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[611] -= amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[447], w_fp[534], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[369] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[531] -= amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[709] += amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[590], w_fp[534], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[369] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[445] += amp_sv[0];
    jamp_sv[531] -= amp_sv[0];
    jamp_sv[611] += amp_sv[0];
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[709] += amp_sv[0];

    // *** DIAGRAM 9550 OF 15495 ***
    // Wavefunction(s) for diagram number 9550
    // (none)
    // Amplitude(s) for diagram number 9550
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[77], w_fp[2], w_fp[561], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[39] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[445] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[77], w_fp[2], w_fp[447], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[77], w_fp[2], w_fp[590], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[39] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[445] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 92 );
#endif
#endif
  }

  //--------------------------------------------------------------------------


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup956( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 29 );
    retrieveWf( wfs, w_cx, nevt, 62 );
    retrieveWf( wfs, w_cx, nevt, 86 );
    retrieveWf( wfs, w_cx, nevt, 100 );
    retrieveWf( wfs, w_cx, nevt, 110 );
    retrieveWf( wfs, w_cx, nevt, 118 );
    retrieveWf( wfs, w_cx, nevt, 122 );
    retrieveWf( wfs, w_cx, nevt, 157 );
    retrieveWf( wfs, w_cx, nevt, 245 );
    retrieveWf( wfs, w_cx, nevt, 480 );
    retrieveWf( wfs, w_cx, nevt, 522 );
    retrieveWf( wfs, w_cx, nevt, 534 );
    retrieveWf( wfs, w_cx, nevt, 550 );
    retrieveWf( wfs, w_cx, nevt, 558 );
    retrieveWf( wfs, w_cx, nevt, 576 );
    retrieveWf( wfs, w_cx, nevt, 662 );
#endif
#endif

    // *** DIAGRAM 9551 OF 15495 ***
    // Wavefunction(s) for diagram number 9551
    // (none)
    // Amplitude(s) for diagram number 9551
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[550], w_fp[534], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[41] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    jamp_sv[331] -= amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    jamp_sv[547] += amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[480], w_fp[534], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[331] -= amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[445] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[459] -= amp_sv[0];
    jamp_sv[461] += amp_sv[0];
    jamp_sv[533] -= amp_sv[0];
    jamp_sv[547] += amp_sv[0];
    jamp_sv[611] += amp_sv[0];
    jamp_sv[621] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[558], w_fp[534], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[369] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[445] += amp_sv[0];
    jamp_sv[459] -= amp_sv[0];
    jamp_sv[461] += amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[533] -= amp_sv[0];
    jamp_sv[611] += amp_sv[0];
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[709] += amp_sv[0];

    // *** DIAGRAM 9552 OF 15495 ***
    // Wavefunction(s) for diagram number 9552
    // (none)
    // Amplitude(s) for diagram number 9552
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[29], w_fp[2], w_fp[550], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[41] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[29], w_fp[2], w_fp[480], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[107] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[533] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[611] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[29], w_fp[2], w_fp[558], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[41] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[533] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[611] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9553 OF 15495 ***
    // Wavefunction(s) for diagram number 9553
    // (none)
    // Amplitude(s) for diagram number 9553
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[118], w_fp[100], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[328] += amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[550] += amp_sv[0];

    // *** DIAGRAM 9554 OF 15495 ***
    // Wavefunction(s) for diagram number 9554
    // (none)
    // Amplitude(s) for diagram number 9554
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[122], w_fp[86], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[472] += amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];

    // *** DIAGRAM 9555 OF 15495 ***
    // Wavefunction(s) for diagram number 9555
    // (none)
    // Amplitude(s) for diagram number 9555
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[2], w_fp[62], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[328] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[544] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[688] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[694] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9556 OF 15495 ***
    // Wavefunction(s) for diagram number 9556
    // (none)
    // Amplitude(s) for diagram number 9556
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[576], w_fp[534], w_fp[100], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[39] += amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[459] -= amp_sv[0];
    jamp_sv[461] += amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[533] -= amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];

    // *** DIAGRAM 9557 OF 15495 ***
    // Wavefunction(s) for diagram number 9557
    // (none)
    // Amplitude(s) for diagram number 9557
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[245], w_fp[2], w_fp[576], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[39] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[41] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[533] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9558 OF 15495 ***
    // Wavefunction(s) for diagram number 9558
    // (none)
    // Amplitude(s) for diagram number 9558
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[157], w_fp[122], w_fp[576], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[459] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[688] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[694] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9559 OF 15495 ***
    // Wavefunction(s) for diagram number 9559
    // (none)
    // Amplitude(s) for diagram number 9559
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[522], w_fp[534], w_fp[86], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[69] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[459] -= amp_sv[0];
    jamp_sv[461] += amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[533] -= amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];

    // *** DIAGRAM 9560 OF 15495 ***
    // Wavefunction(s) for diagram number 9560
    // (none)
    // Amplitude(s) for diagram number 9560
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[110], w_fp[2], w_fp[522], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[459] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];

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
  diagramgroup957( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 62 );
    retrieveWf( wfs, w_cx, nevt, 98 );
    retrieveWf( wfs, w_cx, nevt, 103 );
    retrieveWf( wfs, w_cx, nevt, 104 );
    retrieveWf( wfs, w_cx, nevt, 110 );
    retrieveWf( wfs, w_cx, nevt, 118 );
    retrieveWf( wfs, w_cx, nevt, 122 );
    retrieveWf( wfs, w_cx, nevt, 157 );
    retrieveWf( wfs, w_cx, nevt, 245 );
    retrieveWf( wfs, w_cx, nevt, 460 );
    retrieveWf( wfs, w_cx, nevt, 466 );
    retrieveWf( wfs, w_cx, nevt, 522 );
    retrieveWf( wfs, w_cx, nevt, 534 );
    retrieveWf( wfs, w_cx, nevt, 544 );
    retrieveWf( wfs, w_cx, nevt, 611 );
    retrieveWf( wfs, w_cx, nevt, 663 );
    retrieveWf( wfs, w_cx, nevt, 665 );
#endif
#endif

    // *** DIAGRAM 9561 OF 15495 ***
    // Wavefunction(s) for diagram number 9561
    // (none)
    // Amplitude(s) for diagram number 9561
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[157], w_fp[118], w_fp[522], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[328] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[533] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[544] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9562 OF 15495 ***
    // Wavefunction(s) for diagram number 9562
    // (none)
    // Amplitude(s) for diagram number 9562
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[534], w_fp[62], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[39] += amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[111] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];

    // *** DIAGRAM 9563 OF 15495 ***
    // Wavefunction(s) for diagram number 9563
    // (none)
    // Amplitude(s) for diagram number 9563
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[110], w_fp[122], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[459] += amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];

    // *** DIAGRAM 9564 OF 15495 ***
    // Wavefunction(s) for diagram number 9564
    // (none)
    // Amplitude(s) for diagram number 9564
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[245], w_fp[118], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[315] += amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[531] -= amp_sv[0];
    jamp_sv[533] += amp_sv[0];

    // *** DIAGRAM 9565 OF 15495 ***
    // Wavefunction(s) for diagram number 9565
    // (none)
    // Amplitude(s) for diagram number 9565
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[157], w_fp[2], w_fp[466], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[69] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[157], w_fp[2], w_fp[460], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[69] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[459] -= amp_sv[0];
    jamp_sv[461] += amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[533] -= amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[157], w_fp[2], w_fp[611], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[39] += amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[459] -= amp_sv[0];
    jamp_sv[461] += amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[533] -= amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];

    // *** DIAGRAM 9566 OF 15495 ***
    // Wavefunction(s) for diagram number 9566
    // (none)
    // Amplitude(s) for diagram number 9566
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[665], w_fp[98], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9567 OF 15495 ***
    // Wavefunction(s) for diagram number 9567
    // (none)
    // Amplitude(s) for diagram number 9567
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[663], w_fp[98], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[352] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[664] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9568 OF 15495 ***
    // Wavefunction(s) for diagram number 9568
    // (none)
    // Amplitude(s) for diagram number 9568
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[103], w_fp[6], w_fp[544], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[352] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[574] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[664] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9569 OF 15495 ***
    // Wavefunction(s) for diagram number 9569
    // (none)
    // Amplitude(s) for diagram number 9569
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[663], w_fp[2], w_fp[103], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[352] += amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[664] -= amp_sv[0];

    // *** DIAGRAM 9570 OF 15495 ***
    // Wavefunction(s) for diagram number 9570
    // (none)
    // Amplitude(s) for diagram number 9570
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[104], w_fp[5], w_fp[544], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[454] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];

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
  diagramgroup958( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 29 );
    retrieveWf( wfs, w_cx, nevt, 39 );
    retrieveWf( wfs, w_cx, nevt, 61 );
    retrieveWf( wfs, w_cx, nevt, 89 );
    retrieveWf( wfs, w_cx, nevt, 101 );
    retrieveWf( wfs, w_cx, nevt, 103 );
    retrieveWf( wfs, w_cx, nevt, 104 );
    retrieveWf( wfs, w_cx, nevt, 518 );
    retrieveWf( wfs, w_cx, nevt, 534 );
    retrieveWf( wfs, w_cx, nevt, 612 );
    retrieveWf( wfs, w_cx, nevt, 613 );
    retrieveWf( wfs, w_cx, nevt, 662 );
    retrieveWf( wfs, w_cx, nevt, 665 );
#endif
#endif

    // *** DIAGRAM 9571 OF 15495 ***
    // Wavefunction(s) for diagram number 9571
    // (none)
    // Amplitude(s) for diagram number 9571
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[665], w_fp[2], w_fp[104], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[358] += amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[670] -= amp_sv[0];

    // *** DIAGRAM 9572 OF 15495 ***
    // Wavefunction(s) for diagram number 9572
    // (none)
    // Amplitude(s) for diagram number 9572
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[2], w_fp[101], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[352] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[454] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[574] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[664] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[2], w_fp[89], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[454] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[2], w_fp[61], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[352] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[574] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[664] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9573 OF 15495 ***
    // Wavefunction(s) for diagram number 9573
    // (none)
    // Amplitude(s) for diagram number 9573
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[518], w_fp[534], w_fp[5], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[45] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[339] -= amp_sv[0];
    jamp_sv[341] += amp_sv[0];
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[555] += amp_sv[0];
    jamp_sv[557] -= amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[518], w_fp[534], w_fp[5], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[45] += amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[339] -= amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[427] += amp_sv[0];
    jamp_sv[469] -= amp_sv[0];
    jamp_sv[491] -= amp_sv[0];
    jamp_sv[501] += amp_sv[0];
    jamp_sv[547] += amp_sv[0];
    jamp_sv[555] += amp_sv[0];
    jamp_sv[557] -= amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[589] -= amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[518], w_fp[534], w_fp[5], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[47] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[341] -= amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[427] += amp_sv[0];
    jamp_sv[435] += amp_sv[0];
    jamp_sv[437] -= amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[469] -= amp_sv[0];
    jamp_sv[491] -= amp_sv[0];
    jamp_sv[501] += amp_sv[0];
    jamp_sv[547] += amp_sv[0];
    jamp_sv[589] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];

    // *** DIAGRAM 9574 OF 15495 ***
    // Wavefunction(s) for diagram number 9574
    // (none)
    // Amplitude(s) for diagram number 9574
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[534], w_fp[6], w_fp[612], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[45] += amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[339] -= amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[427] += amp_sv[0];
    jamp_sv[469] -= amp_sv[0];
    jamp_sv[491] -= amp_sv[0];
    jamp_sv[501] += amp_sv[0];
    jamp_sv[547] += amp_sv[0];
    jamp_sv[555] += amp_sv[0];
    jamp_sv[557] -= amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[589] -= amp_sv[0];
    jamp_sv[651] += amp_sv[0];

    // *** DIAGRAM 9575 OF 15495 ***
    // Wavefunction(s) for diagram number 9575
    // (none)
    // Amplitude(s) for diagram number 9575
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[534], w_fp[5], w_fp[613], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[47] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[341] -= amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[427] += amp_sv[0];
    jamp_sv[435] += amp_sv[0];
    jamp_sv[437] -= amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[469] -= amp_sv[0];
    jamp_sv[491] -= amp_sv[0];
    jamp_sv[501] += amp_sv[0];
    jamp_sv[547] += amp_sv[0];
    jamp_sv[589] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];

    // *** DIAGRAM 9576 OF 15495 ***
    // Wavefunction(s) for diagram number 9576
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[518], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[245] );
    // Amplitude(s) for diagram number 9576
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[29], w_fp[245], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[47] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[341] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];

    // *** DIAGRAM 9577 OF 15495 ***
    // Wavefunction(s) for diagram number 9577
    // (none)
    // Amplitude(s) for diagram number 9577
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[29], w_fp[2], w_fp[613], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[47] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[491] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[501] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[653] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9578 OF 15495 ***
    // Wavefunction(s) for diagram number 9578
    // (none)
    // Amplitude(s) for diagram number 9578
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[39], w_fp[245], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[45] += amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[339] -= amp_sv[0];
    jamp_sv[651] += amp_sv[0];

    // *** DIAGRAM 9579 OF 15495 ***
    // Wavefunction(s) for diagram number 9579
    // (none)
    // Amplitude(s) for diagram number 9579
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[39], w_fp[2], w_fp[612], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[45] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[427] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[469] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9580 OF 15495 ***
    // Wavefunction(s) for diagram number 9580
    // (none)
    // Amplitude(s) for diagram number 9580
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[534], w_fp[103], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[45] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[69] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[87] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[534], w_fp[103], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[45] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[69] += amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[349] -= amp_sv[0];
    jamp_sv[427] += amp_sv[0];
    jamp_sv[469] -= amp_sv[0];
    jamp_sv[491] -= amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[501] += amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[661] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[534], w_fp[103], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[83] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[349] -= amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[427] += amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[469] -= amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[491] -= amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[501] += amp_sv[0];
    jamp_sv[661] += amp_sv[0];
    jamp_sv[664] -= amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 245 );
#endif
#endif
  }

  //--------------------------------------------------------------------------


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup959( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 29 );
    retrieveWf( wfs, w_cx, nevt, 61 );
    retrieveWf( wfs, w_cx, nevt, 89 );
    retrieveWf( wfs, w_cx, nevt, 98 );
    retrieveWf( wfs, w_cx, nevt, 101 );
    retrieveWf( wfs, w_cx, nevt, 103 );
    retrieveWf( wfs, w_cx, nevt, 104 );
    retrieveWf( wfs, w_cx, nevt, 534 );
    retrieveWf( wfs, w_cx, nevt, 615 );
    retrieveWf( wfs, w_cx, nevt, 616 );
    retrieveWf( wfs, w_cx, nevt, 670 );
    retrieveWf( wfs, w_cx, nevt, 697 );
#endif
#endif

    // *** DIAGRAM 9581 OF 15495 ***
    // Wavefunction(s) for diagram number 9581
    // (none)
    // Amplitude(s) for diagram number 9581
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[103], w_fp[6], w_fp[697], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[45] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[69] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[87] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[664] += amp_sv[0];

    // *** DIAGRAM 9582 OF 15495 ***
    // Wavefunction(s) for diagram number 9582
    // (none)
    // Amplitude(s) for diagram number 9582
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[534], w_fp[6], w_fp[615], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[45] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[69] += amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[349] -= amp_sv[0];
    jamp_sv[427] += amp_sv[0];
    jamp_sv[469] -= amp_sv[0];
    jamp_sv[491] -= amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[501] += amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[661] += amp_sv[0];

    // *** DIAGRAM 9583 OF 15495 ***
    // Wavefunction(s) for diagram number 9583
    // (none)
    // Amplitude(s) for diagram number 9583
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[534], w_fp[104], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[47] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[69] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[534], w_fp[104], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[47] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[355] -= amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[547] += amp_sv[0];
    jamp_sv[589] -= amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[534], w_fp[104], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[59] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[355] -= amp_sv[0];
    jamp_sv[358] += amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[547] += amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[589] -= amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[670] -= amp_sv[0];

    // *** DIAGRAM 9584 OF 15495 ***
    // Wavefunction(s) for diagram number 9584
    // (none)
    // Amplitude(s) for diagram number 9584
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[104], w_fp[5], w_fp[697], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[47] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[69] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[670] += amp_sv[0];

    // *** DIAGRAM 9585 OF 15495 ***
    // Wavefunction(s) for diagram number 9585
    // (none)
    // Amplitude(s) for diagram number 9585
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[534], w_fp[5], w_fp[616], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[47] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[355] -= amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[547] += amp_sv[0];
    jamp_sv[589] -= amp_sv[0];
    jamp_sv[667] += amp_sv[0];

    // *** DIAGRAM 9586 OF 15495 ***
    // Wavefunction(s) for diagram number 9586
    // (none)
    // Amplitude(s) for diagram number 9586
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[534], w_fp[101], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[45] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[87] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[358] += amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[534], w_fp[89], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[358] += amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[534], w_fp[61], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[568] += amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[664] -= amp_sv[0];

    // *** DIAGRAM 9587 OF 15495 ***
    // Wavefunction(s) for diagram number 9587
    // (none)
    // Amplitude(s) for diagram number 9587
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[670], w_fp[98], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9588 OF 15495 ***
    // Wavefunction(s) for diagram number 9588
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[98], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[110] );
    // Amplitude(s) for diagram number 9588
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[29], w_fp[110], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[341] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[653] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9589 OF 15495 ***
    // Wavefunction(s) for diagram number 9589
    // (none)
    // Amplitude(s) for diagram number 9589
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[670], w_fp[2], w_fp[104], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[355] += amp_sv[0];
    jamp_sv[547] -= amp_sv[0];
    jamp_sv[589] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];

    // *** DIAGRAM 9590 OF 15495 ***
    // Wavefunction(s) for diagram number 9590
    // (none)
    // Amplitude(s) for diagram number 9590
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[29], w_fp[2], w_fp[616], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[47] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[93] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 110 );
#endif
#endif
  }

  //--------------------------------------------------------------------------


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup960( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 29 );
    retrieveWf( wfs, w_cx, nevt, 39 );
    retrieveWf( wfs, w_cx, nevt, 98 );
    retrieveWf( wfs, w_cx, nevt, 102 );
    retrieveWf( wfs, w_cx, nevt, 103 );
    retrieveWf( wfs, w_cx, nevt, 110 );
    retrieveWf( wfs, w_cx, nevt, 113 );
    retrieveWf( wfs, w_cx, nevt, 244 );
    retrieveWf( wfs, w_cx, nevt, 471 );
    retrieveWf( wfs, w_cx, nevt, 486 );
    retrieveWf( wfs, w_cx, nevt, 529 );
    retrieveWf( wfs, w_cx, nevt, 534 );
    retrieveWf( wfs, w_cx, nevt, 547 );
    retrieveWf( wfs, w_cx, nevt, 562 );
    retrieveWf( wfs, w_cx, nevt, 582 );
    retrieveWf( wfs, w_cx, nevt, 615 );
    retrieveWf( wfs, w_cx, nevt, 662 );
    retrieveWf( wfs, w_cx, nevt, 673 );
#endif
#endif

    // *** DIAGRAM 9591 OF 15495 ***
    // Wavefunction(s) for diagram number 9591
    // (none)
    // Amplitude(s) for diagram number 9591
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[673], w_fp[98], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9592 OF 15495 ***
    // Wavefunction(s) for diagram number 9592
    // (none)
    // Amplitude(s) for diagram number 9592
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[39], w_fp[110], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9593 OF 15495 ***
    // Wavefunction(s) for diagram number 9593
    // (none)
    // Amplitude(s) for diagram number 9593
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[673], w_fp[2], w_fp[103], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[349] += amp_sv[0];
    jamp_sv[427] -= amp_sv[0];
    jamp_sv[469] += amp_sv[0];
    jamp_sv[661] -= amp_sv[0];

    // *** DIAGRAM 9594 OF 15495 ***
    // Wavefunction(s) for diagram number 9594
    // (none)
    // Amplitude(s) for diagram number 9594
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[39], w_fp[2], w_fp[615], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[45] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[427] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[469] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9595 OF 15495 ***
    // Wavefunction(s) for diagram number 9595
    // (none)
    // Amplitude(s) for diagram number 9595
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[562], w_fp[534], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[45] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[69] += amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[349] -= amp_sv[0];
    jamp_sv[427] += amp_sv[0];
    jamp_sv[469] -= amp_sv[0];
    jamp_sv[491] -= amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[501] += amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[661] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[471], w_fp[534], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[69] += amp_sv[0];
    jamp_sv[339] += amp_sv[0];
    jamp_sv[349] -= amp_sv[0];
    jamp_sv[371] += amp_sv[0];
    jamp_sv[381] -= amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[547] -= amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[555] -= amp_sv[0];
    jamp_sv[557] += amp_sv[0];
    jamp_sv[589] += amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[651] -= amp_sv[0];
    jamp_sv[661] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[582], w_fp[534], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    jamp_sv[339] += amp_sv[0];
    jamp_sv[371] += amp_sv[0];
    jamp_sv[381] -= amp_sv[0];
    jamp_sv[427] -= amp_sv[0];
    jamp_sv[469] += amp_sv[0];
    jamp_sv[491] += amp_sv[0];
    jamp_sv[501] -= amp_sv[0];
    jamp_sv[547] -= amp_sv[0];
    jamp_sv[555] -= amp_sv[0];
    jamp_sv[557] += amp_sv[0];
    jamp_sv[568] += amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[589] += amp_sv[0];
    jamp_sv[651] -= amp_sv[0];

    // *** DIAGRAM 9596 OF 15495 ***
    // Wavefunction(s) for diagram number 9596
    // (none)
    // Amplitude(s) for diagram number 9596
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[39], w_fp[2], w_fp[562], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[45] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[427] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[469] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[39], w_fp[2], w_fp[471], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[39], w_fp[2], w_fp[582], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[45] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[427] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[469] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9597 OF 15495 ***
    // Wavefunction(s) for diagram number 9597
    // (none)
    // Amplitude(s) for diagram number 9597
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[486], w_fp[534], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[47] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[355] -= amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[547] += amp_sv[0];
    jamp_sv[589] -= amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[547], w_fp[534], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[341] += amp_sv[0];
    jamp_sv[355] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[427] -= amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[469] += amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[491] += amp_sv[0];
    jamp_sv[501] -= amp_sv[0];
    jamp_sv[653] -= amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[529], w_fp[534], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[341] += amp_sv[0];
    jamp_sv[371] += amp_sv[0];
    jamp_sv[381] -= amp_sv[0];
    jamp_sv[427] -= amp_sv[0];
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[469] += amp_sv[0];
    jamp_sv[491] += amp_sv[0];
    jamp_sv[501] -= amp_sv[0];
    jamp_sv[547] -= amp_sv[0];
    jamp_sv[589] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];

    // *** DIAGRAM 9598 OF 15495 ***
    // Wavefunction(s) for diagram number 9598
    // (none)
    // Amplitude(s) for diagram number 9598
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[29], w_fp[2], w_fp[486], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[47] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[93] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[29], w_fp[2], w_fp[547], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[83] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[93] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[491] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[501] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[653] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[29], w_fp[2], w_fp[529], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[47] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[491] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[501] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[653] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9599 OF 15495 ***
    // Wavefunction(s) for diagram number 9599
    // (none)
    // Amplitude(s) for diagram number 9599
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[98], w_fp[113], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[352] += amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[670] += amp_sv[0];

    // *** DIAGRAM 9600 OF 15495 ***
    // Wavefunction(s) for diagram number 9600
    // (none)
    // Amplitude(s) for diagram number 9600
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[244], w_fp[102], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[448] += amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[574] += amp_sv[0];

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
