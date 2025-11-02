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
  diagramgroup991( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 169 );
    retrieveWf( wfs, w_cx, nevt, 174 );
    retrieveWf( wfs, w_cx, nevt, 176 );
    retrieveWf( wfs, w_cx, nevt, 191 );
    retrieveWf( wfs, w_cx, nevt, 193 );
    retrieveWf( wfs, w_cx, nevt, 194 );
    retrieveWf( wfs, w_cx, nevt, 203 );
    retrieveWf( wfs, w_cx, nevt, 258 );
    retrieveWf( wfs, w_cx, nevt, 265 );
    retrieveWf( wfs, w_cx, nevt, 272 );
    retrieveWf( wfs, w_cx, nevt, 273 );
    retrieveWf( wfs, w_cx, nevt, 274 );
    retrieveWf( wfs, w_cx, nevt, 283 );
    retrieveWf( wfs, w_cx, nevt, 662 );
    retrieveWf( wfs, w_cx, nevt, 676 );
    retrieveWf( wfs, w_cx, nevt, 709 );
    retrieveWf( wfs, w_cx, nevt, 710 );
    retrieveWf( wfs, w_cx, nevt, 711 );
    retrieveWf( wfs, w_cx, nevt, 712 );
    retrieveWf( wfs, w_cx, nevt, 713 );
    retrieveWf( wfs, w_cx, nevt, 714 );
#endif
#endif

    // *** DIAGRAM 9901 OF 15495 ***
    // Wavefunction(s) for diagram number 9901
    // (none)
    // Amplitude(s) for diagram number 9901
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[193], w_fp[709], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[456] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[458] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[460] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[471] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[477] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[193], w_fp[710], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[460] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[464] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[471] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[475] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[193], w_fp[711], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[456] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[458] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[464] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[475] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[477] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9902 OF 15495 ***
    // Wavefunction(s) for diagram number 9902
    // (none)
    // Amplitude(s) for diagram number 9902
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[712], w_fp[194], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] += amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[380] += amp_sv[0];
    jamp_sv[394] -= amp_sv[0];
    jamp_sv[418] += amp_sv[0];
    jamp_sv[432] -= amp_sv[0];
    jamp_sv[434] += amp_sv[0];
    jamp_sv[436] += amp_sv[0];
    jamp_sv[437] -= amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[464] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[713], w_fp[194], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[380] += amp_sv[0];
    jamp_sv[391] += amp_sv[0];
    jamp_sv[394] -= amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[418] += amp_sv[0];
    jamp_sv[436] += amp_sv[0];
    jamp_sv[437] -= amp_sv[0];
    jamp_sv[440] -= amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[451] += amp_sv[0];
    jamp_sv[456] += amp_sv[0];
    jamp_sv[458] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[714], w_fp[194], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[391] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[432] += amp_sv[0];
    jamp_sv[434] -= amp_sv[0];
    jamp_sv[440] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[451] += amp_sv[0];
    jamp_sv[453] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[456] += amp_sv[0];
    jamp_sv[458] -= amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];

    // *** DIAGRAM 9903 OF 15495 ***
    // Wavefunction(s) for diagram number 9903
    // (none)
    // Amplitude(s) for diagram number 9903
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[191], w_fp[712], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[432] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[434] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[436] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[437] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[453] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[191], w_fp[713], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[436] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[437] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[451] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[191], w_fp[714], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[432] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[434] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[451] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[453] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9904 OF 15495 ***
    // Wavefunction(s) for diagram number 9904
    // (none)
    // Amplitude(s) for diagram number 9904
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[272], w_fp[194], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] -= amp_sv[0];
    jamp_sv[361] += amp_sv[0];
    jamp_sv[366] += amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[382] -= amp_sv[0];
    jamp_sv[383] += amp_sv[0];
    jamp_sv[393] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[417] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[477] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[273], w_fp[194], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] += amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[372] -= amp_sv[0];
    jamp_sv[374] += amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[380] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[274], w_fp[194], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[372] -= amp_sv[0];
    jamp_sv[374] += amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[380] += amp_sv[0];
    jamp_sv[382] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];

    // *** DIAGRAM 9905 OF 15495 ***
    // Wavefunction(s) for diagram number 9905
    // (none)
    // Amplitude(s) for diagram number 9905
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[203], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[456] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[458] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[464] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9906 OF 15495 ***
    // Wavefunction(s) for diagram number 9906
    // (none)
    // Amplitude(s) for diagram number 9906
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[174], w_fp[193], w_fp[662], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[456] += amp_sv[0];
    jamp_sv[458] -= amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];

    // *** DIAGRAM 9907 OF 15495 ***
    // Wavefunction(s) for diagram number 9907
    // (none)
    // Amplitude(s) for diagram number 9907
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[176], w_fp[169], w_fp[662], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] += amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[415] += amp_sv[0];

    // *** DIAGRAM 9908 OF 15495 ***
    // Wavefunction(s) for diagram number 9908
    // (none)
    // Amplitude(s) for diagram number 9908
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[283], w_fp[676], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[378] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9909 OF 15495 ***
    // Wavefunction(s) for diagram number 9909
    // (none)
    // Amplitude(s) for diagram number 9909
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[174], w_fp[676], w_fp[265], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] += amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[380] += amp_sv[0];

    // *** DIAGRAM 9910 OF 15495 ***
    // Wavefunction(s) for diagram number 9910
    // (none)
    // Amplitude(s) for diagram number 9910
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[176], w_fp[676], w_fp[258], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];

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
  diagramgroup992( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 169 );
    retrieveWf( wfs, w_cx, nevt, 174 );
    retrieveWf( wfs, w_cx, nevt, 176 );
    retrieveWf( wfs, w_cx, nevt, 179 );
    retrieveWf( wfs, w_cx, nevt, 181 );
    retrieveWf( wfs, w_cx, nevt, 191 );
    retrieveWf( wfs, w_cx, nevt, 193 );
    retrieveWf( wfs, w_cx, nevt, 203 );
    retrieveWf( wfs, w_cx, nevt, 205 );
    retrieveWf( wfs, w_cx, nevt, 258 );
    retrieveWf( wfs, w_cx, nevt, 265 );
    retrieveWf( wfs, w_cx, nevt, 283 );
    retrieveWf( wfs, w_cx, nevt, 489 );
    retrieveWf( wfs, w_cx, nevt, 594 );
    retrieveWf( wfs, w_cx, nevt, 662 );
    retrieveWf( wfs, w_cx, nevt, 712 );
    retrieveWf( wfs, w_cx, nevt, 713 );
    retrieveWf( wfs, w_cx, nevt, 714 );
#endif
#endif

    // *** DIAGRAM 9911 OF 15495 ***
    // Wavefunction(s) for diagram number 9911
    // (none)
    // Amplitude(s) for diagram number 9911
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[594], w_fp[489], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[394] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9912 OF 15495 ***
    // Wavefunction(s) for diagram number 9912
    // (none)
    // Amplitude(s) for diagram number 9912
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[594], w_fp[169], w_fp[265], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[394] += amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];

    // *** DIAGRAM 9913 OF 15495 ***
    // Wavefunction(s) for diagram number 9913
    // (none)
    // Amplitude(s) for diagram number 9913
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[594], w_fp[193], w_fp[258], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[464] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9914 OF 15495 ***
    // Wavefunction(s) for diagram number 9914
    // (none)
    // Amplitude(s) for diagram number 9914
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[176], w_fp[489], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[391] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9915 OF 15495 ***
    // Wavefunction(s) for diagram number 9915
    // (none)
    // Amplitude(s) for diagram number 9915
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[283], w_fp[193], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[456] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[458] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9916 OF 15495 ***
    // Wavefunction(s) for diagram number 9916
    // (none)
    // Amplitude(s) for diagram number 9916
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[265], w_fp[203], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[464] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9917 OF 15495 ***
    // Wavefunction(s) for diagram number 9917
    // (none)
    // Amplitude(s) for diagram number 9917
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[174], w_fp[169], w_fp[712], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[464] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[174], w_fp[169], w_fp[713], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[378] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[456] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[458] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[174], w_fp[169], w_fp[714], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[456] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[458] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[464] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9918 OF 15495 ***
    // Wavefunction(s) for diagram number 9918
    // (none)
    // Amplitude(s) for diagram number 9918
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[205], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[390] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[432] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[434] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9919 OF 15495 ***
    // Wavefunction(s) for diagram number 9919
    // (none)
    // Amplitude(s) for diagram number 9919
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[179], w_fp[191], w_fp[662], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[432] += amp_sv[0];
    jamp_sv[434] -= amp_sv[0];
    jamp_sv[440] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];

    // *** DIAGRAM 9920 OF 15495 ***
    // Wavefunction(s) for diagram number 9920
    // (none)
    // Amplitude(s) for diagram number 9920
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[181], w_fp[169], w_fp[662], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[414] += amp_sv[0];

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
  diagramgroup993( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 169 );
    retrieveWf( wfs, w_cx, nevt, 179 );
    retrieveWf( wfs, w_cx, nevt, 181 );
    retrieveWf( wfs, w_cx, nevt, 191 );
    retrieveWf( wfs, w_cx, nevt, 205 );
    retrieveWf( wfs, w_cx, nevt, 258 );
    retrieveWf( wfs, w_cx, nevt, 263 );
    retrieveWf( wfs, w_cx, nevt, 285 );
    retrieveWf( wfs, w_cx, nevt, 489 );
    retrieveWf( wfs, w_cx, nevt, 596 );
    retrieveWf( wfs, w_cx, nevt, 676 );
    retrieveWf( wfs, w_cx, nevt, 709 );
    retrieveWf( wfs, w_cx, nevt, 710 );
    retrieveWf( wfs, w_cx, nevt, 711 );
#endif
#endif

    // *** DIAGRAM 9921 OF 15495 ***
    // Wavefunction(s) for diagram number 9921
    // (none)
    // Amplitude(s) for diagram number 9921
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[285], w_fp[676], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[372] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9922 OF 15495 ***
    // Wavefunction(s) for diagram number 9922
    // (none)
    // Amplitude(s) for diagram number 9922
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[179], w_fp[676], w_fp[263], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[372] -= amp_sv[0];
    jamp_sv[374] += amp_sv[0];

    // *** DIAGRAM 9923 OF 15495 ***
    // Wavefunction(s) for diagram number 9923
    // (none)
    // Amplitude(s) for diagram number 9923
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[181], w_fp[676], w_fp[258], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9924 OF 15495 ***
    // Wavefunction(s) for diagram number 9924
    // (none)
    // Amplitude(s) for diagram number 9924
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[596], w_fp[489], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9925 OF 15495 ***
    // Wavefunction(s) for diagram number 9925
    // (none)
    // Amplitude(s) for diagram number 9925
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[596], w_fp[169], w_fp[263], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[392] += amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[440] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];

    // *** DIAGRAM 9926 OF 15495 ***
    // Wavefunction(s) for diagram number 9926
    // (none)
    // Amplitude(s) for diagram number 9926
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[596], w_fp[191], w_fp[258], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[440] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9927 OF 15495 ***
    // Wavefunction(s) for diagram number 9927
    // (none)
    // Amplitude(s) for diagram number 9927
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[181], w_fp[489], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[390] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9928 OF 15495 ***
    // Wavefunction(s) for diagram number 9928
    // (none)
    // Amplitude(s) for diagram number 9928
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[285], w_fp[191], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[432] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[434] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9929 OF 15495 ***
    // Wavefunction(s) for diagram number 9929
    // (none)
    // Amplitude(s) for diagram number 9929
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[263], w_fp[205], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9930 OF 15495 ***
    // Wavefunction(s) for diagram number 9930
    // (none)
    // Amplitude(s) for diagram number 9930
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[179], w_fp[169], w_fp[709], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[179], w_fp[169], w_fp[710], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[372] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[390] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[432] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[434] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[179], w_fp[169], w_fp[711], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[390] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[432] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[434] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] -= cxtype( 0, 1 ) * amp_sv[0];

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
  diagramgroup994( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 84 );
    retrieveWf( wfs, w_cx, nevt, 169 );
    retrieveWf( wfs, w_cx, nevt, 188 );
    retrieveWf( wfs, w_cx, nevt, 194 );
    retrieveWf( wfs, w_cx, nevt, 210 );
    retrieveWf( wfs, w_cx, nevt, 258 );
    retrieveWf( wfs, w_cx, nevt, 275 );
    retrieveWf( wfs, w_cx, nevt, 288 );
    retrieveWf( wfs, w_cx, nevt, 489 );
    retrieveWf( wfs, w_cx, nevt, 516 );
    retrieveWf( wfs, w_cx, nevt, 662 );
    retrieveWf( wfs, w_cx, nevt, 676 );
#endif
#endif

    // *** DIAGRAM 9931 OF 15495 ***
    // Wavefunction(s) for diagram number 9931
    // (none)
    // Amplitude(s) for diagram number 9931
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[194], w_fp[84], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[391] += amp_sv[0];
    jamp_sv[414] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[451] += amp_sv[0];
    jamp_sv[453] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[475] -= amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];

    // *** DIAGRAM 9932 OF 15495 ***
    // Wavefunction(s) for diagram number 9932
    // (none)
    // Amplitude(s) for diagram number 9932
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[210], w_fp[662], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[450] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[451] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[453] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[475] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[477] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9933 OF 15495 ***
    // Wavefunction(s) for diagram number 9933
    // (none)
    // Amplitude(s) for diagram number 9933
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[188], w_fp[169], w_fp[662], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[390] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9934 OF 15495 ***
    // Wavefunction(s) for diagram number 9934
    // (none)
    // Amplitude(s) for diagram number 9934
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[275], w_fp[676], w_fp[84], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[376] += amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[382] -= amp_sv[0];
    jamp_sv[383] += amp_sv[0];

    // *** DIAGRAM 9935 OF 15495 ***
    // Wavefunction(s) for diagram number 9935
    // (none)
    // Amplitude(s) for diagram number 9935
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[676], w_fp[288], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9936 OF 15495 ***
    // Wavefunction(s) for diagram number 9936
    // (none)
    // Amplitude(s) for diagram number 9936
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[188], w_fp[676], w_fp[258], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[367] += amp_sv[0];

    // *** DIAGRAM 9937 OF 15495 ***
    // Wavefunction(s) for diagram number 9937
    // (none)
    // Amplitude(s) for diagram number 9937
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[489], w_fp[516], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[390] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[419] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9938 OF 15495 ***
    // Wavefunction(s) for diagram number 9938
    // (none)
    // Amplitude(s) for diagram number 9938
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[275], w_fp[169], w_fp[516], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[376] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[451] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[475] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9939 OF 15495 ***
    // Wavefunction(s) for diagram number 9939
    // (none)
    // Amplitude(s) for diagram number 9939
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[516], w_fp[258], w_fp[194], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[376] -= amp_sv[0];
    jamp_sv[377] += amp_sv[0];
    jamp_sv[382] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[390] += amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[395] += amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[415] += amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[419] -= amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[451] -= amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[475] += amp_sv[0];

    // *** DIAGRAM 9940 OF 15495 ***
    // Wavefunction(s) for diagram number 9940
    // (none)
    // Amplitude(s) for diagram number 9940
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[188], w_fp[489], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[390] += amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[415] += amp_sv[0];

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
  diagramgroup995( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 169 );
    retrieveWf( wfs, w_cx, nevt, 194 );
    retrieveWf( wfs, w_cx, nevt, 210 );
    retrieveWf( wfs, w_cx, nevt, 212 );
    retrieveWf( wfs, w_cx, nevt, 213 );
    retrieveWf( wfs, w_cx, nevt, 214 );
    retrieveWf( wfs, w_cx, nevt, 275 );
    retrieveWf( wfs, w_cx, nevt, 288 );
    retrieveWf( wfs, w_cx, nevt, 544 );
    retrieveWf( wfs, w_cx, nevt, 662 );
    retrieveWf( wfs, w_cx, nevt, 680 );
    retrieveWf( wfs, w_cx, nevt, 724 );
    retrieveWf( wfs, w_cx, nevt, 725 );
    retrieveWf( wfs, w_cx, nevt, 726 );
    retrieveWf( wfs, w_cx, nevt, 741 );
#endif
#endif

    // *** DIAGRAM 9941 OF 15495 ***
    // Wavefunction(s) for diagram number 9941
    // (none)
    // Amplitude(s) for diagram number 9941
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[275], w_fp[210], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[450] += amp_sv[0];
    jamp_sv[451] -= amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[475] += amp_sv[0];

    // *** DIAGRAM 9942 OF 15495 ***
    // Wavefunction(s) for diagram number 9942
    // (none)
    // Amplitude(s) for diagram number 9942
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[288], w_fp[194], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] -= amp_sv[0];
    jamp_sv[361] += amp_sv[0];
    jamp_sv[366] += amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[382] -= amp_sv[0];
    jamp_sv[383] += amp_sv[0];
    jamp_sv[393] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[417] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[477] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];

    // *** DIAGRAM 9943 OF 15495 ***
    // Wavefunction(s) for diagram number 9943
    // (none)
    // Amplitude(s) for diagram number 9943
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[169], w_fp[726], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] -= amp_sv[0];
    jamp_sv[361] += amp_sv[0];
    jamp_sv[366] += amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[382] -= amp_sv[0];
    jamp_sv[383] += amp_sv[0];
    jamp_sv[393] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[417] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[477] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[169], w_fp[725], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[376] += amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[382] -= amp_sv[0];
    jamp_sv[383] += amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[391] += amp_sv[0];
    jamp_sv[393] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[414] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[417] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[451] += amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[475] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[169], w_fp[724], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[391] += amp_sv[0];
    jamp_sv[414] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[451] += amp_sv[0];
    jamp_sv[453] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[475] -= amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];

    // *** DIAGRAM 9944 OF 15495 ***
    // Wavefunction(s) for diagram number 9944
    // (none)
    // Amplitude(s) for diagram number 9944
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[214], w_fp[5], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[487] += amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[511] += amp_sv[0];
    jamp_sv[534] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[571] += amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[595] -= amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[214], w_fp[5], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[534] += amp_sv[0];
    jamp_sv[552] -= amp_sv[0];
    jamp_sv[554] += amp_sv[0];
    jamp_sv[560] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[595] -= amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[214], w_fp[5], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[535] += amp_sv[0];
    jamp_sv[552] -= amp_sv[0];
    jamp_sv[554] += amp_sv[0];
    jamp_sv[560] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[571] -= amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];

    // *** DIAGRAM 9945 OF 15495 ***
    // Wavefunction(s) for diagram number 9945
    // (none)
    // Amplitude(s) for diagram number 9945
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[214], w_fp[7], w_fp[544], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[534] += amp_sv[0];
    jamp_sv[552] -= amp_sv[0];
    jamp_sv[554] += amp_sv[0];
    jamp_sv[560] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[595] -= amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];

    // *** DIAGRAM 9946 OF 15495 ***
    // Wavefunction(s) for diagram number 9946
    // (none)
    // Amplitude(s) for diagram number 9946
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[214], w_fp[5], w_fp[680], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[535] += amp_sv[0];
    jamp_sv[552] -= amp_sv[0];
    jamp_sv[554] += amp_sv[0];
    jamp_sv[560] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[571] -= amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];

    // *** DIAGRAM 9947 OF 15495 ***
    // Wavefunction(s) for diagram number 9947
    // (none)
    // Amplitude(s) for diagram number 9947
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[741], w_fp[212], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[570] += amp_sv[0];
    jamp_sv[571] -= amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];

    // *** DIAGRAM 9948 OF 15495 ***
    // Wavefunction(s) for diagram number 9948
    // (none)
    // Amplitude(s) for diagram number 9948
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[212], w_fp[680], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[552] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[571] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9949 OF 15495 ***
    // Wavefunction(s) for diagram number 9949
    // (none)
    // Amplitude(s) for diagram number 9949
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[741], w_fp[213], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[594] += amp_sv[0];
    jamp_sv[595] -= amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];

    // *** DIAGRAM 9950 OF 15495 ***
    // Wavefunction(s) for diagram number 9950
    // (none)
    // Amplitude(s) for diagram number 9950
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[213], w_fp[544], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[576] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] += cxtype( 0, 1 ) * amp_sv[0];

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
  diagramgroup996( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 212 );
    retrieveWf( wfs, w_cx, nevt, 213 );
    retrieveWf( wfs, w_cx, nevt, 260 );
    retrieveWf( wfs, w_cx, nevt, 265 );
    retrieveWf( wfs, w_cx, nevt, 269 );
    retrieveWf( wfs, w_cx, nevt, 270 );
    retrieveWf( wfs, w_cx, nevt, 271 );
    retrieveWf( wfs, w_cx, nevt, 275 );
    retrieveWf( wfs, w_cx, nevt, 540 );
    retrieveWf( wfs, w_cx, nevt, 684 );
    retrieveWf( wfs, w_cx, nevt, 685 );
    retrieveWf( wfs, w_cx, nevt, 687 );
    retrieveWf( wfs, w_cx, nevt, 743 );
#endif
#endif

    // *** DIAGRAM 9951 OF 15495 ***
    // Wavefunction(s) for diagram number 9951
    // (none)
    // Amplitude(s) for diagram number 9951
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[275], w_fp[684], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[496] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9952 OF 15495 ***
    // Wavefunction(s) for diagram number 9952
    // (none)
    // Amplitude(s) for diagram number 9952
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[275], w_fp[685], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[502] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9953 OF 15495 ***
    // Wavefunction(s) for diagram number 9953
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[540], COUPs[1], 1.0, 0., 0., w_fp[489] );
    // Amplitude(s) for diagram number 9953
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[260], w_fp[7], w_fp[489], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[486] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[498] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[500] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[502] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9954 OF 15495 ***
    // Wavefunction(s) for diagram number 9954
    // (none)
    // Amplitude(s) for diagram number 9954
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[685], w_fp[260], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[498] += amp_sv[0];
    jamp_sv[500] -= amp_sv[0];
    jamp_sv[502] -= amp_sv[0];
    jamp_sv[503] += amp_sv[0];

    // *** DIAGRAM 9955 OF 15495 ***
    // Wavefunction(s) for diagram number 9955
    // (none)
    // Amplitude(s) for diagram number 9955
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[265], w_fp[5], w_fp[489], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[498] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[500] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9956 OF 15495 ***
    // Wavefunction(s) for diagram number 9956
    // (none)
    // Amplitude(s) for diagram number 9956
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[684], w_fp[265], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[492] += amp_sv[0];
    jamp_sv[494] -= amp_sv[0];
    jamp_sv[496] -= amp_sv[0];
    jamp_sv[497] += amp_sv[0];

    // *** DIAGRAM 9957 OF 15495 ***
    // Wavefunction(s) for diagram number 9957
    // (none)
    // Amplitude(s) for diagram number 9957
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[540], w_fp[269], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[481] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[486] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[502] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[540], w_fp[270], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[498] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[500] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[540], w_fp[271], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[486] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[498] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[500] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[502] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9958 OF 15495 ***
    // Wavefunction(s) for diagram number 9958
    // (none)
    // Amplitude(s) for diagram number 9958
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[743], w_fp[212], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[570] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[571] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9959 OF 15495 ***
    // Wavefunction(s) for diagram number 9959
    // (none)
    // Amplitude(s) for diagram number 9959
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[275], w_fp[687], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[556] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9960 OF 15495 ***
    // Wavefunction(s) for diagram number 9960
    // (none)
    // Amplitude(s) for diagram number 9960
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[743], w_fp[213], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];

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
  diagramgroup997( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 212 );
    retrieveWf( wfs, w_cx, nevt, 213 );
    retrieveWf( wfs, w_cx, nevt, 214 );
    retrieveWf( wfs, w_cx, nevt, 260 );
    retrieveWf( wfs, w_cx, nevt, 265 );
    retrieveWf( wfs, w_cx, nevt, 275 );
    retrieveWf( wfs, w_cx, nevt, 688 );
    retrieveWf( wfs, w_cx, nevt, 700 );
    retrieveWf( wfs, w_cx, nevt, 705 );
#endif
#endif

    // *** DIAGRAM 9961 OF 15495 ***
    // Wavefunction(s) for diagram number 9961
    // (none)
    // Amplitude(s) for diagram number 9961
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[275], w_fp[688], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[580] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9962 OF 15495 ***
    // Wavefunction(s) for diagram number 9962
    // (none)
    // Amplitude(s) for diagram number 9962
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[260], w_fp[214], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[492] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[536] += amp_sv[0];
    jamp_sv[560] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    jamp_sv[580] += amp_sv[0];
    jamp_sv[581] -= amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[260], w_fp[214], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[492] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[502] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[260], w_fp[214], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[502] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[536] -= amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[560] -= amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[576] += amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[580] -= amp_sv[0];
    jamp_sv[581] += amp_sv[0];

    // *** DIAGRAM 9963 OF 15495 ***
    // Wavefunction(s) for diagram number 9963
    // (none)
    // Amplitude(s) for diagram number 9963
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[214], w_fp[7], w_fp[700], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[492] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[536] += amp_sv[0];
    jamp_sv[560] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    jamp_sv[580] += amp_sv[0];
    jamp_sv[581] -= amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];

    // *** DIAGRAM 9964 OF 15495 ***
    // Wavefunction(s) for diagram number 9964
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[214], COUPs[0], 1.0, 0., 0., w_fp[745] );
    // Amplitude(s) for diagram number 9964
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[260], w_fp[7], w_fp[745], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[492] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[502] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];

    // *** DIAGRAM 9965 OF 15495 ***
    // Wavefunction(s) for diagram number 9965
    // (none)
    // Amplitude(s) for diagram number 9965
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[213], w_fp[700], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[576] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9966 OF 15495 ***
    // Wavefunction(s) for diagram number 9966
    // (none)
    // Amplitude(s) for diagram number 9966
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[688], w_fp[260], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[576] += amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[580] -= amp_sv[0];
    jamp_sv[581] += amp_sv[0];

    // *** DIAGRAM 9967 OF 15495 ***
    // Wavefunction(s) for diagram number 9967
    // (none)
    // Amplitude(s) for diagram number 9967
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[265], w_fp[214], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[514] -= amp_sv[0];
    jamp_sv[538] += amp_sv[0];
    jamp_sv[552] -= amp_sv[0];
    jamp_sv[554] += amp_sv[0];
    jamp_sv[556] += amp_sv[0];
    jamp_sv[557] -= amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[265], w_fp[214], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[492] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[496] += amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[265], w_fp[214], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[492] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[496] += amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[552] += amp_sv[0];
    jamp_sv[554] -= amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[557] += amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];

    // *** DIAGRAM 9968 OF 15495 ***
    // Wavefunction(s) for diagram number 9968
    // (none)
    // Amplitude(s) for diagram number 9968
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[214], w_fp[5], w_fp[705], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[514] -= amp_sv[0];
    jamp_sv[538] += amp_sv[0];
    jamp_sv[552] -= amp_sv[0];
    jamp_sv[554] += amp_sv[0];
    jamp_sv[556] += amp_sv[0];
    jamp_sv[557] -= amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];

    // *** DIAGRAM 9969 OF 15495 ***
    // Wavefunction(s) for diagram number 9969
    // (none)
    // Amplitude(s) for diagram number 9969
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[265], w_fp[5], w_fp[745], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[492] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[496] += amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];

    // *** DIAGRAM 9970 OF 15495 ***
    // Wavefunction(s) for diagram number 9970
    // (none)
    // Amplitude(s) for diagram number 9970
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[212], w_fp[705], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[552] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] += cxtype( 0, 1 ) * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 745 );
#endif
#endif
  }

  //--------------------------------------------------------------------------


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup998( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 168 );
    retrieveWf( wfs, w_cx, nevt, 171 );
    retrieveWf( wfs, w_cx, nevt, 197 );
    retrieveWf( wfs, w_cx, nevt, 212 );
    retrieveWf( wfs, w_cx, nevt, 213 );
    retrieveWf( wfs, w_cx, nevt, 214 );
    retrieveWf( wfs, w_cx, nevt, 219 );
    retrieveWf( wfs, w_cx, nevt, 265 );
    retrieveWf( wfs, w_cx, nevt, 269 );
    retrieveWf( wfs, w_cx, nevt, 270 );
    retrieveWf( wfs, w_cx, nevt, 271 );
    retrieveWf( wfs, w_cx, nevt, 281 );
    retrieveWf( wfs, w_cx, nevt, 540 );
    retrieveWf( wfs, w_cx, nevt, 662 );
    retrieveWf( wfs, w_cx, nevt, 687 );
    retrieveWf( wfs, w_cx, nevt, 706 );
    retrieveWf( wfs, w_cx, nevt, 707 );
    retrieveWf( wfs, w_cx, nevt, 708 );
    retrieveWf( wfs, w_cx, nevt, 712 );
    retrieveWf( wfs, w_cx, nevt, 713 );
    retrieveWf( wfs, w_cx, nevt, 714 );
#endif
#endif

    // *** DIAGRAM 9971 OF 15495 ***
    // Wavefunction(s) for diagram number 9971
    // (none)
    // Amplitude(s) for diagram number 9971
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[687], w_fp[265], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[552] += amp_sv[0];
    jamp_sv[554] -= amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[557] += amp_sv[0];

    // *** DIAGRAM 9972 OF 15495 ***
    // Wavefunction(s) for diagram number 9972
    // (none)
    // Amplitude(s) for diagram number 9972
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[706], w_fp[214], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[492] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[536] += amp_sv[0];
    jamp_sv[560] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    jamp_sv[580] += amp_sv[0];
    jamp_sv[581] -= amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[707], w_fp[214], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[492] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[510] += amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[534] -= amp_sv[0];
    jamp_sv[536] += amp_sv[0];
    jamp_sv[552] += amp_sv[0];
    jamp_sv[554] -= amp_sv[0];
    jamp_sv[580] += amp_sv[0];
    jamp_sv[581] -= amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[595] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[708], w_fp[214], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[486] += amp_sv[0];
    jamp_sv[510] += amp_sv[0];
    jamp_sv[534] -= amp_sv[0];
    jamp_sv[552] += amp_sv[0];
    jamp_sv[554] -= amp_sv[0];
    jamp_sv[560] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    jamp_sv[576] += amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[595] += amp_sv[0];
    jamp_sv[597] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];

    // *** DIAGRAM 9973 OF 15495 ***
    // Wavefunction(s) for diagram number 9973
    // (none)
    // Amplitude(s) for diagram number 9973
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[213], w_fp[706], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[576] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[213], w_fp[707], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[580] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[213], w_fp[708], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[576] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9974 OF 15495 ***
    // Wavefunction(s) for diagram number 9974
    // (none)
    // Amplitude(s) for diagram number 9974
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[712], w_fp[214], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[514] -= amp_sv[0];
    jamp_sv[538] += amp_sv[0];
    jamp_sv[552] -= amp_sv[0];
    jamp_sv[554] += amp_sv[0];
    jamp_sv[556] += amp_sv[0];
    jamp_sv[557] -= amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[713], w_fp[214], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[511] += amp_sv[0];
    jamp_sv[514] -= amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[538] += amp_sv[0];
    jamp_sv[556] += amp_sv[0];
    jamp_sv[557] -= amp_sv[0];
    jamp_sv[560] -= amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[571] += amp_sv[0];
    jamp_sv[576] += amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[714], w_fp[214], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[487] += amp_sv[0];
    jamp_sv[511] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[552] += amp_sv[0];
    jamp_sv[554] -= amp_sv[0];
    jamp_sv[560] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[571] += amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[576] += amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];

    // *** DIAGRAM 9975 OF 15495 ***
    // Wavefunction(s) for diagram number 9975
    // (none)
    // Amplitude(s) for diagram number 9975
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[212], w_fp[712], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[552] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[212], w_fp[713], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[556] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[571] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[212], w_fp[714], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[552] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[571] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9976 OF 15495 ***
    // Wavefunction(s) for diagram number 9976
    // (none)
    // Amplitude(s) for diagram number 9976
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[269], w_fp[214], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[481] += amp_sv[0];
    jamp_sv[486] += amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[496] += amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[502] -= amp_sv[0];
    jamp_sv[503] += amp_sv[0];
    jamp_sv[513] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[537] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[597] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[270], w_fp[214], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[492] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[496] += amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[271], w_fp[214], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[492] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[502] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];

    // *** DIAGRAM 9977 OF 15495 ***
    // Wavefunction(s) for diagram number 9977
    // (none)
    // Amplitude(s) for diagram number 9977
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[219], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[576] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9978 OF 15495 ***
    // Wavefunction(s) for diagram number 9978
    // (none)
    // Amplitude(s) for diagram number 9978
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[168], w_fp[213], w_fp[662], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[576] += amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];

    // *** DIAGRAM 9979 OF 15495 ***
    // Wavefunction(s) for diagram number 9979
    // (none)
    // Amplitude(s) for diagram number 9979
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[171], w_fp[197], w_fp[662], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[535] += amp_sv[0];

    // *** DIAGRAM 9980 OF 15495 ***
    // Wavefunction(s) for diagram number 9980
    // (none)
    // Amplitude(s) for diagram number 9980
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[281], w_fp[540], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[498] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[500] -= cxtype( 0, 1 ) * amp_sv[0];

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
  diagramgroup999( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 168 );
    retrieveWf( wfs, w_cx, nevt, 171 );
    retrieveWf( wfs, w_cx, nevt, 197 );
    retrieveWf( wfs, w_cx, nevt, 213 );
    retrieveWf( wfs, w_cx, nevt, 219 );
    retrieveWf( wfs, w_cx, nevt, 220 );
    retrieveWf( wfs, w_cx, nevt, 258 );
    retrieveWf( wfs, w_cx, nevt, 265 );
    retrieveWf( wfs, w_cx, nevt, 281 );
    retrieveWf( wfs, w_cx, nevt, 491 );
    retrieveWf( wfs, w_cx, nevt, 527 );
    retrieveWf( wfs, w_cx, nevt, 540 );
    retrieveWf( wfs, w_cx, nevt, 662 );
    retrieveWf( wfs, w_cx, nevt, 712 );
    retrieveWf( wfs, w_cx, nevt, 713 );
    retrieveWf( wfs, w_cx, nevt, 714 );
#endif
#endif

    // *** DIAGRAM 9981 OF 15495 ***
    // Wavefunction(s) for diagram number 9981
    // (none)
    // Amplitude(s) for diagram number 9981
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[168], w_fp[540], w_fp[265], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];

    // *** DIAGRAM 9982 OF 15495 ***
    // Wavefunction(s) for diagram number 9982
    // (none)
    // Amplitude(s) for diagram number 9982
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[171], w_fp[540], w_fp[258], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9983 OF 15495 ***
    // Wavefunction(s) for diagram number 9983
    // (none)
    // Amplitude(s) for diagram number 9983
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[527], w_fp[491], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[514] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9984 OF 15495 ***
    // Wavefunction(s) for diagram number 9984
    // (none)
    // Amplitude(s) for diagram number 9984
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[527], w_fp[197], w_fp[265], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[514] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];

    // *** DIAGRAM 9985 OF 15495 ***
    // Wavefunction(s) for diagram number 9985
    // (none)
    // Amplitude(s) for diagram number 9985
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[527], w_fp[213], w_fp[258], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[584] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9986 OF 15495 ***
    // Wavefunction(s) for diagram number 9986
    // (none)
    // Amplitude(s) for diagram number 9986
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[171], w_fp[491], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[511] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9987 OF 15495 ***
    // Wavefunction(s) for diagram number 9987
    // (none)
    // Amplitude(s) for diagram number 9987
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[281], w_fp[213], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[576] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9988 OF 15495 ***
    // Wavefunction(s) for diagram number 9988
    // (none)
    // Amplitude(s) for diagram number 9988
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[265], w_fp[219], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[498] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[500] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9989 OF 15495 ***
    // Wavefunction(s) for diagram number 9989
    // (none)
    // Amplitude(s) for diagram number 9989
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[168], w_fp[197], w_fp[712], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[498] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[500] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[168], w_fp[197], w_fp[713], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[498] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[500] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[576] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[168], w_fp[197], w_fp[714], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[576] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9990 OF 15495 ***
    // Wavefunction(s) for diagram number 9990
    // (none)
    // Amplitude(s) for diagram number 9990
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[220], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[486] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] -= cxtype( 0, 1 ) * amp_sv[0];

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
  diagramgroup1000( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 179 );
    retrieveWf( wfs, w_cx, nevt, 180 );
    retrieveWf( wfs, w_cx, nevt, 197 );
    retrieveWf( wfs, w_cx, nevt, 212 );
    retrieveWf( wfs, w_cx, nevt, 258 );
    retrieveWf( wfs, w_cx, nevt, 260 );
    retrieveWf( wfs, w_cx, nevt, 285 );
    retrieveWf( wfs, w_cx, nevt, 491 );
    retrieveWf( wfs, w_cx, nevt, 540 );
    retrieveWf( wfs, w_cx, nevt, 596 );
    retrieveWf( wfs, w_cx, nevt, 662 );
#endif
#endif

    // *** DIAGRAM 9991 OF 15495 ***
    // Wavefunction(s) for diagram number 9991
    // (none)
    // Amplitude(s) for diagram number 9991
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[179], w_fp[212], w_fp[662], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[552] += amp_sv[0];
    jamp_sv[554] -= amp_sv[0];
    jamp_sv[560] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];

    // *** DIAGRAM 9992 OF 15495 ***
    // Wavefunction(s) for diagram number 9992
    // (none)
    // Amplitude(s) for diagram number 9992
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[180], w_fp[197], w_fp[662], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[534] += amp_sv[0];

    // *** DIAGRAM 9993 OF 15495 ***
    // Wavefunction(s) for diagram number 9993
    // (none)
    // Amplitude(s) for diagram number 9993
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[285], w_fp[540], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[492] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9994 OF 15495 ***
    // Wavefunction(s) for diagram number 9994
    // (none)
    // Amplitude(s) for diagram number 9994
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[179], w_fp[540], w_fp[260], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[492] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];

    // *** DIAGRAM 9995 OF 15495 ***
    // Wavefunction(s) for diagram number 9995
    // (none)
    // Amplitude(s) for diagram number 9995
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[180], w_fp[540], w_fp[258], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[486] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9996 OF 15495 ***
    // Wavefunction(s) for diagram number 9996
    // (none)
    // Amplitude(s) for diagram number 9996
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[596], w_fp[491], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9997 OF 15495 ***
    // Wavefunction(s) for diagram number 9997
    // (none)
    // Amplitude(s) for diagram number 9997
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[596], w_fp[197], w_fp[260], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[512] += amp_sv[0];
    jamp_sv[536] -= amp_sv[0];
    jamp_sv[560] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];

    // *** DIAGRAM 9998 OF 15495 ***
    // Wavefunction(s) for diagram number 9998
    // (none)
    // Amplitude(s) for diagram number 9998
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[596], w_fp[212], w_fp[258], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[560] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9999 OF 15495 ***
    // Wavefunction(s) for diagram number 9999
    // (none)
    // Amplitude(s) for diagram number 9999
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[180], w_fp[491], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10000 OF 15495 ***
    // Wavefunction(s) for diagram number 10000
    // (none)
    // Amplitude(s) for diagram number 10000
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[285], w_fp[212], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[552] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];

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
