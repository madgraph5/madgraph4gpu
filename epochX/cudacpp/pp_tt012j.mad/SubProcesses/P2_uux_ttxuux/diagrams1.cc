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
  constexpr int nw6 = CPPProcess::nw6;       // dimensions of each wavefunction (HELAS KEK 91-11): e.g. 6 for e+ e- -> mu+ mu- (fermions and vectors)
  constexpr int npar = CPPProcess::npar;     // #particles in total (external = initial + final): e.g. 4 for e+ e- -> mu+ mu-
  constexpr int nwf = CPPProcess::nwf;       // #wavefunctions = #external (npar) + #internal: e.g. 5 for e+ e- -> mu+ mu- (1 internal is gamma or Z)
  constexpr int ncolor = CPPProcess::ncolor; // the number of leading colors

  using Parameters_sm_dependentCouplings::ndcoup;   // #couplings that vary event by event (depend on running alphas QCD)
  using Parameters_sm_independentCouplings::nicoup; // #couplings that are fixed for all events (do not depend on running alphas QCD)

#ifdef __CUDACC__ // this must be __CUDACC__ (not MGONGPUCPP_GPUIMPL)
#pragma nv_diagnostic push
#pragma nv_diag_suppress 177 // e.g. <<warning #177-D: variable "nevt" was declared but never referenced>>
#endif
  constexpr int nIPD = CPPProcess::nIPD; // SM independent parameters
  constexpr int nIPC = CPPProcess::nIPC; // SM independent couplings
#ifdef __CUDACC__ // this must be __CUDACC__ (not MGONGPUCPP_GPUIMPL)
#pragma nv_diagnostic pop
#endif

  //--------------------------------------------------------------------------

  __global__ void
  diagramgroup1( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
#ifdef MGONGPUCPP_GPUIMPL
                 fptype* jamps,                  // output jamps[ncolor*2*nevt]
                 const fptype* couplings,        // input: dependent couplings[nevt*ndcoup*2] for all events
#else
                 cxtype_sv* jamp_sv,             // output jamps[ncolor*2*neppV]
                 const fptype** COUPs,           // input: dependent and independent COUPs[nxcoup] for this event page
#endif
                 const unsigned int* channelIds, // input: channelIds[nevt] for GPU or SCALAR channelId[0] for C++ (1 to #diagrams, 0 to disable SDE)
                 fptype* numerators,             // input/output: multichannel numerators[nevtORneppV], add helicity ihel
                 fptype* denominators,           // input/output: multichannel denominators[nevtORneppV], add helicity ihel
                 const fptype* cIPC,             // input: GPU __device__ or GPU host address of cIPC
                 const fptype* cIPD,             // input: GPU __device__ or GPU host address of cIPD
                 const short* cHelFlat,          // input: GPU __device__ or GPU host address of cHel
                 const fptype* momenta,          // input: momenta[npar*4*nevtORneppV]
                 const int ihel )                // input: helicity (0 to ncomb)
  {
    // A uniform interface for diagramgroupXXX including channelIDs, numerators and denominators is used also #ifndef MGONGPU_SUPPORTS_MULTICHANNEL
    // In that case, however, the boilerplate code asserts that all three pointers all nullptr as a sanity check
#include "diagrams_boilerplate.h"
#ifdef MGONGPUCPP_GPUIMPL
    using M_ACCESS = DeviceAccessMomenta; // non-trivial access: buffer includes all events
#else
    using M_ACCESS = HostAccessMomenta; // non-trivial access: buffer includes all events
#endif

#ifdef MGONGPUCPP_GPUIMPL
    // *** RETRIEVE WAVEFUNCTIONS FROM PREVIOUS DIAGRAM GROUPS ***
    // (none)
#endif

    // Reinterpret the flat array pointer for helicities as a multidimensional array pointer
    const short (*cHel)[npar] = reinterpret_cast<const short(*)[npar]>( cHelFlat );

    // *** DIAGRAM 1 OF 14 ***
    // Wavefunction(s) for diagram number 1
    ixxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][0], +1, w_fp[0], 0 );
    oxxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][1], -1, w_fp[1], 1 );
    oxxxxx<M_ACCESS, W_ACCESS>( momenta, cIPD[0], cHel[ihel][2], +1, w_fp[2], 2 );
    ixxxxx<M_ACCESS, W_ACCESS>( momenta, cIPD[0], cHel[ihel][3], -1, w_fp[3], 3 );
    oxxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][4], +1, w_fp[4], 4 );
    ixxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][5], -1, w_fp[5], 5 );
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], COUPs[1], 1.0, 0., 0., w_fp[6] );
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[5], w_fp[4], COUPs[1], 1.0, 0., 0., w_fp[7] );
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[6], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[8] );
    // Amplitude(s) for diagram number 1
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[8], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 1 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += 1. / 36. * amp_sv[0];
    jamp_sv[1] -= 1. / 12. * amp_sv[0];
    jamp_sv[2] -= 1. / 12. * amp_sv[0];
    jamp_sv[4] += 1. / 4. * amp_sv[0];

    // *** DIAGRAM 2 OF 14 ***
    // Wavefunction(s) for diagram number 2
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[6], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[8] );
    // Amplitude(s) for diagram number 2
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[2], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 2 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += 1. / 36. * amp_sv[0];
    jamp_sv[1] -= 1. / 12. * amp_sv[0];
    jamp_sv[2] -= 1. / 12. * amp_sv[0];
    jamp_sv[3] += 1. / 4. * amp_sv[0];

    // *** DIAGRAM 3 OF 14 ***
    // Wavefunction(s) for diagram number 3
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[2], COUPs[1], 1.0, 0., 0., w_fp[8] );
    // Amplitude(s) for diagram number 3
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[6], w_fp[7], w_fp[8], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 3 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[3] += 1. / 4. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] -= 1. / 4. * cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4 OF 14 ***
    // Wavefunction(s) for diagram number 4
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[5], w_fp[6], COUPs[1], 1.0, 0., 0., w_fp[9] );
    // Amplitude(s) for diagram number 4
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[4], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 4 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += 1. / 36. * amp_sv[0];
    jamp_sv[1] -= 1. / 12. * amp_sv[0];
    jamp_sv[4] += 1. / 4. * amp_sv[0];
    jamp_sv[5] -= 1. / 12. * amp_sv[0];

    // *** DIAGRAM 5 OF 14 ***
    // Wavefunction(s) for diagram number 5
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[4], w_fp[6], COUPs[1], 1.0, 0., 0., w_fp[9] );
    // Amplitude(s) for diagram number 5
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[9], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 5 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += 1. / 36. * amp_sv[0];
    jamp_sv[1] -= 1. / 12. * amp_sv[0];
    jamp_sv[3] += 1. / 4. * amp_sv[0];
    jamp_sv[5] -= 1. / 12. * amp_sv[0];

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE JAMPS ***
    for( int icol = 0; icol < ncolor; icol++ )
      J_ACCESS::kernelAccessIcol( jamps, icol ) = jamp_sv[icol]; // set jamps
#endif

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE WAVEFUNCTIONS FOR NEXT DIAGRAM GROUPS ***
    //for( int iwf = 0; iwf < nwf; iwf++ ) storeWf( wfs, w_cx, nevt, iwf );
    storeWf( wfs, w_cx, nevt, 0 );
    storeWf( wfs, w_cx, nevt, 1 );
    storeWf( wfs, w_cx, nevt, 2 );
    storeWf( wfs, w_cx, nevt, 3 );
    storeWf( wfs, w_cx, nevt, 4 );
    storeWf( wfs, w_cx, nevt, 5 );
    storeWf( wfs, w_cx, nevt, 6 );
    storeWf( wfs, w_cx, nevt, 7 );
    storeWf( wfs, w_cx, nevt, 8 );
    storeWf( wfs, w_cx, nevt, 9 );
#endif
  }

  //--------------------------------------------------------------------------

  __global__ void
  diagramgroup2( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
#ifdef MGONGPUCPP_GPUIMPL
                 fptype* jamps,                  // output jamps[ncolor*2*nevt]
                 const fptype* couplings,        // input: dependent couplings[nevt*ndcoup*2] for all events
#else
                 cxtype_sv* jamp_sv,             // output jamps[ncolor*2*neppV]
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
    // *** RETRIEVE WAVEFUNCTIONS FROM PREVIOUS DIAGRAM GROUPS ***
    //for( int iwf = 0; iwf < nwf; iwf++ ) retrieveWf( wfs, w_cx, nevt, iwf );
    retrieveWf( wfs, w_cx, nevt, 0 );
    retrieveWf( wfs, w_cx, nevt, 1 );
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 8 );
#endif

    // *** DIAGRAM 6 OF 14 ***
    // Wavefunction(s) for diagram number 6
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[4], COUPs[1], 1.0, 0., 0., w_fp[9] );
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[5], w_fp[1], COUPs[1], 1.0, 0., 0., w_fp[6] );
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[9], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[10] );
    // Amplitude(s) for diagram number 6
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[10], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 6 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[2] -= 1. / 4. * amp_sv[0];
    jamp_sv[3] += 1. / 12. * amp_sv[0];
    jamp_sv[4] += 1. / 12. * amp_sv[0];
    jamp_sv[5] -= 1. / 36. * amp_sv[0];

    // *** DIAGRAM 7 OF 14 ***
    // Wavefunction(s) for diagram number 7
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[9], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[10] );
    // Amplitude(s) for diagram number 7
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[10], w_fp[2], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 7 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[1] -= 1. / 4. * amp_sv[0];
    jamp_sv[3] += 1. / 12. * amp_sv[0];
    jamp_sv[4] += 1. / 12. * amp_sv[0];
    jamp_sv[5] -= 1. / 36. * amp_sv[0];

    // *** DIAGRAM 8 OF 14 ***
    // Wavefunction(s) for diagram number 8
    // (none)
    // Amplitude(s) for diagram number 8
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[6], w_fp[8], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 8 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[1] -= 1. / 4. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] += 1. / 4. * cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9 OF 14 ***
    // Wavefunction(s) for diagram number 9
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[5], w_fp[9], COUPs[1], 1.0, 0., 0., w_fp[10] );
    // Amplitude(s) for diagram number 9
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[10], w_fp[1], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 9 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += 1. / 12. * amp_sv[0];
    jamp_sv[2] -= 1. / 4. * amp_sv[0];
    jamp_sv[3] += 1. / 12. * amp_sv[0];
    jamp_sv[5] -= 1. / 36. * amp_sv[0];

    // *** DIAGRAM 10 OF 14 ***
    // Wavefunction(s) for diagram number 10
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[9], COUPs[1], 1.0, 0., 0., w_fp[10] );
    // Amplitude(s) for diagram number 10
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[10], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 10 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += 1. / 12. * amp_sv[0];
    jamp_sv[1] -= 1. / 4. * amp_sv[0];
    jamp_sv[3] += 1. / 12. * amp_sv[0];
    jamp_sv[5] -= 1. / 36. * amp_sv[0];

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE JAMPS ***
    for( int icol = 0; icol < ncolor; icol++ )
      J_ACCESS::kernelAccessIcol( jamps, icol ) += jamp_sv[icol]; // update jamps
#endif

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE WAVEFUNCTIONS FOR NEXT DIAGRAM GROUPS ***
    //for( int iwf = 0; iwf < nwf; iwf++ ) storeWf( wfs, w_cx, nevt, iwf );
    storeWf( wfs, w_cx, nevt, 6 );
    storeWf( wfs, w_cx, nevt, 9 );
    storeWf( wfs, w_cx, nevt, 10 );
#endif
  }

  //--------------------------------------------------------------------------

  __global__ void
  diagramgroup3( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
#ifdef MGONGPUCPP_GPUIMPL
                 fptype* jamps,                  // output jamps[ncolor*2*nevt]
                 const fptype* couplings,        // input: dependent couplings[nevt*ndcoup*2] for all events
#else
                 cxtype_sv* jamp_sv,             // output jamps[ncolor*2*neppV]
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
    // *** RETRIEVE WAVEFUNCTIONS FROM PREVIOUS DIAGRAM GROUPS ***
    //for( int iwf = 0; iwf < nwf; iwf++ ) retrieveWf( wfs, w_cx, nevt, iwf );
    retrieveWf( wfs, w_cx, nevt, 0 );
    retrieveWf( wfs, w_cx, nevt, 1 );
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 8 );
#endif

    // *** DIAGRAM 11 OF 14 ***
    // Wavefunction(s) for diagram number 11
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[6], COUPs[1], 1.0, 0., 0., w_fp[10] );
    // Amplitude(s) for diagram number 11
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[10], w_fp[4], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 11 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += 1. / 12. * amp_sv[0];
    jamp_sv[1] -= 1. / 4. * amp_sv[0];
    jamp_sv[4] += 1. / 12. * amp_sv[0];
    jamp_sv[5] -= 1. / 36. * amp_sv[0];

    // *** DIAGRAM 12 OF 14 ***
    // Wavefunction(s) for diagram number 12
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[8], COUPs[1], 1.0, 0., 0., w_fp[10] );
    // Amplitude(s) for diagram number 12
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[10], w_fp[4], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 12 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += 1. / 12. * amp_sv[0];
    jamp_sv[2] -= 1. / 4. * amp_sv[0];
    jamp_sv[4] += 1. / 12. * amp_sv[0];
    jamp_sv[5] -= 1. / 36. * amp_sv[0];

    // *** DIAGRAM 13 OF 14 ***
    // Wavefunction(s) for diagram number 13
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[7], COUPs[1], 1.0, 0., 0., w_fp[6] );
    // Amplitude(s) for diagram number 13
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[6], w_fp[1], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 13 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += 1. / 36. * amp_sv[0];
    jamp_sv[2] -= 1. / 12. * amp_sv[0];
    jamp_sv[3] += 1. / 4. * amp_sv[0];
    jamp_sv[5] -= 1. / 12. * amp_sv[0];

    // *** DIAGRAM 14 OF 14 ***
    // Wavefunction(s) for diagram number 14
    // (none)
    // Amplitude(s) for diagram number 14
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[10], w_fp[1], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 14 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += 1. / 36. * amp_sv[0];
    jamp_sv[2] -= 1. / 12. * amp_sv[0];
    jamp_sv[4] += 1. / 4. * amp_sv[0];
    jamp_sv[5] -= 1. / 12. * amp_sv[0];

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE JAMPS ***
    for( int icol = 0; icol < ncolor; icol++ )
      J_ACCESS::kernelAccessIcol( jamps, icol ) += jamp_sv[icol]; // update jamps
#endif

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE WAVEFUNCTIONS FOR NEXT DIAGRAM GROUPS ***
    // (none)
#endif
  }

  //--------------------------------------------------------------------------
}
