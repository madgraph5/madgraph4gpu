// Copyright (C) 2020-2025 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Sep 2025) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2025) for the MG5aMC CUDACPP plugin.

#include "GpuRuntime.h"
#include "HelAmps_SMEFTsim_topU3l_MwScheme_UFO.h"
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

  using Parameters_SMEFTsim_topU3l_MwScheme_UFO_dependentCouplings::ndcoup;   // #couplings that vary event by event (depend on running alphas QCD)
  using Parameters_SMEFTsim_topU3l_MwScheme_UFO_independentCouplings::nicoup; // #couplings that are fixed for all events (do not depend on running alphas QCD)

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

    // *** DIAGRAM 1 OF 72 ***
    // Wavefunction(s) for diagram number 1
    vxxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][0], -1, w_fp[0], 0 );
    vxxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][1], -1, w_fp[1], 1 );
    oxxxxx<M_ACCESS, W_ACCESS>( momenta, cIPD[0], cHel[ihel][2], +1, w_fp[2], 2 );
    ixxxxx<M_ACCESS, W_ACCESS>( momenta, cIPD[0], cHel[ihel][3], -1, w_fp[3], 3 );
    oxxxxx<M_ACCESS, W_ACCESS>( momenta, cIPD[0], cHel[ihel][4], +1, w_fp[4], 4 );
    ixxxxx<M_ACCESS, W_ACCESS>( momenta, cIPD[0], cHel[ihel][5], -1, w_fp[5], 5 );
    VVV5P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], COUPs[0], 1.0, 0., 0., w_fp[6] );
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[2], COUPs[1], 1.0, 0., 0., w_fp[7] );
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[4], w_fp[6], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[8] );
    // Amplitude(s) for diagram number 1
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[8], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[3] += 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] -= 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2 OF 72 ***
    // Wavefunction(s) for diagram number 2
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[5], w_fp[6], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[9] );
    // Amplitude(s) for diagram number 2
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[4], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[3] += 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] -= 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3 OF 72 ***
    // Wavefunction(s) for diagram number 3
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[5], w_fp[4], COUPs[1], 1.0, 0., 0., w_fp[10] );
    // Amplitude(s) for diagram number 3
    VVV5_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[6], w_fp[7], w_fp[10], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += 1. / 2. * amp_sv[0];
    jamp_sv[2] -= 1. / 2. * amp_sv[0];
    jamp_sv[9] -= 1. / 2. * amp_sv[0];
    jamp_sv[10] += 1. / 2. * amp_sv[0];

    // *** DIAGRAM 4 OF 72 ***
    // Wavefunction(s) for diagram number 4
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[5], w_fp[2], COUPs[1], 1.0, 0., 0., w_fp[11] );
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[6], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[12] );
    // Amplitude(s) for diagram number 4
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[4], w_fp[11], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] -= 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] += 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5 OF 72 ***
    // Wavefunction(s) for diagram number 5
    // (none)
    // Amplitude(s) for diagram number 5
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[8], w_fp[11], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] -= 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[3] += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] += 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 10 );
    storeWf( wfs, w_cx, nevt, 11 );
    storeWf( wfs, w_cx, nevt, 12 );
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
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 9 );
    retrieveWf( wfs, w_cx, nevt, 10 );
    retrieveWf( wfs, w_cx, nevt, 11 );
    retrieveWf( wfs, w_cx, nevt, 12 );
#endif

    // *** DIAGRAM 6 OF 72 ***
    // Wavefunction(s) for diagram number 6
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[4], COUPs[1], 1.0, 0., 0., w_fp[8] );
    // Amplitude(s) for diagram number 6
    VVV5_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[6], w_fp[11], w_fp[8], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= 1. / 2. * amp_sv[0];
    jamp_sv[3] += 1. / 2. * amp_sv[0];
    jamp_sv[8] += 1. / 2. * amp_sv[0];
    jamp_sv[11] -= 1. / 2. * amp_sv[0];

    // *** DIAGRAM 7 OF 72 ***
    // Wavefunction(s) for diagram number 7
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[6], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[13] );
    // Amplitude(s) for diagram number 7
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[13], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[1] -= 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] += 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8 OF 72 ***
    // Wavefunction(s) for diagram number 8
    // (none)
    // Amplitude(s) for diagram number 8
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[2], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[3] += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] += 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9 OF 72 ***
    // Wavefunction(s) for diagram number 9
    // (none)
    // Amplitude(s) for diagram number 9
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[13], w_fp[10], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[1] -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] -= 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10 OF 72 ***
    // Wavefunction(s) for diagram number 10
    // (none)
    // Amplitude(s) for diagram number 10
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[2], w_fp[10], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] -= 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE JAMPS ***
    for( int icol = 0; icol < ncolor; icol++ )
      J_ACCESS::kernelAccessIcol( jamps, icol ) += jamp_sv[icol]; // update jamps
#endif

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE WAVEFUNCTIONS FOR NEXT DIAGRAM GROUPS ***
    //for( int iwf = 0; iwf < nwf; iwf++ ) storeWf( wfs, w_cx, nevt, iwf );
    storeWf( wfs, w_cx, nevt, 8 );
    storeWf( wfs, w_cx, nevt, 13 );
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
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 10 );
#endif

    // *** DIAGRAM 11 OF 72 ***
    // Wavefunction(s) for diagram number 11
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[12] );
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[13] );
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[5], w_fp[12], COUPs[1], 1.0, 0., 0., w_fp[9] );
    // Amplitude(s) for diagram number 11
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[4], w_fp[9], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= 1. / 2. * amp_sv[0];
    jamp_sv[5] += 1. / 6. * amp_sv[0];

    // *** DIAGRAM 12 OF 72 ***
    // Wavefunction(s) for diagram number 12
    // (none)
    // Amplitude(s) for diagram number 12
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[12], w_fp[10], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= 1. / 6. * amp_sv[0];
    jamp_sv[5] += 1. / 2. * amp_sv[0];

    // *** DIAGRAM 13 OF 72 ***
    // Wavefunction(s) for diagram number 13
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[4], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[6] );
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[12], COUPs[1], 1.0, 0., 0., w_fp[14] );
    // Amplitude(s) for diagram number 13
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[6], w_fp[14], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= 1. / 6. * amp_sv[0];
    jamp_sv[5] += 1. / 2. * amp_sv[0];

    // *** DIAGRAM 14 OF 72 ***
    // Wavefunction(s) for diagram number 14
    // (none)
    // Amplitude(s) for diagram number 14
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[6], w_fp[9], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= 1. / 2. * amp_sv[0];
    jamp_sv[5] += 1. / 6. * amp_sv[0];

    // *** DIAGRAM 15 OF 72 ***
    // Wavefunction(s) for diagram number 15
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[5], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[15] );
    // Amplitude(s) for diagram number 15
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[4], w_fp[14], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += 1. / 2. * amp_sv[0];
    jamp_sv[4] -= 1. / 6. * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 12 );
    storeWf( wfs, w_cx, nevt, 13 );
    storeWf( wfs, w_cx, nevt, 14 );
    storeWf( wfs, w_cx, nevt, 15 );
#endif
  }

  //--------------------------------------------------------------------------

  __global__ void
  diagramgroup4( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 1 );
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 8 );
    retrieveWf( wfs, w_cx, nevt, 9 );
    retrieveWf( wfs, w_cx, nevt, 10 );
    retrieveWf( wfs, w_cx, nevt, 12 );
    retrieveWf( wfs, w_cx, nevt, 14 );
    retrieveWf( wfs, w_cx, nevt, 15 );
#endif

    // *** DIAGRAM 16 OF 72 ***
    // Wavefunction(s) for diagram number 16
    // (none)
    // Amplitude(s) for diagram number 16
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[12], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += 1. / 6. * amp_sv[0];
    jamp_sv[4] -= 1. / 2. * amp_sv[0];

    // *** DIAGRAM 17 OF 72 ***
    // Wavefunction(s) for diagram number 17
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[12], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[16] );
    // Amplitude(s) for diagram number 17
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[16], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= 1. / 2. * amp_sv[0];
    jamp_sv[1] += 1. / 6. * amp_sv[0];

    // *** DIAGRAM 18 OF 72 ***
    // Wavefunction(s) for diagram number 18
    // (none)
    // Amplitude(s) for diagram number 18
    VVV5_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[8], w_fp[9], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 19 OF 72 ***
    // Wavefunction(s) for diagram number 19
    // (none)
    // Amplitude(s) for diagram number 19
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[16], w_fp[10], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= 1. / 6. * amp_sv[0];
    jamp_sv[1] += 1. / 2. * amp_sv[0];

    // *** DIAGRAM 20 OF 72 ***
    // Wavefunction(s) for diagram number 20
    // (none)
    // Amplitude(s) for diagram number 20
    VVV5_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[10], w_fp[14], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE JAMPS ***
    for( int icol = 0; icol < ncolor; icol++ )
      J_ACCESS::kernelAccessIcol( jamps, icol ) += jamp_sv[icol]; // update jamps
#endif

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE WAVEFUNCTIONS FOR NEXT DIAGRAM GROUPS ***
    //for( int iwf = 0; iwf < nwf; iwf++ ) storeWf( wfs, w_cx, nevt, iwf );
    storeWf( wfs, w_cx, nevt, 16 );
#endif
  }

  //--------------------------------------------------------------------------

  __global__ void
  diagramgroup5( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 10 );
    retrieveWf( wfs, w_cx, nevt, 11 );
    retrieveWf( wfs, w_cx, nevt, 15 );
#endif

    // *** DIAGRAM 21 OF 72 ***
    // Wavefunction(s) for diagram number 21
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[14] );
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[16] );
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[14], w_fp[4], COUPs[1], 1.0, 0., 0., w_fp[9] );
    // Amplitude(s) for diagram number 21
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[16], w_fp[9], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += 1. / 6. * amp_sv[0];
    jamp_sv[8] -= 1. / 2. * amp_sv[0];

    // *** DIAGRAM 22 OF 72 ***
    // Wavefunction(s) for diagram number 22
    // (none)
    // Amplitude(s) for diagram number 22
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[14], w_fp[16], w_fp[10], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += 1. / 2. * amp_sv[0];
    jamp_sv[8] -= 1. / 6. * amp_sv[0];

    // *** DIAGRAM 23 OF 72 ***
    // Wavefunction(s) for diagram number 23
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[14], w_fp[2], COUPs[1], 1.0, 0., 0., w_fp[12] );
    // Amplitude(s) for diagram number 23
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[6], w_fp[12], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= 1. / 6. * amp_sv[0];
    jamp_sv[10] += 1. / 2. * amp_sv[0];

    // *** DIAGRAM 24 OF 72 ***
    // Wavefunction(s) for diagram number 24
    // (none)
    // Amplitude(s) for diagram number 24
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[14], w_fp[6], w_fp[11], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= 1. / 2. * amp_sv[0];
    jamp_sv[10] += 1. / 6. * amp_sv[0];

    // *** DIAGRAM 25 OF 72 ***
    // Wavefunction(s) for diagram number 25
    // (none)
    // Amplitude(s) for diagram number 25
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[4], w_fp[12], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= 1. / 6. * amp_sv[0];
    jamp_sv[6] += 1. / 2. * amp_sv[0];

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE JAMPS ***
    for( int icol = 0; icol < ncolor; icol++ )
      J_ACCESS::kernelAccessIcol( jamps, icol ) += jamp_sv[icol]; // update jamps
#endif

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE WAVEFUNCTIONS FOR NEXT DIAGRAM GROUPS ***
    //for( int iwf = 0; iwf < nwf; iwf++ ) storeWf( wfs, w_cx, nevt, iwf );
    storeWf( wfs, w_cx, nevt, 9 );
    storeWf( wfs, w_cx, nevt, 12 );
    storeWf( wfs, w_cx, nevt, 14 );
    storeWf( wfs, w_cx, nevt, 16 );
#endif
  }

  //--------------------------------------------------------------------------

  __global__ void
  diagramgroup6( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 1 );
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 9 );
    retrieveWf( wfs, w_cx, nevt, 10 );
    retrieveWf( wfs, w_cx, nevt, 11 );
    retrieveWf( wfs, w_cx, nevt, 12 );
    retrieveWf( wfs, w_cx, nevt, 14 );
    retrieveWf( wfs, w_cx, nevt, 15 );
#endif

    // *** DIAGRAM 26 OF 72 ***
    // Wavefunction(s) for diagram number 26
    // (none)
    // Amplitude(s) for diagram number 26
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[2], w_fp[9], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= 1. / 2. * amp_sv[0];
    jamp_sv[6] += 1. / 6. * amp_sv[0];

    // *** DIAGRAM 27 OF 72 ***
    // Wavefunction(s) for diagram number 27
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[14], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[17] );
    // Amplitude(s) for diagram number 27
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[17], w_fp[4], w_fp[11], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] -= 1. / 2. * amp_sv[0];
    jamp_sv[10] += 1. / 6. * amp_sv[0];

    // *** DIAGRAM 28 OF 72 ***
    // Wavefunction(s) for diagram number 28
    // (none)
    // Amplitude(s) for diagram number 28
    VVV5_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[11], w_fp[9], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 29 OF 72 ***
    // Wavefunction(s) for diagram number 29
    // (none)
    // Amplitude(s) for diagram number 29
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[17], w_fp[2], w_fp[10], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] -= 1. / 6. * amp_sv[0];
    jamp_sv[10] += 1. / 2. * amp_sv[0];

    // *** DIAGRAM 30 OF 72 ***
    // Wavefunction(s) for diagram number 30
    // (none)
    // Amplitude(s) for diagram number 30
    VVV5_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[10], w_fp[12], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE JAMPS ***
    for( int icol = 0; icol < ncolor; icol++ )
      J_ACCESS::kernelAccessIcol( jamps, icol ) += jamp_sv[icol]; // update jamps
#endif

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE WAVEFUNCTIONS FOR NEXT DIAGRAM GROUPS ***
    //for( int iwf = 0; iwf < nwf; iwf++ ) storeWf( wfs, w_cx, nevt, iwf );
    storeWf( wfs, w_cx, nevt, 17 );
#endif
  }

  //--------------------------------------------------------------------------

  __global__ void
  diagramgroup7( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 11 );
    retrieveWf( wfs, w_cx, nevt, 13 );
    retrieveWf( wfs, w_cx, nevt, 15 );
    retrieveWf( wfs, w_cx, nevt, 16 );
#endif

    // *** DIAGRAM 31 OF 72 ***
    // Wavefunction(s) for diagram number 31
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[4], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[12] );
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[12], COUPs[1], 1.0, 0., 0., w_fp[17] );
    // Amplitude(s) for diagram number 31
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[16], w_fp[17], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += 1. / 6. * amp_sv[0];
    jamp_sv[7] -= 1. / 2. * amp_sv[0];

    // *** DIAGRAM 32 OF 72 ***
    // Wavefunction(s) for diagram number 32
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[5], w_fp[12], COUPs[1], 1.0, 0., 0., w_fp[9] );
    // Amplitude(s) for diagram number 32
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[16], w_fp[9], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += 1. / 2. * amp_sv[0];
    jamp_sv[7] -= 1. / 6. * amp_sv[0];

    // *** DIAGRAM 33 OF 72 ***
    // Wavefunction(s) for diagram number 33
    // (none)
    // Amplitude(s) for diagram number 33
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[2], w_fp[9], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += 1. / 2. * amp_sv[0];
    jamp_sv[7] -= 1. / 6. * amp_sv[0];

    // *** DIAGRAM 34 OF 72 ***
    // Wavefunction(s) for diagram number 34
    // (none)
    // Amplitude(s) for diagram number 34
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[12], w_fp[11], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += 1. / 6. * amp_sv[0];
    jamp_sv[7] -= 1. / 2. * amp_sv[0];

    // *** DIAGRAM 35 OF 72 ***
    // Wavefunction(s) for diagram number 35
    // (none)
    // Amplitude(s) for diagram number 35
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[2], w_fp[17], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= 1. / 2. * amp_sv[0];
    jamp_sv[6] += 1. / 6. * amp_sv[0];

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE JAMPS ***
    for( int icol = 0; icol < ncolor; icol++ )
      J_ACCESS::kernelAccessIcol( jamps, icol ) += jamp_sv[icol]; // update jamps
#endif

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE WAVEFUNCTIONS FOR NEXT DIAGRAM GROUPS ***
    //for( int iwf = 0; iwf < nwf; iwf++ ) storeWf( wfs, w_cx, nevt, iwf );
    storeWf( wfs, w_cx, nevt, 9 );
    storeWf( wfs, w_cx, nevt, 12 );
    storeWf( wfs, w_cx, nevt, 17 );
#endif
  }

  //--------------------------------------------------------------------------

  __global__ void
  diagramgroup8( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 1 );
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 9 );
    retrieveWf( wfs, w_cx, nevt, 11 );
    retrieveWf( wfs, w_cx, nevt, 12 );
    retrieveWf( wfs, w_cx, nevt, 15 );
    retrieveWf( wfs, w_cx, nevt, 17 );
#endif

    // *** DIAGRAM 36 OF 72 ***
    // Wavefunction(s) for diagram number 36
    // (none)
    // Amplitude(s) for diagram number 36
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[12], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= 1. / 6. * amp_sv[0];
    jamp_sv[6] += 1. / 2. * amp_sv[0];

    // *** DIAGRAM 37 OF 72 ***
    // Wavefunction(s) for diagram number 37
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[12], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[14] );
    // Amplitude(s) for diagram number 37
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[14], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += 1. / 2. * amp_sv[0];
    jamp_sv[3] -= 1. / 6. * amp_sv[0];

    // *** DIAGRAM 38 OF 72 ***
    // Wavefunction(s) for diagram number 38
    // (none)
    // Amplitude(s) for diagram number 38
    VVV5_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[7], w_fp[9], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 39 OF 72 ***
    // Wavefunction(s) for diagram number 39
    // (none)
    // Amplitude(s) for diagram number 39
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[14], w_fp[11], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += 1. / 6. * amp_sv[0];
    jamp_sv[3] -= 1. / 2. * amp_sv[0];

    // *** DIAGRAM 40 OF 72 ***
    // Wavefunction(s) for diagram number 40
    // (none)
    // Amplitude(s) for diagram number 40
    VVV5_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[11], w_fp[17], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE JAMPS ***
    for( int icol = 0; icol < ncolor; icol++ )
      J_ACCESS::kernelAccessIcol( jamps, icol ) += jamp_sv[icol]; // update jamps
#endif

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE WAVEFUNCTIONS FOR NEXT DIAGRAM GROUPS ***
    //for( int iwf = 0; iwf < nwf; iwf++ ) storeWf( wfs, w_cx, nevt, iwf );
    storeWf( wfs, w_cx, nevt, 14 );
#endif
  }

  //--------------------------------------------------------------------------

  __global__ void
  diagramgroup9( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 8 );
    retrieveWf( wfs, w_cx, nevt, 13 );
    retrieveWf( wfs, w_cx, nevt, 16 );
#endif

    // *** DIAGRAM 41 OF 72 ***
    // Wavefunction(s) for diagram number 41
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[5], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[17] );
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[17], w_fp[4], COUPs[1], 1.0, 0., 0., w_fp[14] );
    // Amplitude(s) for diagram number 41
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[16], w_fp[14], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] -= 1. / 6. * amp_sv[0];
    jamp_sv[9] += 1. / 2. * amp_sv[0];

    // *** DIAGRAM 42 OF 72 ***
    // Wavefunction(s) for diagram number 42
    // (none)
    // Amplitude(s) for diagram number 42
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[17], w_fp[16], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] -= 1. / 2. * amp_sv[0];
    jamp_sv[9] += 1. / 6. * amp_sv[0];

    // *** DIAGRAM 43 OF 72 ***
    // Wavefunction(s) for diagram number 43
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[17], w_fp[2], COUPs[1], 1.0, 0., 0., w_fp[9] );
    // Amplitude(s) for diagram number 43
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[4], w_fp[9], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += 1. / 6. * amp_sv[0];
    jamp_sv[7] -= 1. / 2. * amp_sv[0];

    // *** DIAGRAM 44 OF 72 ***
    // Wavefunction(s) for diagram number 44
    // (none)
    // Amplitude(s) for diagram number 44
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[2], w_fp[14], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += 1. / 2. * amp_sv[0];
    jamp_sv[7] -= 1. / 6. * amp_sv[0];

    // *** DIAGRAM 45 OF 72 ***
    // Wavefunction(s) for diagram number 45
    // (none)
    // Amplitude(s) for diagram number 45
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[6], w_fp[9], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += 1. / 6. * amp_sv[0];
    jamp_sv[11] -= 1. / 2. * amp_sv[0];

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE JAMPS ***
    for( int icol = 0; icol < ncolor; icol++ )
      J_ACCESS::kernelAccessIcol( jamps, icol ) += jamp_sv[icol]; // update jamps
#endif

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE WAVEFUNCTIONS FOR NEXT DIAGRAM GROUPS ***
    //for( int iwf = 0; iwf < nwf; iwf++ ) storeWf( wfs, w_cx, nevt, iwf );
    storeWf( wfs, w_cx, nevt, 9 );
    storeWf( wfs, w_cx, nevt, 14 );
    storeWf( wfs, w_cx, nevt, 17 );
#endif
  }

  //--------------------------------------------------------------------------

  __global__ void
  diagramgroup10( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 1 );
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 8 );
    retrieveWf( wfs, w_cx, nevt, 9 );
    retrieveWf( wfs, w_cx, nevt, 14 );
    retrieveWf( wfs, w_cx, nevt, 17 );
#endif

    // *** DIAGRAM 46 OF 72 ***
    // Wavefunction(s) for diagram number 46
    // (none)
    // Amplitude(s) for diagram number 46
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[17], w_fp[6], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += 1. / 2. * amp_sv[0];
    jamp_sv[11] -= 1. / 6. * amp_sv[0];

    // *** DIAGRAM 47 OF 72 ***
    // Wavefunction(s) for diagram number 47
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[17], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[12] );
    // Amplitude(s) for diagram number 47
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[4], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += 1. / 2. * amp_sv[0];
    jamp_sv[11] -= 1. / 6. * amp_sv[0];

    // *** DIAGRAM 48 OF 72 ***
    // Wavefunction(s) for diagram number 48
    // (none)
    // Amplitude(s) for diagram number 48
    VVV5_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[7], w_fp[14], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 49 OF 72 ***
    // Wavefunction(s) for diagram number 49
    // (none)
    // Amplitude(s) for diagram number 49
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[2], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += 1. / 6. * amp_sv[0];
    jamp_sv[11] -= 1. / 2. * amp_sv[0];

    // *** DIAGRAM 50 OF 72 ***
    // Wavefunction(s) for diagram number 50
    // (none)
    // Amplitude(s) for diagram number 50
    VVV5_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[8], w_fp[9], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE JAMPS ***
    for( int icol = 0; icol < ncolor; icol++ )
      J_ACCESS::kernelAccessIcol( jamps, icol ) += jamp_sv[icol]; // update jamps
#endif

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE WAVEFUNCTIONS FOR NEXT DIAGRAM GROUPS ***
    //for( int iwf = 0; iwf < nwf; iwf++ ) storeWf( wfs, w_cx, nevt, iwf );
    storeWf( wfs, w_cx, nevt, 12 );
#endif
  }

  //--------------------------------------------------------------------------

  __global__ void
  diagramgroup11( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 8 );
    retrieveWf( wfs, w_cx, nevt, 10 );
    retrieveWf( wfs, w_cx, nevt, 11 );
    retrieveWf( wfs, w_cx, nevt, 13 );
    retrieveWf( wfs, w_cx, nevt, 16 );
#endif

    // *** DIAGRAM 51 OF 72 ***
    // Wavefunction(s) for diagram number 51
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[16], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[9] );
    // Amplitude(s) for diagram number 51
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[9], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] -= 1. / 2. * amp_sv[0];
    jamp_sv[9] += 1. / 6. * amp_sv[0];

    // *** DIAGRAM 52 OF 72 ***
    // Wavefunction(s) for diagram number 52
    VVV5P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[8], COUPs[0], 1.0, 0., 0., w_fp[12] );
    // Amplitude(s) for diagram number 52
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[16], w_fp[12], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 53 OF 72 ***
    // Wavefunction(s) for diagram number 53
    // (none)
    // Amplitude(s) for diagram number 53
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[9], w_fp[10], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] -= 1. / 6. * amp_sv[0];
    jamp_sv[9] += 1. / 2. * amp_sv[0];

    // *** DIAGRAM 54 OF 72 ***
    // Wavefunction(s) for diagram number 54
    VVV5P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[10], COUPs[0], 1.0, 0., 0., w_fp[9] );
    // Amplitude(s) for diagram number 54
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[16], w_fp[9], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 55 OF 72 ***
    // Wavefunction(s) for diagram number 55
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[13], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[16] );
    // Amplitude(s) for diagram number 55
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[4], w_fp[11], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= 1. / 2. * amp_sv[0];
    jamp_sv[2] += 1. / 6. * amp_sv[0];

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE JAMPS ***
    for( int icol = 0; icol < ncolor; icol++ )
      J_ACCESS::kernelAccessIcol( jamps, icol ) += jamp_sv[icol]; // update jamps
#endif

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE WAVEFUNCTIONS FOR NEXT DIAGRAM GROUPS ***
    //for( int iwf = 0; iwf < nwf; iwf++ ) storeWf( wfs, w_cx, nevt, iwf );
    storeWf( wfs, w_cx, nevt, 9 );
    storeWf( wfs, w_cx, nevt, 12 );
    storeWf( wfs, w_cx, nevt, 16 );
#endif
  }

  //--------------------------------------------------------------------------

  __global__ void
  diagramgroup12( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 9 );
    retrieveWf( wfs, w_cx, nevt, 10 );
    retrieveWf( wfs, w_cx, nevt, 11 );
    retrieveWf( wfs, w_cx, nevt, 13 );
    retrieveWf( wfs, w_cx, nevt, 16 );
#endif

    // *** DIAGRAM 56 OF 72 ***
    // Wavefunction(s) for diagram number 56
    VVV5P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[11], COUPs[0], 1.0, 0., 0., w_fp[14] );
    // Amplitude(s) for diagram number 56
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[4], w_fp[14], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 57 OF 72 ***
    // Wavefunction(s) for diagram number 57
    // (none)
    // Amplitude(s) for diagram number 57
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[2], w_fp[10], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= 1. / 6. * amp_sv[0];
    jamp_sv[2] += 1. / 2. * amp_sv[0];

    // *** DIAGRAM 58 OF 72 ***
    // Wavefunction(s) for diagram number 58
    // (none)
    // Amplitude(s) for diagram number 58
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[2], w_fp[9], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 59 OF 72 ***
    // Wavefunction(s) for diagram number 59
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[6], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[13] );
    // Amplitude(s) for diagram number 59
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[13], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += 1. / 2. * amp_sv[0];
    jamp_sv[11] -= 1. / 6. * amp_sv[0];

    // *** DIAGRAM 60 OF 72 ***
    // Wavefunction(s) for diagram number 60
    VVV5P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[7], COUPs[0], 1.0, 0., 0., w_fp[16] );
    // Amplitude(s) for diagram number 60
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[6], w_fp[16], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE JAMPS ***
    for( int icol = 0; icol < ncolor; icol++ )
      J_ACCESS::kernelAccessIcol( jamps, icol ) += jamp_sv[icol]; // update jamps
#endif

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE WAVEFUNCTIONS FOR NEXT DIAGRAM GROUPS ***
    //for( int iwf = 0; iwf < nwf; iwf++ ) storeWf( wfs, w_cx, nevt, iwf );
    storeWf( wfs, w_cx, nevt, 13 );
    storeWf( wfs, w_cx, nevt, 14 );
    storeWf( wfs, w_cx, nevt, 16 );
#endif
  }

  //--------------------------------------------------------------------------

  __global__ void
  diagramgroup13( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 8 );
    retrieveWf( wfs, w_cx, nevt, 11 );
    retrieveWf( wfs, w_cx, nevt, 13 );
    retrieveWf( wfs, w_cx, nevt, 14 );
    retrieveWf( wfs, w_cx, nevt, 15 );
    retrieveWf( wfs, w_cx, nevt, 16 );
#endif

    // *** DIAGRAM 61 OF 72 ***
    // Wavefunction(s) for diagram number 61
    // (none)
    // Amplitude(s) for diagram number 61
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[13], w_fp[11], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += 1. / 6. * amp_sv[0];
    jamp_sv[11] -= 1. / 2. * amp_sv[0];

    // *** DIAGRAM 62 OF 72 ***
    // Wavefunction(s) for diagram number 62
    // (none)
    // Amplitude(s) for diagram number 62
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[6], w_fp[14], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 63 OF 72 ***
    // Wavefunction(s) for diagram number 63
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[15], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[6] );
    // Amplitude(s) for diagram number 63
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[6], w_fp[4], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += 1. / 2. * amp_sv[0];
    jamp_sv[3] -= 1. / 6. * amp_sv[0];

    // *** DIAGRAM 64 OF 72 ***
    // Wavefunction(s) for diagram number 64
    // (none)
    // Amplitude(s) for diagram number 64
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[4], w_fp[16], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 65 OF 72 ***
    // Wavefunction(s) for diagram number 65
    // (none)
    // Amplitude(s) for diagram number 65
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[6], w_fp[2], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += 1. / 6. * amp_sv[0];
    jamp_sv[3] -= 1. / 2. * amp_sv[0];

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE JAMPS ***
    for( int icol = 0; icol < ncolor; icol++ )
      J_ACCESS::kernelAccessIcol( jamps, icol ) += jamp_sv[icol]; // update jamps
#endif

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE WAVEFUNCTIONS FOR NEXT DIAGRAM GROUPS ***
    //for( int iwf = 0; iwf < nwf; iwf++ ) storeWf( wfs, w_cx, nevt, iwf );
    storeWf( wfs, w_cx, nevt, 6 );
#endif
  }

  //--------------------------------------------------------------------------

  __global__ void
  diagramgroup14( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 8 );
    retrieveWf( wfs, w_cx, nevt, 9 );
    retrieveWf( wfs, w_cx, nevt, 10 );
    retrieveWf( wfs, w_cx, nevt, 11 );
    retrieveWf( wfs, w_cx, nevt, 12 );
    retrieveWf( wfs, w_cx, nevt, 15 );
    retrieveWf( wfs, w_cx, nevt, 16 );
#endif

    // *** DIAGRAM 66 OF 72 ***
    // Wavefunction(s) for diagram number 66
    // (none)
    // Amplitude(s) for diagram number 66
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[2], w_fp[12], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 67 OF 72 ***
    // Wavefunction(s) for diagram number 67
    // (none)
    // Amplitude(s) for diagram number 67
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[7], w_fp[10], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += 1. / 2. * amp_sv[0];
    jamp_sv[2] -= 1. / 2. * amp_sv[0];
    jamp_sv[9] -= 1. / 2. * amp_sv[0];
    jamp_sv[10] += 1. / 2. * amp_sv[0];
    VVVV9_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[7], w_fp[10], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += 1. / 2. * amp_sv[0];
    jamp_sv[5] -= 1. / 2. * amp_sv[0];
    jamp_sv[6] -= 1. / 2. * amp_sv[0];
    jamp_sv[10] += 1. / 2. * amp_sv[0];
    VVVV10_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[7], w_fp[10], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += 1. / 2. * amp_sv[0];
    jamp_sv[5] -= 1. / 2. * amp_sv[0];
    jamp_sv[6] -= 1. / 2. * amp_sv[0];
    jamp_sv[9] += 1. / 2. * amp_sv[0];

    // *** DIAGRAM 68 OF 72 ***
    // Wavefunction(s) for diagram number 68
    // (none)
    // Amplitude(s) for diagram number 68
    VVV5_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[10], w_fp[16], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += 1. / 2. * amp_sv[0];
    jamp_sv[5] -= 1. / 2. * amp_sv[0];
    jamp_sv[6] -= 1. / 2. * amp_sv[0];
    jamp_sv[10] += 1. / 2. * amp_sv[0];

    // *** DIAGRAM 69 OF 72 ***
    // Wavefunction(s) for diagram number 69
    // (none)
    // Amplitude(s) for diagram number 69
    VVV5_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[7], w_fp[9], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += 1. / 2. * amp_sv[0];
    jamp_sv[5] -= 1. / 2. * amp_sv[0];
    jamp_sv[6] -= 1. / 2. * amp_sv[0];
    jamp_sv[9] += 1. / 2. * amp_sv[0];

    // *** DIAGRAM 70 OF 72 ***
    // Wavefunction(s) for diagram number 70
    // (none)
    // Amplitude(s) for diagram number 70
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[11], w_fp[8], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= 1. / 2. * amp_sv[0];
    jamp_sv[3] += 1. / 2. * amp_sv[0];
    jamp_sv[8] += 1. / 2. * amp_sv[0];
    jamp_sv[11] -= 1. / 2. * amp_sv[0];
    VVVV9_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[11], w_fp[8], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= 1. / 2. * amp_sv[0];
    jamp_sv[4] += 1. / 2. * amp_sv[0];
    jamp_sv[7] += 1. / 2. * amp_sv[0];
    jamp_sv[11] -= 1. / 2. * amp_sv[0];
    VVVV10_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[11], w_fp[8], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= 1. / 2. * amp_sv[0];
    jamp_sv[4] += 1. / 2. * amp_sv[0];
    jamp_sv[7] += 1. / 2. * amp_sv[0];
    jamp_sv[8] -= 1. / 2. * amp_sv[0];

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE JAMPS ***
    for( int icol = 0; icol < ncolor; icol++ )
      J_ACCESS::kernelAccessIcol( jamps, icol ) += jamp_sv[icol]; // update jamps
#endif

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE WAVEFUNCTIONS FOR NEXT DIAGRAM GROUPS ***
    //for( int iwf = 0; iwf < nwf; iwf++ ) storeWf( wfs, w_cx, nevt, iwf );
    // (none)
#endif
  }

  //--------------------------------------------------------------------------

  __global__ void
  diagramgroup15( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 1 );
    retrieveWf( wfs, w_cx, nevt, 8 );
    retrieveWf( wfs, w_cx, nevt, 11 );
    retrieveWf( wfs, w_cx, nevt, 12 );
    retrieveWf( wfs, w_cx, nevt, 14 );
#endif

    // *** DIAGRAM 71 OF 72 ***
    // Wavefunction(s) for diagram number 71
    // (none)
    // Amplitude(s) for diagram number 71
    VVV5_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[8], w_fp[14], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= 1. / 2. * amp_sv[0];
    jamp_sv[4] += 1. / 2. * amp_sv[0];
    jamp_sv[7] += 1. / 2. * amp_sv[0];
    jamp_sv[11] -= 1. / 2. * amp_sv[0];

    // *** DIAGRAM 72 OF 72 ***
    // Wavefunction(s) for diagram number 72
    // (none)
    // Amplitude(s) for diagram number 72
    VVV5_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[11], w_fp[12], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= 1. / 2. * amp_sv[0];
    jamp_sv[4] += 1. / 2. * amp_sv[0];
    jamp_sv[7] += 1. / 2. * amp_sv[0];
    jamp_sv[8] -= 1. / 2. * amp_sv[0];

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
