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
    constexpr int npar = CPPProcess::npar;
    const short (*cHel)[npar] = reinterpret_cast<const short(*)[npar]>( cHelFlat );

    // *** DIAGRAM 1 OF 36 ***
    // Wavefunction(s) for diagram number 1
    vxxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][0], -1, w_fp[0], 0 );
    oxxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][1], -1, w_fp[1], 1 );
    oxxxxx<M_ACCESS, W_ACCESS>( momenta, cIPD[0], cHel[ihel][2], +1, w_fp[2], 2 );
    ixxxxx<M_ACCESS, W_ACCESS>( momenta, cIPD[0], cHel[ihel][3], -1, w_fp[3], 3 );
    vxxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][4], +1, w_fp[4], 4 );
    ixxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][5], -1, w_fp[5], 5 );
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[4], COUPs[0], 1.0, 0., 0., w_fp[6] );
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[2], COUPs[1], 1.0, 0., 0., w_fp[7] );
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[6], COUPs[1], 1.0, 0., 0., w_fp[8] );
    // Amplitude(s) for diagram number 1
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[8], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 1 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[4] -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] += 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] -= 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2 OF 36 ***
    // Wavefunction(s) for diagram number 2
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[5], w_fp[6], COUPs[1], 1.0, 0., 0., w_fp[8] );
    // Amplitude(s) for diagram number 2
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[1], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 2 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[5] += 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] -= 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3 OF 36 ***
    // Wavefunction(s) for diagram number 3
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[5], w_fp[1], COUPs[1], 1.0, 0., 0., w_fp[8] );
    // Amplitude(s) for diagram number 3
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[6], w_fp[7], w_fp[8], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 3 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[4] -= 1. / 2. * amp_sv[0];
    jamp_sv[7] += 1. / 2. * amp_sv[0];
    jamp_sv[8] -= 1. / 2. * amp_sv[0];
    jamp_sv[11] += 1. / 2. * amp_sv[0];

    // *** DIAGRAM 4 OF 36 ***
    // Wavefunction(s) for diagram number 4
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[6], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[9] );
    // Amplitude(s) for diagram number 4
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[9], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 4 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[6] += 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] -= 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5 OF 36 ***
    // Wavefunction(s) for diagram number 5
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[6], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[9] );
    // Amplitude(s) for diagram number 5
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[2], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 5 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[4] -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] += 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] -= 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6 OF 36 ***
    // Wavefunction(s) for diagram number 6
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[9] );
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[4], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[6] );
    // Amplitude(s) for diagram number 6
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[6], w_fp[9], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 6 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[3] += 1. / 2. * amp_sv[0];
    jamp_sv[6] -= 1. / 6. * amp_sv[0];

    // *** DIAGRAM 7 OF 36 ***
    // Wavefunction(s) for diagram number 7
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[4], COUPs[1], 1.0, 0., 0., w_fp[10] );
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[9], COUPs[1], 1.0, 0., 0., w_fp[11] );
    // Amplitude(s) for diagram number 7
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[10], w_fp[11], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 7 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[2] -= 1. / 6. * amp_sv[0];
    jamp_sv[3] += 1. / 2. * amp_sv[0];

    // *** DIAGRAM 8 OF 36 ***
    // Wavefunction(s) for diagram number 8
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[5], w_fp[4], COUPs[1], 1.0, 0., 0., w_fp[12] );
    // Amplitude(s) for diagram number 8
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[1], w_fp[11], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 8 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[2] -= 1. / 6. * amp_sv[0];
    jamp_sv[7] += 1. / 2. * amp_sv[0];

    // *** DIAGRAM 9 OF 36 ***
    // Wavefunction(s) for diagram number 9
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[9], w_fp[4], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[13] );
    // Amplitude(s) for diagram number 9
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[13], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 9 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[6] -= 1. / 6. * amp_sv[0];
    jamp_sv[7] += 1. / 2. * amp_sv[0];

    // *** DIAGRAM 10 OF 36 ***
    // Wavefunction(s) for diagram number 10
    // (none)
    // Amplitude(s) for diagram number 10
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[4], w_fp[8], w_fp[11], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 10 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[3] += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11 OF 36 ***
    // Wavefunction(s) for diagram number 11
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[11] );
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[4], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[13] );
    // Amplitude(s) for diagram number 11
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[11], w_fp[13], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 11 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += 1. / 2. * amp_sv[0];
    jamp_sv[9] -= 1. / 6. * amp_sv[0];

    // *** DIAGRAM 12 OF 36 ***
    // Wavefunction(s) for diagram number 12
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[11], w_fp[2], COUPs[1], 1.0, 0., 0., w_fp[9] );
    // Amplitude(s) for diagram number 12
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[10], w_fp[9], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 12 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[2] -= 1. / 6. * amp_sv[0];
    jamp_sv[11] += 1. / 2. * amp_sv[0];

    // *** DIAGRAM 13 OF 36 ***
    // Wavefunction(s) for diagram number 13
    // (none)
    // Amplitude(s) for diagram number 13
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[1], w_fp[9], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 13 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += 1. / 2. * amp_sv[0];
    jamp_sv[2] -= 1. / 6. * amp_sv[0];

    // *** DIAGRAM 14 OF 36 ***
    // Wavefunction(s) for diagram number 14
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[11], w_fp[4], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[14] );
    // Amplitude(s) for diagram number 14
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[14], w_fp[2], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 14 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[9] -= 1. / 6. * amp_sv[0];
    jamp_sv[11] += 1. / 2. * amp_sv[0];

    // *** DIAGRAM 15 OF 36 ***
    // Wavefunction(s) for diagram number 15
    // (none)
    // Amplitude(s) for diagram number 15
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[4], w_fp[8], w_fp[9], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 15 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 16 OF 36 ***
    // Wavefunction(s) for diagram number 16
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[0], COUPs[1], 1.0, 0., 0., w_fp[9] );
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[5], w_fp[9], COUPs[1], 1.0, 0., 0., w_fp[14] );
    // Amplitude(s) for diagram number 16
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[13], w_fp[14], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 16 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += 1. / 2. * amp_sv[0];
    jamp_sv[1] -= 1. / 6. * amp_sv[0];

    // *** DIAGRAM 17 OF 36 ***
    // Wavefunction(s) for diagram number 17
    // (none)
    // Amplitude(s) for diagram number 17
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[6], w_fp[2], w_fp[14], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 17 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[1] -= 1. / 6. * amp_sv[0];
    jamp_sv[4] += 1. / 2. * amp_sv[0];

    // *** DIAGRAM 18 OF 36 ***
    // Wavefunction(s) for diagram number 18
    // (none)
    // Amplitude(s) for diagram number 18
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[9], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 18 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += 1. / 2. * amp_sv[0];
    jamp_sv[5] -= 1. / 6. * amp_sv[0];

    // *** DIAGRAM 19 OF 36 ***
    // Wavefunction(s) for diagram number 19
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[9], w_fp[4], COUPs[1], 1.0, 0., 0., w_fp[11] );
    // Amplitude(s) for diagram number 19
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[11], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 19 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[4] += 1. / 2. * amp_sv[0];
    jamp_sv[5] -= 1. / 6. * amp_sv[0];

    // *** DIAGRAM 20 OF 36 ***
    // Wavefunction(s) for diagram number 20
    // (none)
    // Amplitude(s) for diagram number 20
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[4], w_fp[7], w_fp[14], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 20 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 21 OF 36 ***
    // Wavefunction(s) for diagram number 21
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[5], w_fp[0], COUPs[1], 1.0, 0., 0., w_fp[14] );
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[14], w_fp[1], COUPs[1], 1.0, 0., 0., w_fp[11] );
    // Amplitude(s) for diagram number 21
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[13], w_fp[11], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 21 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[1] -= 1. / 6. * amp_sv[0];
    jamp_sv[8] += 1. / 2. * amp_sv[0];

    // *** DIAGRAM 22 OF 36 ***
    // Wavefunction(s) for diagram number 22
    // (none)
    // Amplitude(s) for diagram number 22
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[6], w_fp[2], w_fp[11], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 22 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[1] -= 1. / 6. * amp_sv[0];
    jamp_sv[3] += 1. / 2. * amp_sv[0];

    // *** DIAGRAM 23 OF 36 ***
    // Wavefunction(s) for diagram number 23
    // (none)
    // Amplitude(s) for diagram number 23
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[14], w_fp[10], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 23 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[3] += 1. / 2. * amp_sv[0];
    jamp_sv[10] -= 1. / 6. * amp_sv[0];

    // *** DIAGRAM 24 OF 36 ***
    // Wavefunction(s) for diagram number 24
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[14], w_fp[4], COUPs[1], 1.0, 0., 0., w_fp[9] );
    // Amplitude(s) for diagram number 24
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[1], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 24 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[8] += 1. / 2. * amp_sv[0];
    jamp_sv[10] -= 1. / 6. * amp_sv[0];

    // *** DIAGRAM 25 OF 36 ***
    // Wavefunction(s) for diagram number 25
    // (none)
    // Amplitude(s) for diagram number 25
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[4], w_fp[7], w_fp[11], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 25 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[3] -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 26 OF 36 ***
    // Wavefunction(s) for diagram number 26
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[13], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[11] );
    // Amplitude(s) for diagram number 26
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[11], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 26 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[8] += 1. / 2. * amp_sv[0];
    jamp_sv[9] -= 1. / 6. * amp_sv[0];

    // *** DIAGRAM 27 OF 36 ***
    // Wavefunction(s) for diagram number 27
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[8], COUPs[0], 1.0, 0., 0., w_fp[11] );
    // Amplitude(s) for diagram number 27
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[13], w_fp[11], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 27 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 28 OF 36 ***
    // Wavefunction(s) for diagram number 28
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[6], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[13] );
    // Amplitude(s) for diagram number 28
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[2], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 28 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[4] += 1. / 2. * amp_sv[0];
    jamp_sv[6] -= 1. / 6. * amp_sv[0];

    // *** DIAGRAM 29 OF 36 ***
    // Wavefunction(s) for diagram number 29
    // (none)
    // Amplitude(s) for diagram number 29
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[6], w_fp[2], w_fp[11], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 29 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[3] -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 30 OF 36 ***
    // Wavefunction(s) for diagram number 30
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[10], w_fp[0], COUPs[1], 1.0, 0., 0., w_fp[6] );
    // Amplitude(s) for diagram number 30
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[6], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 30 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[10] -= 1. / 6. * amp_sv[0];
    jamp_sv[11] += 1. / 2. * amp_sv[0];

    // *** DIAGRAM 31 OF 36 ***
    // Wavefunction(s) for diagram number 31
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[7], COUPs[0], 1.0, 0., 0., w_fp[6] );
    // Amplitude(s) for diagram number 31
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[10], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 31 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[3] += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 32 OF 36 ***
    // Wavefunction(s) for diagram number 32
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[12], w_fp[0], COUPs[1], 1.0, 0., 0., w_fp[10] );
    // Amplitude(s) for diagram number 32
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[10], w_fp[1], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 32 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[5] -= 1. / 6. * amp_sv[0];
    jamp_sv[7] += 1. / 2. * amp_sv[0];

    // *** DIAGRAM 33 OF 36 ***
    // Wavefunction(s) for diagram number 33
    // (none)
    // Amplitude(s) for diagram number 33
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[1], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 33 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 34 OF 36 ***
    // Wavefunction(s) for diagram number 34
    // (none)
    // Amplitude(s) for diagram number 34
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[4], w_fp[7], w_fp[8], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[4] -= 1. / 2. * amp_sv[0];
    jamp_sv[7] += 1. / 2. * amp_sv[0];
    jamp_sv[8] -= 1. / 2. * amp_sv[0];
    jamp_sv[11] += 1. / 2. * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[4], w_fp[7], w_fp[8], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[0] -= 1. / 2. * amp_sv[0];
    jamp_sv[3] -= 1. / 2. * amp_sv[0];
    jamp_sv[7] += 1. / 2. * amp_sv[0];
    jamp_sv[11] += 1. / 2. * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[4], w_fp[7], w_fp[8], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[0] -= 1. / 2. * amp_sv[0];
    jamp_sv[3] -= 1. / 2. * amp_sv[0];
    jamp_sv[4] += 1. / 2. * amp_sv[0];
    jamp_sv[8] += 1. / 2. * amp_sv[0];

    // *** DIAGRAM 35 OF 36 ***
    // Wavefunction(s) for diagram number 35
    // (none)
    // Amplitude(s) for diagram number 35
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[4], w_fp[8], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 35 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] -= 1. / 2. * amp_sv[0];
    jamp_sv[3] -= 1. / 2. * amp_sv[0];
    jamp_sv[7] += 1. / 2. * amp_sv[0];
    jamp_sv[11] += 1. / 2. * amp_sv[0];

    // *** DIAGRAM 36 OF 36 ***
    // Wavefunction(s) for diagram number 36
    // (none)
    // Amplitude(s) for diagram number 36
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[4], w_fp[7], w_fp[11], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 36 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] -= 1. / 2. * amp_sv[0];
    jamp_sv[3] -= 1. / 2. * amp_sv[0];
    jamp_sv[4] += 1. / 2. * amp_sv[0];
    jamp_sv[8] += 1. / 2. * amp_sv[0];

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE JAMPS ***
    for( int icol = 0; icol < ncolor; icol++ )
      J_ACCESS::kernelAccessIcol( jamps, icol ) = jamp_sv[icol]; // set jamps
#endif

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE WAVEFUNCTIONS FOR NEXT DIAGRAM GROUPS ***
    // (none)
#endif
  }

  //--------------------------------------------------------------------------
}
