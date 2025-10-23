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

    // *** DIAGRAM 1 OF 123 ***
    // Wavefunction(s) for diagram number 1
    vxxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][0], -1, w_fp[0], 0 );
    vxxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][1], -1, w_fp[1], 1 );
    oxxxxx<M_ACCESS, W_ACCESS>( momenta, cIPD[0], cHel[ihel][2], +1, w_fp[2], 2 );
    ixxxxx<M_ACCESS, W_ACCESS>( momenta, cIPD[0], cHel[ihel][3], -1, w_fp[3], 3 );
    vxxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][4], +1, w_fp[4], 4 );
    vxxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][5], +1, w_fp[5], 5 );
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], COUPs[0], 1.0, 0., 0., w_fp[6] );
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[2], COUPs[1], 1.0, 0., 0., w_fp[7] );
    // Amplitude(s) for diagram number 1
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[6], w_fp[7], w_fp[4], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[6], w_fp[7], w_fp[4], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[6], w_fp[7], w_fp[4], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2 OF 123 ***
    // Wavefunction(s) for diagram number 2
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[6], w_fp[4], COUPs[0], 1.0, 0., 0., w_fp[8] );
    // Amplitude(s) for diagram number 2
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[7], w_fp[5], w_fp[8], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 2 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3 OF 123 ***
    // Wavefunction(s) for diagram number 3
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[6], w_fp[5], COUPs[0], 1.0, 0., 0., w_fp[9] );
    // Amplitude(s) for diagram number 3
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[7], w_fp[4], w_fp[9], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 3 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4 OF 123 ***
    // Wavefunction(s) for diagram number 4
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[4], w_fp[5], COUPs[0], 1.0, 0., 0., w_fp[10] );
    // Amplitude(s) for diagram number 4
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[6], w_fp[7], w_fp[10], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 4 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5 OF 123 ***
    // Wavefunction(s) for diagram number 5
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[4], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[11] );
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[6], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[12] );
    // Amplitude(s) for diagram number 5
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[11], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 5 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6 OF 123 ***
    // Wavefunction(s) for diagram number 6
    // (none)
    // Amplitude(s) for diagram number 6
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[11], w_fp[9], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 6 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[12] += amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[17] += amp_sv[0];

    // *** DIAGRAM 7 OF 123 ***
    // Wavefunction(s) for diagram number 7
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[5], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[13] );
    // Amplitude(s) for diagram number 7
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[11], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 7 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8 OF 123 ***
    // Wavefunction(s) for diagram number 8
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[5], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[14] );
    // Amplitude(s) for diagram number 8
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[14], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 8 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9 OF 123 ***
    // Wavefunction(s) for diagram number 9
    // (none)
    // Amplitude(s) for diagram number 9
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[14], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 9 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[18] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];

    // *** DIAGRAM 10 OF 123 ***
    // Wavefunction(s) for diagram number 10
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[4], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[15] );
    // Amplitude(s) for diagram number 10
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[14], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 10 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11 OF 123 ***
    // Wavefunction(s) for diagram number 11
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[6], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[16] );
    // Amplitude(s) for diagram number 11
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[16], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 11 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12 OF 123 ***
    // Wavefunction(s) for diagram number 12
    // (none)
    // Amplitude(s) for diagram number 12
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[2], w_fp[9], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 12 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];

    // *** DIAGRAM 13 OF 123 ***
    // Wavefunction(s) for diagram number 13
    // (none)
    // Amplitude(s) for diagram number 13
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[16], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 13 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14 OF 123 ***
    // Wavefunction(s) for diagram number 14
    // (none)
    // Amplitude(s) for diagram number 14
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[2], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 14 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];

    // *** DIAGRAM 15 OF 123 ***
    // Wavefunction(s) for diagram number 15
    // (none)
    // Amplitude(s) for diagram number 15
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[16], w_fp[10], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 15 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[7] += amp_sv[0];

    // *** DIAGRAM 16 OF 123 ***
    // Wavefunction(s) for diagram number 16
    // (none)
    // Amplitude(s) for diagram number 16
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[2], w_fp[10], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 16 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[16] += amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];

    // *** DIAGRAM 17 OF 123 ***
    // Wavefunction(s) for diagram number 17
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[12] );
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[16] );
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[12], w_fp[4], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[8] );
    // Amplitude(s) for diagram number 17
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[8], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 17 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[3] -= amp_sv[0];

    // *** DIAGRAM 18 OF 123 ***
    // Wavefunction(s) for diagram number 18
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[12], w_fp[5], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[9] );
    // Amplitude(s) for diagram number 18
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[9], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 18 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[5] -= amp_sv[0];

    // *** DIAGRAM 19 OF 123 ***
    // Wavefunction(s) for diagram number 19
    // (none)
    // Amplitude(s) for diagram number 19
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[12], w_fp[10], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 19 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 20 OF 123 ***
    // Wavefunction(s) for diagram number 20
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[4], COUPs[0], 1.0, 0., 0., w_fp[6] );
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[12], COUPs[1], 1.0, 0., 0., w_fp[17] );
    // Amplitude(s) for diagram number 20
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[6], w_fp[5], w_fp[17], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 20 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];

    // *** DIAGRAM 21 OF 123 ***
    // Wavefunction(s) for diagram number 21
    // (none)
    // Amplitude(s) for diagram number 21
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[9], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 21 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 22 OF 123 ***
    // Wavefunction(s) for diagram number 22
    // (none)
    // Amplitude(s) for diagram number 22
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[12], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 22 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 23 OF 123 ***
    // Wavefunction(s) for diagram number 23
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[5], COUPs[0], 1.0, 0., 0., w_fp[18] );
    // Amplitude(s) for diagram number 23
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[18], w_fp[4], w_fp[17], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 23 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[3] += amp_sv[0];
    jamp_sv[4] -= amp_sv[0];

    // *** DIAGRAM 24 OF 123 ***
    // Wavefunction(s) for diagram number 24
    // (none)
    // Amplitude(s) for diagram number 24
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[8], w_fp[18], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 24 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 25 OF 123 ***
    // Wavefunction(s) for diagram number 25
    // (none)
    // Amplitude(s) for diagram number 25
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[12], w_fp[18], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 25 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 26 OF 123 ***
    // Wavefunction(s) for diagram number 26
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[12], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[19] );
    // Amplitude(s) for diagram number 26
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[19], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 26 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[1] -= amp_sv[0];

    // *** DIAGRAM 27 OF 123 ***
    // Wavefunction(s) for diagram number 27
    // (none)
    // Amplitude(s) for diagram number 27
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[9], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 27 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[4] -= amp_sv[0];

    // *** DIAGRAM 28 OF 123 ***
    // Wavefunction(s) for diagram number 28
    // (none)
    // Amplitude(s) for diagram number 28
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[19], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 28 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] -= amp_sv[0];

    // *** DIAGRAM 29 OF 123 ***
    // Wavefunction(s) for diagram number 29
    // (none)
    // Amplitude(s) for diagram number 29
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[8], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 29 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[2] -= amp_sv[0];

    // *** DIAGRAM 30 OF 123 ***
    // Wavefunction(s) for diagram number 30
    // (none)
    // Amplitude(s) for diagram number 30
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[19], w_fp[10], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 30 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 31 OF 123 ***
    // Wavefunction(s) for diagram number 31
    // (none)
    // Amplitude(s) for diagram number 31
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[10], w_fp[17], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 31 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];

    // *** DIAGRAM 32 OF 123 ***
    // Wavefunction(s) for diagram number 32
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[4], w_fp[5], COUPs[2], 1.0, 0., 0., w_fp[17] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[4], w_fp[5], COUPs[2], 1.0, 0., 0., w_fp[19] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[4], w_fp[5], COUPs[2], 1.0, 0., 0., w_fp[8] );
    // Amplitude(s) for diagram number 32
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[12], w_fp[17], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[0] += amp_sv[0];
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[12], w_fp[19], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[2] += amp_sv[0];
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[4] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[12], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[2] += amp_sv[0];
    jamp_sv[4] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];

    // *** DIAGRAM 33 OF 123 ***
    // Wavefunction(s) for diagram number 33
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[12] );
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[9] );
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[12], w_fp[4], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[20] );
    // Amplitude(s) for diagram number 33
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[20], w_fp[9], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 33 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[11] -= amp_sv[0];

    // *** DIAGRAM 34 OF 123 ***
    // Wavefunction(s) for diagram number 34
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[12], w_fp[5], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[21] );
    // Amplitude(s) for diagram number 34
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[9], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 34 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[9] -= amp_sv[0];

    // *** DIAGRAM 35 OF 123 ***
    // Wavefunction(s) for diagram number 35
    // (none)
    // Amplitude(s) for diagram number 35
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[9], w_fp[10], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 35 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 36 OF 123 ***
    // Wavefunction(s) for diagram number 36
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[12], w_fp[2], COUPs[1], 1.0, 0., 0., w_fp[22] );
    // Amplitude(s) for diagram number 36
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[6], w_fp[5], w_fp[22], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 36 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];

    // *** DIAGRAM 37 OF 123 ***
    // Wavefunction(s) for diagram number 37
    // (none)
    // Amplitude(s) for diagram number 37
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[2], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 37 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 38 OF 123 ***
    // Wavefunction(s) for diagram number 38
    // (none)
    // Amplitude(s) for diagram number 38
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[14], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 38 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 39 OF 123 ***
    // Wavefunction(s) for diagram number 39
    // (none)
    // Amplitude(s) for diagram number 39
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[18], w_fp[4], w_fp[22], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 39 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[11] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];

    // *** DIAGRAM 40 OF 123 ***
    // Wavefunction(s) for diagram number 40
    // (none)
    // Amplitude(s) for diagram number 40
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[20], w_fp[2], w_fp[18], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 40 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 41 OF 123 ***
    // Wavefunction(s) for diagram number 41
    // (none)
    // Amplitude(s) for diagram number 41
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[11], w_fp[18], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 41 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 42 OF 123 ***
    // Wavefunction(s) for diagram number 42
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[12], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[23] );
    // Amplitude(s) for diagram number 42
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[23], w_fp[11], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 42 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[17] -= amp_sv[0];

    // *** DIAGRAM 43 OF 123 ***
    // Wavefunction(s) for diagram number 43
    // (none)
    // Amplitude(s) for diagram number 43
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[11], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 43 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[15] -= amp_sv[0];

    // *** DIAGRAM 44 OF 123 ***
    // Wavefunction(s) for diagram number 44
    // (none)
    // Amplitude(s) for diagram number 44
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[23], w_fp[14], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 44 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[23] -= amp_sv[0];

    // *** DIAGRAM 45 OF 123 ***
    // Wavefunction(s) for diagram number 45
    // (none)
    // Amplitude(s) for diagram number 45
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[20], w_fp[14], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 45 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[21] -= amp_sv[0];

    // *** DIAGRAM 46 OF 123 ***
    // Wavefunction(s) for diagram number 46
    // (none)
    // Amplitude(s) for diagram number 46
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[23], w_fp[2], w_fp[10], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 46 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 47 OF 123 ***
    // Wavefunction(s) for diagram number 47
    // (none)
    // Amplitude(s) for diagram number 47
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[10], w_fp[22], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 47 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];

    // *** DIAGRAM 48 OF 123 ***
    // Wavefunction(s) for diagram number 48
    // (none)
    // Amplitude(s) for diagram number 48
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[2], w_fp[17], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[9] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[2], w_fp[19], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[15] += amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[2], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[15] += amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];

    // *** DIAGRAM 49 OF 123 ***
    // Wavefunction(s) for diagram number 49
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[4], COUPs[0], 1.0, 0., 0., w_fp[12] );
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[12], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[22] );
    // Amplitude(s) for diagram number 49
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[22], w_fp[9], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 49 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 50 OF 123 ***
    // Wavefunction(s) for diagram number 50
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[12], w_fp[5], COUPs[0], 1.0, 0., 0., w_fp[23] );
    // Amplitude(s) for diagram number 50
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[9], w_fp[23], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 50 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[6] += amp_sv[0];
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];

    // *** DIAGRAM 51 OF 123 ***
    // Wavefunction(s) for diagram number 51
    // (none)
    // Amplitude(s) for diagram number 51
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[9], w_fp[12], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 51 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 52 OF 123 ***
    // Wavefunction(s) for diagram number 52
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[12], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[20] );
    // Amplitude(s) for diagram number 52
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[20], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 52 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 53 OF 123 ***
    // Wavefunction(s) for diagram number 53
    // (none)
    // Amplitude(s) for diagram number 53
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[2], w_fp[23], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 53 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];

    // *** DIAGRAM 54 OF 123 ***
    // Wavefunction(s) for diagram number 54
    // (none)
    // Amplitude(s) for diagram number 54
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[14], w_fp[12], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 54 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 55 OF 123 ***
    // Wavefunction(s) for diagram number 55
    // (none)
    // Amplitude(s) for diagram number 55
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[20], w_fp[18], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 55 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[13] += amp_sv[0];

    // *** DIAGRAM 56 OF 123 ***
    // Wavefunction(s) for diagram number 56
    // (none)
    // Amplitude(s) for diagram number 56
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[22], w_fp[2], w_fp[18], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 56 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[21] += amp_sv[0];

    // *** DIAGRAM 57 OF 123 ***
    // Wavefunction(s) for diagram number 57
    // (none)
    // Amplitude(s) for diagram number 57
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[18], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 57 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 58 OF 123 ***
    // Wavefunction(s) for diagram number 58
    // (none)
    // Amplitude(s) for diagram number 58
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[1], w_fp[7], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[1], w_fp[7], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[1], w_fp[7], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 59 OF 123 ***
    // Wavefunction(s) for diagram number 59
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[12], w_fp[1], COUPs[0], 1.0, 0., 0., w_fp[21] );
    // Amplitude(s) for diagram number 59
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[7], w_fp[5], w_fp[21], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 59 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 60 OF 123 ***
    // Wavefunction(s) for diagram number 60
    // (none)
    // Amplitude(s) for diagram number 60
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[7], w_fp[23], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 60 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 61 OF 123 ***
    // Wavefunction(s) for diagram number 61
    // (none)
    // Amplitude(s) for diagram number 61
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[14], w_fp[21], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 61 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[19] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];

    // *** DIAGRAM 62 OF 123 ***
    // Wavefunction(s) for diagram number 62
    // (none)
    // Amplitude(s) for diagram number 62
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[22], w_fp[14], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 62 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 63 OF 123 ***
    // Wavefunction(s) for diagram number 63
    // (none)
    // Amplitude(s) for diagram number 63
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[2], w_fp[21], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 63 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[12] -= amp_sv[0];

    // *** DIAGRAM 64 OF 123 ***
    // Wavefunction(s) for diagram number 64
    // (none)
    // Amplitude(s) for diagram number 64
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[20], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 64 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 65 OF 123 ***
    // Wavefunction(s) for diagram number 65
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[5], COUPs[0], 1.0, 0., 0., w_fp[20] );
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[20], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[21] );
    // Amplitude(s) for diagram number 65
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[9], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 65 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 66 OF 123 ***
    // Wavefunction(s) for diagram number 66
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[20], w_fp[4], COUPs[0], 1.0, 0., 0., w_fp[22] );
    // Amplitude(s) for diagram number 66
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[9], w_fp[22], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 66 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[7] += amp_sv[0];
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[9] += amp_sv[0];
    jamp_sv[10] -= amp_sv[0];

    // *** DIAGRAM 67 OF 123 ***
    // Wavefunction(s) for diagram number 67
    // (none)
    // Amplitude(s) for diagram number 67
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[9], w_fp[20], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 67 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 68 OF 123 ***
    // Wavefunction(s) for diagram number 68
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[20], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[23] );
    // Amplitude(s) for diagram number 68
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[23], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 68 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 69 OF 123 ***
    // Wavefunction(s) for diagram number 69
    // (none)
    // Amplitude(s) for diagram number 69
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[2], w_fp[22], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 69 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];

    // *** DIAGRAM 70 OF 123 ***
    // Wavefunction(s) for diagram number 70
    // (none)
    // Amplitude(s) for diagram number 70
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[11], w_fp[20], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 70 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 71 OF 123 ***
    // Wavefunction(s) for diagram number 71
    // (none)
    // Amplitude(s) for diagram number 71
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[23], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 71 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[19] += amp_sv[0];

    // *** DIAGRAM 72 OF 123 ***
    // Wavefunction(s) for diagram number 72
    // (none)
    // Amplitude(s) for diagram number 72
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[2], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 72 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[8] += amp_sv[0];
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[15] += amp_sv[0];

    // *** DIAGRAM 73 OF 123 ***
    // Wavefunction(s) for diagram number 73
    // (none)
    // Amplitude(s) for diagram number 73
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[20], w_fp[6], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 73 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 74 OF 123 ***
    // Wavefunction(s) for diagram number 74
    // (none)
    // Amplitude(s) for diagram number 74
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[20], w_fp[1], w_fp[7], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[20], w_fp[1], w_fp[7], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[20], w_fp[1], w_fp[7], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 75 OF 123 ***
    // Wavefunction(s) for diagram number 75
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[20], w_fp[1], COUPs[0], 1.0, 0., 0., w_fp[12] );
    // Amplitude(s) for diagram number 75
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[7], w_fp[4], w_fp[12], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 75 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 76 OF 123 ***
    // Wavefunction(s) for diagram number 76
    // (none)
    // Amplitude(s) for diagram number 76
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[7], w_fp[22], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 76 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 77 OF 123 ***
    // Wavefunction(s) for diagram number 77
    // (none)
    // Amplitude(s) for diagram number 77
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[11], w_fp[12], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 77 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[13] += amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[15] += amp_sv[0];
    jamp_sv[16] -= amp_sv[0];

    // *** DIAGRAM 78 OF 123 ***
    // Wavefunction(s) for diagram number 78
    // (none)
    // Amplitude(s) for diagram number 78
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[11], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 78 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 79 OF 123 ***
    // Wavefunction(s) for diagram number 79
    // (none)
    // Amplitude(s) for diagram number 79
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[2], w_fp[12], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 79 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[18] -= amp_sv[0];

    // *** DIAGRAM 80 OF 123 ***
    // Wavefunction(s) for diagram number 80
    // (none)
    // Amplitude(s) for diagram number 80
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[23], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 80 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 81 OF 123 ***
    // Wavefunction(s) for diagram number 81
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[9], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[23] );
    // Amplitude(s) for diagram number 81
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[23], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 81 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[7] -= amp_sv[0];

    // *** DIAGRAM 82 OF 123 ***
    // Wavefunction(s) for diagram number 82
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[15], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[12] );
    // Amplitude(s) for diagram number 82
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[9], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 82 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[10] -= amp_sv[0];

    // *** DIAGRAM 83 OF 123 ***
    // Wavefunction(s) for diagram number 83
    // (none)
    // Amplitude(s) for diagram number 83
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[23], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 83 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[6] -= amp_sv[0];

    // *** DIAGRAM 84 OF 123 ***
    // Wavefunction(s) for diagram number 84
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[13], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[21] );
    // Amplitude(s) for diagram number 84
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[9], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 84 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[8] -= amp_sv[0];

    // *** DIAGRAM 85 OF 123 ***
    // Wavefunction(s) for diagram number 85
    // (none)
    // Amplitude(s) for diagram number 85
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[23], w_fp[10], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 85 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 86 OF 123 ***
    // Wavefunction(s) for diagram number 86
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[10], COUPs[0], 1.0, 0., 0., w_fp[23] );
    // Amplitude(s) for diagram number 86
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[9], w_fp[23], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 86 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[6] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];

    // *** DIAGRAM 87 OF 123 ***
    // Wavefunction(s) for diagram number 87
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[16], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[22] );
    // Amplitude(s) for diagram number 87
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[22], w_fp[11], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 87 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[16] -= amp_sv[0];

    // *** DIAGRAM 88 OF 123 ***
    // Wavefunction(s) for diagram number 88
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[11], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[20] );
    // Amplitude(s) for diagram number 88
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[20], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 88 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[13] -= amp_sv[0];

    // *** DIAGRAM 89 OF 123 ***
    // Wavefunction(s) for diagram number 89
    // (none)
    // Amplitude(s) for diagram number 89
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[22], w_fp[14], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 89 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[22] -= amp_sv[0];

    // *** DIAGRAM 90 OF 123 ***
    // Wavefunction(s) for diagram number 90
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[14], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[24] );
    // Amplitude(s) for diagram number 90
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[24], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 90 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[19] -= amp_sv[0];

    // *** DIAGRAM 91 OF 123 ***
    // Wavefunction(s) for diagram number 91
    // (none)
    // Amplitude(s) for diagram number 91
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[22], w_fp[2], w_fp[10], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 91 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 92 OF 123 ***
    // Wavefunction(s) for diagram number 92
    // (none)
    // Amplitude(s) for diagram number 92
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[2], w_fp[23], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 92 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];

    // *** DIAGRAM 93 OF 123 ***
    // Wavefunction(s) for diagram number 93
    // (none)
    // Amplitude(s) for diagram number 93
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[6], w_fp[7], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[6], w_fp[7], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[6], w_fp[7], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 94 OF 123 ***
    // Wavefunction(s) for diagram number 94
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[6], COUPs[0], 1.0, 0., 0., w_fp[22] );
    // Amplitude(s) for diagram number 94
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[7], w_fp[5], w_fp[22], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 94 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 95 OF 123 ***
    // Wavefunction(s) for diagram number 95
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[7], COUPs[0], 1.0, 0., 0., w_fp[25] );
    // Amplitude(s) for diagram number 95
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[6], w_fp[5], w_fp[25], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 95 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 96 OF 123 ***
    // Wavefunction(s) for diagram number 96
    // (none)
    // Amplitude(s) for diagram number 96
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[14], w_fp[22], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 96 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[18] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];

    // *** DIAGRAM 97 OF 123 ***
    // Wavefunction(s) for diagram number 97
    // (none)
    // Amplitude(s) for diagram number 97
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[24], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 97 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 98 OF 123 ***
    // Wavefunction(s) for diagram number 98
    // (none)
    // Amplitude(s) for diagram number 98
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[2], w_fp[22], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 98 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];

    // *** DIAGRAM 99 OF 123 ***
    // Wavefunction(s) for diagram number 99
    // (none)
    // Amplitude(s) for diagram number 99
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[2], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 99 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 100 OF 123 ***
    // Wavefunction(s) for diagram number 100
    // (none)
    // Amplitude(s) for diagram number 100
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[18], w_fp[7], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[18], w_fp[7], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[18], w_fp[7], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 13 );
    storeWf( wfs, w_cx, nevt, 14 );
    storeWf( wfs, w_cx, nevt, 15 );
    storeWf( wfs, w_cx, nevt, 16 );
    storeWf( wfs, w_cx, nevt, 17 );
    storeWf( wfs, w_cx, nevt, 18 );
    storeWf( wfs, w_cx, nevt, 19 );
    storeWf( wfs, w_cx, nevt, 20 );
    storeWf( wfs, w_cx, nevt, 21 );
    storeWf( wfs, w_cx, nevt, 22 );
    storeWf( wfs, w_cx, nevt, 23 );
    storeWf( wfs, w_cx, nevt, 24 );
    storeWf( wfs, w_cx, nevt, 25 );
#endif
  }

  //--------------------------------------------------------------------------
}
