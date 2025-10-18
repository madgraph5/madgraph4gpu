// Copyright (C) 2020-2025 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Sep 2025) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2025) for the MG5aMC CUDACPP plugin.

/* clang-format off */

  //--------------------------------------------------------------------------

  __global__ void
  diagramgroup1( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
                 fptype* jamps,                  // output jamps[ncolor*2*nevtORneppV]
                 const unsigned int* channelIds, // input: channelIds[nevt] for GPU or SCALAR channelId[0] for C++ (1 to #diagrams, 0 to disable SDE)
#ifdef MGONGPUCPP_GPUIMPL
                 const fptype* couplings,        // input: dependent couplings[nevt*ndcoup*2] for all events
#else
                 const fptype** COUPs,           // input: dependent and independent COUPs[nxcoup] for this event page
#endif
                 fptype* numerators,             // input/output: multichannel numerators[nevtORneppV], add helicity ihel
                 fptype* denominators,           // input/output: multichannel denominators[nevtORneppV], add helicity ihel
                 const fptype* momenta,          // input: momenta[npar*4*nevtORneppV]
                 const int ihel )                // input: helicity (0 to ncomb)
  {
    // A uniform interface for diagramgroupXXX including channelIDs, numerators and denominators is used also #ifndef MGONGPU_SUPPORTS_MULTICHANNEL
    // In that case, however, the boilerplate code asserts that all three pointers all nullptr as a sanity check
#include "diagram_boilerplate.h"
#ifdef MGONGPUCPP_GPUIMPL
    using M_ACCESS = DeviceAccessMomenta; // non-trivial access: buffer includes all events
#else
    using M_ACCESS = HostAccessMomenta; // non-trivial access: buffer includes all events
#endif

#ifdef MGONGPUCPP_GPUIMPL
    // *** RETRIEVE WAVEFUNCTIONS FROM PREVIOUS DIAGRAM GROUPS ***
    // (none)
#endif

    // *** DIAGRAM 1 OF 5 ***
    // Wavefunction(s) for diagram number 1
    ixxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][0], +1, w_fp[0], 0 );
    oxxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][1], -1, w_fp[1], 1 );
    oxxxxx<M_ACCESS, W_ACCESS>( momenta, cIPD[0], cHel[ihel][2], +1, w_fp[2], 2 );
    ixxxxx<M_ACCESS, W_ACCESS>( momenta, cIPD[0], cHel[ihel][3], -1, w_fp[3], 3 );
    vxxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][4], +1, w_fp[4], 4 );
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[4], COUPs[1], 1.0, 0., 0., w_fp[5] );
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[2], COUPs[1], 1.0, 0., 0., w_fp[6] );
    // Amplitude(s) for diagram number 1
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[1], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 1 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    J_ACCESS::kernelAccessIcol( jamps, 1 ) -= 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 3 ) += 1. / 6. * amp_sv[0];

    // *** DIAGRAM 2 OF 5 ***
    // Wavefunction(s) for diagram number 2
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[4], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[5] );
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], COUPs[1], 1.0, 0., 0., w_fp[7] );
    // Amplitude(s) for diagram number 2
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[5], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 2 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) += 1. / 6. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 1 ) -= 1. / 2. * amp_sv[0];

    // *** DIAGRAM 3 OF 5 ***
    // Wavefunction(s) for diagram number 3
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[4], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[5] );
    // Amplitude(s) for diagram number 3
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[2], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 3 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) += 1. / 6. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 2 ) -= 1. / 2. * amp_sv[0];

    // *** DIAGRAM 4 OF 5 ***
    // Wavefunction(s) for diagram number 4
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[4], COUPs[1], 1.0, 0., 0., w_fp[5] );
    // Amplitude(s) for diagram number 4
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[5], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 4 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    J_ACCESS::kernelAccessIcol( jamps, 2 ) -= 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 3 ) += 1. / 6. * amp_sv[0];

    // *** DIAGRAM 5 OF 5 ***
    // Wavefunction(s) for diagram number 5
    // (none)
    // Amplitude(s) for diagram number 5
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[4], w_fp[7], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 5 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    J_ACCESS::kernelAccessIcol( jamps, 1 ) += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 2 ) -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE WAVEFUNCTIONS FOR NEXT DIAGRAM GROUPS ***
    // (none)
#endif
  }

  //--------------------------------------------------------------------------

/* clang-format on */
