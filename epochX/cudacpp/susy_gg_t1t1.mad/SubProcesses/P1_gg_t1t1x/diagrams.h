// Copyright (C) 2020-2025 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Sep 2025) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2025) for the MG5aMC CUDACPP plugin.

/* clang-format off */

  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram1( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // A uniform interface for diagramXXX including channelIDs, numerators and denominators is used also #ifndef MGONGPU_SUPPORTS_MULTICHANNEL
    // In that case, however, the boilerplate code asserts that all three pointers all nullptr as a sanity check
#include "diagram_boilerplate.h"
#ifdef MGONGPUCPP_GPUIMPL
    using M_ACCESS = DeviceAccessMomenta; // non-trivial access: buffer includes all events
#else
    using M_ACCESS = HostAccessMomenta; // non-trivial access: buffer includes all events
#endif
    // *** DIAGRAM 1 OF 6 ***
    // Wavefunction(s) for diagram number 1
    vxxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][0], -1, w_fp[0], 0 );
    vxxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][1], -1, w_fp[1], 1 );
    sxxxxx<M_ACCESS, W_ACCESS>( momenta, +1, w_fp[2], 2 );
    sxxxxx<M_ACCESS, W_ACCESS>( momenta, +1, w_fp[3], 3 );
    // Amplitude(s) for diagram number 1
    VVSS1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[3], w_fp[2], COUPs[0], 1.0, &amp_fp[0] );
    J_ACCESS::kernelAccessIcol( jamps, 1 ) += amp_sv[0];
    VVSS1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[3], w_fp[2], COUPs[0], 1.0, &amp_fp[0] );
    J_ACCESS::kernelAccessIcol( jamps, 0 ) += amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram2( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
            fptype* jamps,                  // output jamps[ncolor*2*nevtORneppV]
            const unsigned int* channelIds, // input: channelIds[nevt] for GPU or SCALAR channelId[0] for C++ (1 to #diagrams, 0 to disable SDE)
#ifdef MGONGPUCPP_GPUIMPL
            const fptype* couplings,        // input: dependent couplings[nevt*ndcoup*2] for all events
#else
            const fptype** COUPs,           // input: dependent and independent COUPs[nxcoup] for this event page
#endif
            fptype* numerators,             // input/output: multichannel numerators[nevtORneppV], add helicity ihel
            fptype* denominators )          // input/output: multichannel denominators[nevtORneppV], add helicity ihel
  {
    // A uniform interface for diagramXXX including channelIDs, numerators and denominators is used also #ifndef MGONGPU_SUPPORTS_MULTICHANNEL
    // In that case, however, the boilerplate code asserts that all three pointers all nullptr as a sanity check
#include "diagram_boilerplate.h"
    // *** DIAGRAM 2 OF 6 ***
    // Wavefunction(s) for diagram number 2
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], COUPs[1], 1.0, 0., 0., w_fp[4] );
    // Amplitude(s) for diagram number 2
    VSS1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[4], w_fp[3], w_fp[2], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 2 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 1 ) += cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram3( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
            fptype* jamps,                  // output jamps[ncolor*2*nevtORneppV]
            const unsigned int* channelIds, // input: channelIds[nevt] for GPU or SCALAR channelId[0] for C++ (1 to #diagrams, 0 to disable SDE)
#ifdef MGONGPUCPP_GPUIMPL
            const fptype* couplings,        // input: dependent couplings[nevt*ndcoup*2] for all events
#else
            const fptype** COUPs,           // input: dependent and independent COUPs[nxcoup] for this event page
#endif
            fptype* numerators,             // input/output: multichannel numerators[nevtORneppV], add helicity ihel
            fptype* denominators )          // input/output: multichannel denominators[nevtORneppV], add helicity ihel
  {
    // A uniform interface for diagramXXX including channelIDs, numerators and denominators is used also #ifndef MGONGPU_SUPPORTS_MULTICHANNEL
    // In that case, however, the boilerplate code asserts that all three pointers all nullptr as a sanity check
#include "diagram_boilerplate.h"
    // *** DIAGRAM 3 OF 6 ***
    // Wavefunction(s) for diagram number 3
    VSS1_2<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[2], COUPs[2], 1.0, cIPD[0], cIPD[1], w_fp[4] );
    // Amplitude(s) for diagram number 3
    VSS1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[3], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 3 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) += amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram4( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
            fptype* jamps,                  // output jamps[ncolor*2*nevtORneppV]
            const unsigned int* channelIds, // input: channelIds[nevt] for GPU or SCALAR channelId[0] for C++ (1 to #diagrams, 0 to disable SDE)
#ifdef MGONGPUCPP_GPUIMPL
            const fptype* couplings,        // input: dependent couplings[nevt*ndcoup*2] for all events
#else
            const fptype** COUPs,           // input: dependent and independent COUPs[nxcoup] for this event page
#endif
            fptype* numerators,             // input/output: multichannel numerators[nevtORneppV], add helicity ihel
            fptype* denominators )          // input/output: multichannel denominators[nevtORneppV], add helicity ihel
  {
    // A uniform interface for diagramXXX including channelIDs, numerators and denominators is used also #ifndef MGONGPU_SUPPORTS_MULTICHANNEL
    // In that case, however, the boilerplate code asserts that all three pointers all nullptr as a sanity check
#include "diagram_boilerplate.h"
    // *** DIAGRAM 4 OF 6 ***
    // Wavefunction(s) for diagram number 4
    VSS1_3<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[2], COUPs[3], -1.0, cIPD[2], cIPD[3], w_fp[4] );
    // Amplitude(s) for diagram number 4
    VSS1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[3], w_fp[4], COUPs[3], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 4 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) += amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram5( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
            fptype* jamps,                  // output jamps[ncolor*2*nevtORneppV]
            const unsigned int* channelIds, // input: channelIds[nevt] for GPU or SCALAR channelId[0] for C++ (1 to #diagrams, 0 to disable SDE)
#ifdef MGONGPUCPP_GPUIMPL
            const fptype* couplings,        // input: dependent couplings[nevt*ndcoup*2] for all events
#else
            const fptype** COUPs,           // input: dependent and independent COUPs[nxcoup] for this event page
#endif
            fptype* numerators,             // input/output: multichannel numerators[nevtORneppV], add helicity ihel
            fptype* denominators )          // input/output: multichannel denominators[nevtORneppV], add helicity ihel
  {
    // A uniform interface for diagramXXX including channelIDs, numerators and denominators is used also #ifndef MGONGPU_SUPPORTS_MULTICHANNEL
    // In that case, however, the boilerplate code asserts that all three pointers all nullptr as a sanity check
#include "diagram_boilerplate.h"
    // *** DIAGRAM 5 OF 6 ***
    // Wavefunction(s) for diagram number 5
    VSS1_3<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[3], COUPs[2], 1.0, cIPD[0], cIPD[1], w_fp[4] );
    // Amplitude(s) for diagram number 5
    VSS1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[4], w_fp[2], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 5 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    J_ACCESS::kernelAccessIcol( jamps, 1 ) += amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram6( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
            fptype* jamps,                  // output jamps[ncolor*2*nevtORneppV]
            const unsigned int* channelIds, // input: channelIds[nevt] for GPU or SCALAR channelId[0] for C++ (1 to #diagrams, 0 to disable SDE)
#ifdef MGONGPUCPP_GPUIMPL
            const fptype* couplings,        // input: dependent couplings[nevt*ndcoup*2] for all events
#else
            const fptype** COUPs,           // input: dependent and independent COUPs[nxcoup] for this event page
#endif
            fptype* numerators,             // input/output: multichannel numerators[nevtORneppV], add helicity ihel
            fptype* denominators )          // input/output: multichannel denominators[nevtORneppV], add helicity ihel
  {
    // A uniform interface for diagramXXX including channelIDs, numerators and denominators is used also #ifndef MGONGPU_SUPPORTS_MULTICHANNEL
    // In that case, however, the boilerplate code asserts that all three pointers all nullptr as a sanity check
#include "diagram_boilerplate.h"
    // *** DIAGRAM 6 OF 6 ***
    // Wavefunction(s) for diagram number 6
    VSS1_3<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[3], COUPs[3], 1.0, cIPD[2], cIPD[3], w_fp[4] );
    // Amplitude(s) for diagram number 6
    VSS1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[2], w_fp[4], COUPs[3], -1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 6 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    J_ACCESS::kernelAccessIcol( jamps, 1 ) += amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

/* clang-format on */
