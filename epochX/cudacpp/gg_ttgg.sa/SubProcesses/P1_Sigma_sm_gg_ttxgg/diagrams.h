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
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 1 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 6 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 7 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 16 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 17 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 22 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 23 ) -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[6], w_fp[7], w_fp[4], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 6 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 12 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 14 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 18 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 20 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 22 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 23 ) -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[6], w_fp[7], w_fp[4], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 1 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 7 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 12 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 14 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 16 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 17 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 18 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 20 ) += cxtype( 0, 1 ) * amp_sv[0];
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
    // *** DIAGRAM 2 OF 123 ***
    // Wavefunction(s) for diagram number 2
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[6], w_fp[4], COUPs[0], 1.0, 0., 0., w_fp[8] );
    // Amplitude(s) for diagram number 2
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[7], w_fp[5], w_fp[8], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 6 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 12 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 14 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 18 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 20 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 22 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 23 ) -= cxtype( 0, 1 ) * amp_sv[0];
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
    // *** DIAGRAM 3 OF 123 ***
    // Wavefunction(s) for diagram number 3
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[6], w_fp[5], COUPs[0], 1.0, 0., 0., w_fp[9] );
    // Amplitude(s) for diagram number 3
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[7], w_fp[4], w_fp[9], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 1 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 7 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 12 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 14 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 16 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 17 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 18 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 20 ) += cxtype( 0, 1 ) * amp_sv[0];
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
    // *** DIAGRAM 4 OF 123 ***
    // Wavefunction(s) for diagram number 4
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[4], w_fp[5], COUPs[0], 1.0, 0., 0., w_fp[10] );
    // Amplitude(s) for diagram number 4
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[6], w_fp[7], w_fp[10], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 1 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 6 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 7 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 16 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 17 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 22 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 23 ) -= cxtype( 0, 1 ) * amp_sv[0];
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
    // *** DIAGRAM 5 OF 123 ***
    // Wavefunction(s) for diagram number 5
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[4], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[11] );
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[6], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[12] );
    // Amplitude(s) for diagram number 5
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[11], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 16 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 17 ) -= cxtype( 0, 1 ) * amp_sv[0];
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
    // *** DIAGRAM 6 OF 123 ***
    // Wavefunction(s) for diagram number 6
    // (none)
    // Amplitude(s) for diagram number 6
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[11], w_fp[9], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 12 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 14 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 16 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 17 ) += amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram7( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 7 OF 123 ***
    // Wavefunction(s) for diagram number 7
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[5], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[13] );
    // Amplitude(s) for diagram number 7
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[11], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 12 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 14 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram8( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 8 OF 123 ***
    // Wavefunction(s) for diagram number 8
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[5], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[14] );
    // Amplitude(s) for diagram number 8
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[14], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 22 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 23 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram9( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 9 OF 123 ***
    // Wavefunction(s) for diagram number 9
    // (none)
    // Amplitude(s) for diagram number 9
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[14], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 18 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 20 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 22 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 23 ) += amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram10( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 10 OF 123 ***
    // Wavefunction(s) for diagram number 10
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[4], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[15] );
    // Amplitude(s) for diagram number 10
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[14], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 18 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 20 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram11( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 11 OF 123 ***
    // Wavefunction(s) for diagram number 11
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[6], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[16] );
    // Amplitude(s) for diagram number 11
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[16], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 1 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 7 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram12( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 12 OF 123 ***
    // Wavefunction(s) for diagram number 12
    // (none)
    // Amplitude(s) for diagram number 12
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[2], w_fp[9], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 1 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 7 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 18 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 20 ) += amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram13( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 13 OF 123 ***
    // Wavefunction(s) for diagram number 13
    // (none)
    // Amplitude(s) for diagram number 13
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[16], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 6 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram14( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 14 OF 123 ***
    // Wavefunction(s) for diagram number 14
    // (none)
    // Amplitude(s) for diagram number 14
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[2], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 6 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 12 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 14 ) += amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram15( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 15 OF 123 ***
    // Wavefunction(s) for diagram number 15
    // (none)
    // Amplitude(s) for diagram number 15
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[16], w_fp[10], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 1 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 6 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 7 ) += amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram16( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 16 OF 123 ***
    // Wavefunction(s) for diagram number 16
    // (none)
    // Amplitude(s) for diagram number 16
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[2], w_fp[10], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 16 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 17 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 22 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 23 ) += amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram17( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 17 OF 123 ***
    // Wavefunction(s) for diagram number 17
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[12] );
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[16] );
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[12], w_fp[4], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[8] );
    // Amplitude(s) for diagram number 17
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[8], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 3 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram18( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 18 OF 123 ***
    // Wavefunction(s) for diagram number 18
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[12], w_fp[5], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[9] );
    // Amplitude(s) for diagram number 18
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[9], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 5 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram19( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 19 OF 123 ***
    // Wavefunction(s) for diagram number 19
    // (none)
    // Amplitude(s) for diagram number 19
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[12], w_fp[10], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 3 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 5 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram20( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 20 OF 123 ***
    // Wavefunction(s) for diagram number 20
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[4], COUPs[0], 1.0, 0., 0., w_fp[6] );
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[12], COUPs[1], 1.0, 0., 0., w_fp[17] );
    // Amplitude(s) for diagram number 20
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[6], w_fp[5], w_fp[17], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 2 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 4 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 5 ) += amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram21( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 21 OF 123 ***
    // Wavefunction(s) for diagram number 21
    // (none)
    // Amplitude(s) for diagram number 21
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[9], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 4 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 5 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram22( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 22 OF 123 ***
    // Wavefunction(s) for diagram number 22
    // (none)
    // Amplitude(s) for diagram number 22
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[12], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 2 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram23( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 23 OF 123 ***
    // Wavefunction(s) for diagram number 23
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[5], COUPs[0], 1.0, 0., 0., w_fp[18] );
    // Amplitude(s) for diagram number 23
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[18], w_fp[4], w_fp[17], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 1 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 2 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 3 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 4 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram24( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 24 OF 123 ***
    // Wavefunction(s) for diagram number 24
    // (none)
    // Amplitude(s) for diagram number 24
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[8], w_fp[18], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 2 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 3 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram25( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 25 OF 123 ***
    // Wavefunction(s) for diagram number 25
    // (none)
    // Amplitude(s) for diagram number 25
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[12], w_fp[18], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 1 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 4 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram26( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 26 OF 123 ***
    // Wavefunction(s) for diagram number 26
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[12], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[19] );
    // Amplitude(s) for diagram number 26
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[19], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 1 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram27( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 27 OF 123 ***
    // Wavefunction(s) for diagram number 27
    // (none)
    // Amplitude(s) for diagram number 27
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[9], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 4 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram28( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 28 OF 123 ***
    // Wavefunction(s) for diagram number 28
    // (none)
    // Amplitude(s) for diagram number 28
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[19], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram29( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 29 OF 123 ***
    // Wavefunction(s) for diagram number 29
    // (none)
    // Amplitude(s) for diagram number 29
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[8], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 2 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram30( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 30 OF 123 ***
    // Wavefunction(s) for diagram number 30
    // (none)
    // Amplitude(s) for diagram number 30
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[19], w_fp[10], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 1 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram31( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 31 OF 123 ***
    // Wavefunction(s) for diagram number 31
    // (none)
    // Amplitude(s) for diagram number 31
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[10], w_fp[17], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 1 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 3 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 5 ) += amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram32( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 32 OF 123 ***
    // Wavefunction(s) for diagram number 32
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[4], w_fp[5], COUPs[2], 1.0, 0., 0., w_fp[17] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[4], w_fp[5], COUPs[2], 1.0, 0., 0., w_fp[19] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[4], w_fp[5], COUPs[2], 1.0, 0., 0., w_fp[8] );
    // Amplitude(s) for diagram number 32
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[12], w_fp[17], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 1 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 3 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 5 ) += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[12], w_fp[19], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 1 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 2 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 3 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 4 ) += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[12], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 2 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 4 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 5 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram33( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 33 OF 123 ***
    // Wavefunction(s) for diagram number 33
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[12] );
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[9] );
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[12], w_fp[4], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[20] );
    // Amplitude(s) for diagram number 33
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[20], w_fp[9], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 11 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram34( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 34 OF 123 ***
    // Wavefunction(s) for diagram number 34
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[12], w_fp[5], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[21] );
    // Amplitude(s) for diagram number 34
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[9], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 9 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram35( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 35 OF 123 ***
    // Wavefunction(s) for diagram number 35
    // (none)
    // Amplitude(s) for diagram number 35
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[9], w_fp[10], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 9 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram36( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 36 OF 123 ***
    // Wavefunction(s) for diagram number 36
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[12], w_fp[2], COUPs[1], 1.0, 0., 0., w_fp[22] );
    // Amplitude(s) for diagram number 36
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[6], w_fp[5], w_fp[22], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 9 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 15 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 21 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 23 ) += amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram37( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 37 OF 123 ***
    // Wavefunction(s) for diagram number 37
    // (none)
    // Amplitude(s) for diagram number 37
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[2], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 9 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 15 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram38( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 38 OF 123 ***
    // Wavefunction(s) for diagram number 38
    // (none)
    // Amplitude(s) for diagram number 38
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[14], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 21 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 23 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram39( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 39 OF 123 ***
    // Wavefunction(s) for diagram number 39
    // (none)
    // Amplitude(s) for diagram number 39
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[18], w_fp[4], w_fp[22], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 11 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 15 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 17 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 21 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram40( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 40 OF 123 ***
    // Wavefunction(s) for diagram number 40
    // (none)
    // Amplitude(s) for diagram number 40
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[20], w_fp[2], w_fp[18], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 11 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 21 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram41( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 41 OF 123 ***
    // Wavefunction(s) for diagram number 41
    // (none)
    // Amplitude(s) for diagram number 41
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[11], w_fp[18], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 15 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 17 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram42( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 42 OF 123 ***
    // Wavefunction(s) for diagram number 42
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[12], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[23] );
    // Amplitude(s) for diagram number 42
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[23], w_fp[11], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 17 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram43( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 43 OF 123 ***
    // Wavefunction(s) for diagram number 43
    // (none)
    // Amplitude(s) for diagram number 43
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[11], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 15 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram44( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 44 OF 123 ***
    // Wavefunction(s) for diagram number 44
    // (none)
    // Amplitude(s) for diagram number 44
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[23], w_fp[14], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 23 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram45( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 45 OF 123 ***
    // Wavefunction(s) for diagram number 45
    // (none)
    // Amplitude(s) for diagram number 45
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[20], w_fp[14], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 21 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram46( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 46 OF 123 ***
    // Wavefunction(s) for diagram number 46
    // (none)
    // Amplitude(s) for diagram number 46
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[23], w_fp[2], w_fp[10], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 17 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 23 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram47( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 47 OF 123 ***
    // Wavefunction(s) for diagram number 47
    // (none)
    // Amplitude(s) for diagram number 47
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[10], w_fp[22], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 9 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 17 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 23 ) += amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram48( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 48 OF 123 ***
    // Wavefunction(s) for diagram number 48
    // (none)
    // Amplitude(s) for diagram number 48
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[2], w_fp[17], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 9 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 17 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 23 ) += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[2], w_fp[19], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 11 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 15 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 17 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 21 ) += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[2], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 9 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 15 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 21 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 23 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram49( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 49 OF 123 ***
    // Wavefunction(s) for diagram number 49
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[4], COUPs[0], 1.0, 0., 0., w_fp[12] );
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[12], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[22] );
    // Amplitude(s) for diagram number 49
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[22], w_fp[9], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 10 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram50( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 50 OF 123 ***
    // Wavefunction(s) for diagram number 50
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[12], w_fp[5], COUPs[0], 1.0, 0., 0., w_fp[23] );
    // Amplitude(s) for diagram number 50
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[9], w_fp[23], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 6 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 8 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 10 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) += amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram51( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 51 OF 123 ***
    // Wavefunction(s) for diagram number 51
    // (none)
    // Amplitude(s) for diagram number 51
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[9], w_fp[12], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 6 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 8 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram52( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 52 OF 123 ***
    // Wavefunction(s) for diagram number 52
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[12], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[20] );
    // Amplitude(s) for diagram number 52
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[20], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 3 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 13 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram53( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 53 OF 123 ***
    // Wavefunction(s) for diagram number 53
    // (none)
    // Amplitude(s) for diagram number 53
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[2], w_fp[23], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 3 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 13 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 19 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 22 ) += amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram54( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 54 OF 123 ***
    // Wavefunction(s) for diagram number 54
    // (none)
    // Amplitude(s) for diagram number 54
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[14], w_fp[12], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 19 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 22 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram55( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 55 OF 123 ***
    // Wavefunction(s) for diagram number 55
    // (none)
    // Amplitude(s) for diagram number 55
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[20], w_fp[18], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 2 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 3 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 12 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 13 ) += amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram56( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 56 OF 123 ***
    // Wavefunction(s) for diagram number 56
    // (none)
    // Amplitude(s) for diagram number 56
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[22], w_fp[2], w_fp[18], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 10 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 20 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 21 ) += amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram57( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 57 OF 123 ***
    // Wavefunction(s) for diagram number 57
    // (none)
    // Amplitude(s) for diagram number 57
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[18], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 2 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 3 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 10 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 12 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 13 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 20 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 21 ) += cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram58( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 58 OF 123 ***
    // Wavefunction(s) for diagram number 58
    // (none)
    // Amplitude(s) for diagram number 58
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[1], w_fp[7], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 2 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 6 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 8 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 12 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 19 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 20 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 21 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 22 ) += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[1], w_fp[7], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 2 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 3 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 10 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 12 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 13 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 20 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 21 ) -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[1], w_fp[7], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 3 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 6 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 8 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 10 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 13 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 19 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 22 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram59( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 59 OF 123 ***
    // Wavefunction(s) for diagram number 59
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[12], w_fp[1], COUPs[0], 1.0, 0., 0., w_fp[21] );
    // Amplitude(s) for diagram number 59
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[7], w_fp[5], w_fp[21], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 2 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 6 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 8 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 12 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 19 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 20 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 21 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 22 ) += cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram60( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 60 OF 123 ***
    // Wavefunction(s) for diagram number 60
    // (none)
    // Amplitude(s) for diagram number 60
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[7], w_fp[23], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 3 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 6 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 8 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 10 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 13 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 19 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 22 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram61( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 61 OF 123 ***
    // Wavefunction(s) for diagram number 61
    // (none)
    // Amplitude(s) for diagram number 61
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[14], w_fp[21], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 19 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 20 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 21 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 22 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram62( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 62 OF 123 ***
    // Wavefunction(s) for diagram number 62
    // (none)
    // Amplitude(s) for diagram number 62
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[22], w_fp[14], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 20 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 21 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram63( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 63 OF 123 ***
    // Wavefunction(s) for diagram number 63
    // (none)
    // Amplitude(s) for diagram number 63
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[2], w_fp[21], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 2 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 6 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 8 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 12 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram64( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 64 OF 123 ***
    // Wavefunction(s) for diagram number 64
    // (none)
    // Amplitude(s) for diagram number 64
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[20], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 2 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 12 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram65( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 65 OF 123 ***
    // Wavefunction(s) for diagram number 65
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[5], COUPs[0], 1.0, 0., 0., w_fp[20] );
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[20], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[21] );
    // Amplitude(s) for diagram number 65
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[9], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 8 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 9 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram66( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 66 OF 123 ***
    // Wavefunction(s) for diagram number 66
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[20], w_fp[4], COUPs[0], 1.0, 0., 0., w_fp[22] );
    // Amplitude(s) for diagram number 66
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[9], w_fp[22], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 7 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 8 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 9 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 10 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram67( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 67 OF 123 ***
    // Wavefunction(s) for diagram number 67
    // (none)
    // Amplitude(s) for diagram number 67
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[9], w_fp[20], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 7 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 10 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram68( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 68 OF 123 ***
    // Wavefunction(s) for diagram number 68
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[20], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[23] );
    // Amplitude(s) for diagram number 68
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[23], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 5 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 19 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram69( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 69 OF 123 ***
    // Wavefunction(s) for diagram number 69
    // (none)
    // Amplitude(s) for diagram number 69
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[2], w_fp[22], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 5 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 13 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 16 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 19 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram70( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 70 OF 123 ***
    // Wavefunction(s) for diagram number 70
    // (none)
    // Amplitude(s) for diagram number 70
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[11], w_fp[20], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 13 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 16 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram71( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 71 OF 123 ***
    // Wavefunction(s) for diagram number 71
    // (none)
    // Amplitude(s) for diagram number 71
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[23], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 4 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 5 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 18 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 19 ) += amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram72( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 72 OF 123 ***
    // Wavefunction(s) for diagram number 72
    // (none)
    // Amplitude(s) for diagram number 72
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[2], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 8 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 9 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 14 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 15 ) += amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram73( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 73 OF 123 ***
    // Wavefunction(s) for diagram number 73
    // (none)
    // Amplitude(s) for diagram number 73
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[20], w_fp[6], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 4 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 5 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 8 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 9 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 14 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 15 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 18 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 19 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram74( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 74 OF 123 ***
    // Wavefunction(s) for diagram number 74
    // (none)
    // Amplitude(s) for diagram number 74
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[20], w_fp[1], w_fp[7], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 4 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 7 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 10 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 13 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 14 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 15 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 16 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 18 ) -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[20], w_fp[1], w_fp[7], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 4 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 5 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 8 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 9 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 14 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 15 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 18 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 19 ) += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[20], w_fp[1], w_fp[7], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 5 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 7 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 8 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 9 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 10 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 13 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 16 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 19 ) += cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram75( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 75 OF 123 ***
    // Wavefunction(s) for diagram number 75
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[20], w_fp[1], COUPs[0], 1.0, 0., 0., w_fp[12] );
    // Amplitude(s) for diagram number 75
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[7], w_fp[4], w_fp[12], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 4 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 7 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 10 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 13 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 14 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 15 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 16 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 18 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram76( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 76 OF 123 ***
    // Wavefunction(s) for diagram number 76
    // (none)
    // Amplitude(s) for diagram number 76
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[7], w_fp[22], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 5 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 7 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 8 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 9 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 10 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 13 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 16 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 19 ) += cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram77( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 77 OF 123 ***
    // Wavefunction(s) for diagram number 77
    // (none)
    // Amplitude(s) for diagram number 77
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[11], w_fp[12], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 13 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 14 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 15 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 16 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram78( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 78 OF 123 ***
    // Wavefunction(s) for diagram number 78
    // (none)
    // Amplitude(s) for diagram number 78
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[11], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 14 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 15 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram79( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 79 OF 123 ***
    // Wavefunction(s) for diagram number 79
    // (none)
    // Amplitude(s) for diagram number 79
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[2], w_fp[12], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 4 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 7 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 10 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 18 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram80( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 80 OF 123 ***
    // Wavefunction(s) for diagram number 80
    // (none)
    // Amplitude(s) for diagram number 80
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[23], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 4 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 18 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram81( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 81 OF 123 ***
    // Wavefunction(s) for diagram number 81
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[9], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[23] );
    // Amplitude(s) for diagram number 81
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[23], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 7 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram82( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 82 OF 123 ***
    // Wavefunction(s) for diagram number 82
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[15], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[12] );
    // Amplitude(s) for diagram number 82
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[9], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 10 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram83( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 83 OF 123 ***
    // Wavefunction(s) for diagram number 83
    // (none)
    // Amplitude(s) for diagram number 83
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[23], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 6 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram84( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 84 OF 123 ***
    // Wavefunction(s) for diagram number 84
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[13], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[21] );
    // Amplitude(s) for diagram number 84
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[9], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 8 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram85( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 85 OF 123 ***
    // Wavefunction(s) for diagram number 85
    // (none)
    // Amplitude(s) for diagram number 85
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[23], w_fp[10], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 6 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 7 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram86( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 86 OF 123 ***
    // Wavefunction(s) for diagram number 86
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[10], COUPs[0], 1.0, 0., 0., w_fp[23] );
    // Amplitude(s) for diagram number 86
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[9], w_fp[23], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 6 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 7 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 9 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) += amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram87( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 87 OF 123 ***
    // Wavefunction(s) for diagram number 87
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[16], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[22] );
    // Amplitude(s) for diagram number 87
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[22], w_fp[11], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 16 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram88( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 88 OF 123 ***
    // Wavefunction(s) for diagram number 88
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[11], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[20] );
    // Amplitude(s) for diagram number 88
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[20], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 13 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram89( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 89 OF 123 ***
    // Wavefunction(s) for diagram number 89
    // (none)
    // Amplitude(s) for diagram number 89
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[22], w_fp[14], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 22 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram90( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 90 OF 123 ***
    // Wavefunction(s) for diagram number 90
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[14], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[24] );
    // Amplitude(s) for diagram number 90
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[24], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 19 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram91( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 91 OF 123 ***
    // Wavefunction(s) for diagram number 91
    // (none)
    // Amplitude(s) for diagram number 91
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[22], w_fp[2], w_fp[10], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 16 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 22 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram92( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 92 OF 123 ***
    // Wavefunction(s) for diagram number 92
    // (none)
    // Amplitude(s) for diagram number 92
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[2], w_fp[23], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 3 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 5 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 16 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 22 ) += amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram93( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 93 OF 123 ***
    // Wavefunction(s) for diagram number 93
    // (none)
    // Amplitude(s) for diagram number 93
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[6], w_fp[7], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 2 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 8 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 14 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 18 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 19 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 21 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 23 ) -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[6], w_fp[7], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 2 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 4 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 5 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 9 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 15 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 21 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 23 ) -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[6], w_fp[7], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 4 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 5 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 8 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 9 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 14 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 15 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 18 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 19 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram94( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 94 OF 123 ***
    // Wavefunction(s) for diagram number 94
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[6], COUPs[0], 1.0, 0., 0., w_fp[22] );
    // Amplitude(s) for diagram number 94
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[7], w_fp[5], w_fp[22], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 2 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 8 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 14 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 18 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 19 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 21 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 23 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram95( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 95 OF 123 ***
    // Wavefunction(s) for diagram number 95
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[7], COUPs[0], 1.0, 0., 0., w_fp[25] );
    // Amplitude(s) for diagram number 95
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[6], w_fp[5], w_fp[25], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 2 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 4 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 5 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 9 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 15 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 21 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 23 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram96( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 96 OF 123 ***
    // Wavefunction(s) for diagram number 96
    // (none)
    // Amplitude(s) for diagram number 96
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[14], w_fp[22], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 18 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 19 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 21 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 23 ) += amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram97( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 97 OF 123 ***
    // Wavefunction(s) for diagram number 97
    // (none)
    // Amplitude(s) for diagram number 97
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[24], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 18 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 19 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram98( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 98 OF 123 ***
    // Wavefunction(s) for diagram number 98
    // (none)
    // Amplitude(s) for diagram number 98
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[2], w_fp[22], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 2 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 8 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 14 ) += amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram99( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 99 OF 123 ***
    // Wavefunction(s) for diagram number 99
    // (none)
    // Amplitude(s) for diagram number 99
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[2], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 8 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 14 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram100( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 100 OF 123 ***
    // Wavefunction(s) for diagram number 100
    // (none)
    // Amplitude(s) for diagram number 100
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[18], w_fp[7], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 1 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 4 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 10 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 12 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 13 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 15 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 17 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 20 ) += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[18], w_fp[7], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 1 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 2 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 3 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 4 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 15 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 17 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 21 ) += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[18], w_fp[7], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 2 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 3 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 10 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 12 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 13 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 20 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 21 ) += cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram101( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 101 OF 123 ***
    // Wavefunction(s) for diagram number 101
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[18], COUPs[0], 1.0, 0., 0., w_fp[6] );
    // Amplitude(s) for diagram number 101
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[7], w_fp[4], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 1 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 4 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 10 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 12 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 13 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 15 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 17 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 20 ) += cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram102( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 102 OF 123 ***
    // Wavefunction(s) for diagram number 102
    // (none)
    // Amplitude(s) for diagram number 102
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[18], w_fp[4], w_fp[25], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 1 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 2 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 3 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 4 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 15 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 17 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 21 ) += cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram103( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 103 OF 123 ***
    // Wavefunction(s) for diagram number 103
    // (none)
    // Amplitude(s) for diagram number 103
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[11], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 12 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 13 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 15 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 17 ) += amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram104( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 104 OF 123 ***
    // Wavefunction(s) for diagram number 104
    // (none)
    // Amplitude(s) for diagram number 104
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[20], w_fp[18], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 12 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 13 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram105( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 105 OF 123 ***
    // Wavefunction(s) for diagram number 105
    // (none)
    // Amplitude(s) for diagram number 105
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[2], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 1 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 4 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 10 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 20 ) += amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram106( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 106 OF 123 ***
    // Wavefunction(s) for diagram number 106
    // (none)
    // Amplitude(s) for diagram number 106
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[2], w_fp[18], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 10 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 20 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram107( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 107 OF 123 ***
    // Wavefunction(s) for diagram number 107
    // (none)
    // Amplitude(s) for diagram number 107
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[7], w_fp[10], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 1 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 6 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 7 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 16 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 17 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 22 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 23 ) -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[7], w_fp[10], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 1 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 3 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 5 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 9 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 17 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 23 ) -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[7], w_fp[10], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 3 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 5 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 6 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 7 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 9 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 16 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 22 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram108( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 108 OF 123 ***
    // Wavefunction(s) for diagram number 108
    // (none)
    // Amplitude(s) for diagram number 108
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[10], w_fp[25], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 1 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 3 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 5 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 9 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 17 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 23 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram109( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 109 OF 123 ***
    // Wavefunction(s) for diagram number 109
    // (none)
    // Amplitude(s) for diagram number 109
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[7], w_fp[23], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 3 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 5 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 6 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 7 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 9 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 16 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 22 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram110( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 110 OF 123 ***
    // Wavefunction(s) for diagram number 110
    // (none)
    // Amplitude(s) for diagram number 110
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[20], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 12 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram111( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 111 OF 123 ***
    // Wavefunction(s) for diagram number 111
    // (none)
    // Amplitude(s) for diagram number 111
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[11], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 14 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram112( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 112 OF 123 ***
    // Wavefunction(s) for diagram number 112
    // (none)
    // Amplitude(s) for diagram number 112
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[24], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 18 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram113( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 113 OF 123 ***
    // Wavefunction(s) for diagram number 113
    // (none)
    // Amplitude(s) for diagram number 113
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[14], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 20 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram114( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 114 OF 123 ***
    // Wavefunction(s) for diagram number 114
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[4], COUPs[2], 1.0, 0., 0., w_fp[12] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[4], COUPs[2], 1.0, 0., 0., w_fp[24] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[4], COUPs[2], 1.0, 0., 0., w_fp[21] );
    // Amplitude(s) for diagram number 114
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[7], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 2 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 8 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 14 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 18 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 19 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 21 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 23 ) -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[24], w_fp[7], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 2 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 6 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 8 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 12 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 19 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 20 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 21 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 22 ) -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[7], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 6 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 12 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 14 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 18 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 20 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 22 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 23 ) += cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram115( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 115 OF 123 ***
    // Wavefunction(s) for diagram number 115
    // (none)
    // Amplitude(s) for diagram number 115
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[14], w_fp[12], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 18 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 19 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 21 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 23 ) += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[14], w_fp[24], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 19 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 20 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 21 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 22 ) += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[14], w_fp[21], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 18 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 20 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 22 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 23 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram116( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 116 OF 123 ***
    // Wavefunction(s) for diagram number 116
    // (none)
    // Amplitude(s) for diagram number 116
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[2], w_fp[12], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 2 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 8 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 14 ) += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[2], w_fp[24], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 2 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 6 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 8 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 12 ) += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[2], w_fp[21], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 6 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 12 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 14 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram117( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 117 OF 123 ***
    // Wavefunction(s) for diagram number 117
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[5], COUPs[2], 1.0, 0., 0., w_fp[21] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[5], COUPs[2], 1.0, 0., 0., w_fp[13] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[5], COUPs[2], 1.0, 0., 0., w_fp[24] );
    // Amplitude(s) for diagram number 117
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[7], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 1 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 4 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 10 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 12 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 13 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 15 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 17 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 20 ) += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[7], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 4 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 7 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 10 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 13 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 14 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 15 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 16 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 18 ) += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[24], w_fp[7], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 1 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 7 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 12 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 14 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 16 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 17 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 18 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 20 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram118( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 118 OF 123 ***
    // Wavefunction(s) for diagram number 118
    // (none)
    // Amplitude(s) for diagram number 118
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[11], w_fp[21], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 12 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 13 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 15 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 17 ) += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[11], w_fp[13], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 13 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 14 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 15 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 16 ) += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[11], w_fp[24], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 12 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 14 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 16 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 17 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram119( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 119 OF 123 ***
    // Wavefunction(s) for diagram number 119
    // (none)
    // Amplitude(s) for diagram number 119
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[2], w_fp[21], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 1 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 4 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 10 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 20 ) += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[2], w_fp[13], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 4 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 7 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 10 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 18 ) += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[2], w_fp[24], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 1 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 7 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 18 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 20 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram120( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 120 OF 123 ***
    // Wavefunction(s) for diagram number 120
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[4], w_fp[5], COUPs[2], 1.0, 0., 0., w_fp[24] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[4], w_fp[5], COUPs[2], 1.0, 0., 0., w_fp[15] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[4], w_fp[5], COUPs[2], 1.0, 0., 0., w_fp[13] );
    // Amplitude(s) for diagram number 120
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[9], w_fp[24], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 6 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 7 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 9 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[9], w_fp[15], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 7 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 8 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 9 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 10 ) += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[9], w_fp[13], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 6 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 8 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 10 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram121( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 121 OF 123 ***
    // Wavefunction(s) for diagram number 121
    // (none)
    // Amplitude(s) for diagram number 121
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[2], w_fp[24], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 3 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 5 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 16 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 22 ) += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[2], w_fp[15], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 5 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 13 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 16 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 19 ) += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[2], w_fp[13], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 3 ) -= amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 13 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 19 ) += amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 22 ) -= amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram122( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 122 OF 123 ***
    // Wavefunction(s) for diagram number 122
    // (none)
    // Amplitude(s) for diagram number 122
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[24], w_fp[1], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 3 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 5 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 6 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 7 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 9 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 16 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 22 ) -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[1], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 5 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 7 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 8 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 9 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 10 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 13 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 16 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 19 ) -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[1], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 3 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 6 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 8 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 10 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 13 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 19 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 22 ) += cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ INLINE void
  diagram123( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    // *** DIAGRAM 123 OF 123 ***
    // Wavefunction(s) for diagram number 123
    // (none)
    // Amplitude(s) for diagram number 123
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[17], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 1 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 3 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 5 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 9 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 17 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 23 ) += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[19], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 1 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 2 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 3 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 4 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 15 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 17 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 21 ) += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[8], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 2 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 4 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 5 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 9 ) -= cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 15 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 21 ) += cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 23 ) -= cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

/* clang-format on */
