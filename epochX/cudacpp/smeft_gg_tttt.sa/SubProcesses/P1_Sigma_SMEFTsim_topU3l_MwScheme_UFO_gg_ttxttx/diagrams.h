// Copyright (C) 2020-2025 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Sep 2025) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2025) for the MG5aMC CUDACPP plugin.

/* clang-format off */

  //--------------------------------------------------------------------------

  __global__ void
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
    J_ACCESS::kernelAccessIcol( jamps, 2 ) -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 3 ) += 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 10 ) += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) -= 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 2 OF 72 ***
    // Wavefunction(s) for diagram number 2
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[5], w_fp[6], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[9] );
    // Amplitude(s) for diagram number 2
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[4], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 1 ) -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 3 ) += 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 9 ) += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) -= 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 3 OF 72 ***
    // Wavefunction(s) for diagram number 3
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[5], w_fp[4], COUPs[1], 1.0, 0., 0., w_fp[10] );
    // Amplitude(s) for diagram number 3
    VVV5_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[6], w_fp[7], w_fp[10], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 1 ) += 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 2 ) -= 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 9 ) -= 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 10 ) += 1. / 2. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 4 OF 72 ***
    // Wavefunction(s) for diagram number 4
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[5], w_fp[2], COUPs[1], 1.0, 0., 0., w_fp[11] );
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[6], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[12] );
    // Amplitude(s) for diagram number 4
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[4], w_fp[11], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 2 ) -= 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 8 ) -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 10 ) += 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 5 OF 72 ***
    // Wavefunction(s) for diagram number 5
    // (none)
    // Amplitude(s) for diagram number 5
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[8], w_fp[11], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 2 ) -= 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 3 ) += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 10 ) += 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 6 OF 72 ***
    // Wavefunction(s) for diagram number 6
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[4], COUPs[1], 1.0, 0., 0., w_fp[8] );
    // Amplitude(s) for diagram number 6
    VVV5_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[6], w_fp[11], w_fp[8], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) -= 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 3 ) += 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 8 ) += 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) -= 1. / 2. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 7 OF 72 ***
    // Wavefunction(s) for diagram number 7
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[6], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[13] );
    // Amplitude(s) for diagram number 7
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[13], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 1 ) -= 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 8 ) -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 9 ) += 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 8 OF 72 ***
    // Wavefunction(s) for diagram number 8
    // (none)
    // Amplitude(s) for diagram number 8
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[2], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 1 ) -= 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 3 ) += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 9 ) += 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 9 OF 72 ***
    // Wavefunction(s) for diagram number 9
    // (none)
    // Amplitude(s) for diagram number 9
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[13], w_fp[10], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) += 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 1 ) -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 8 ) -= 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 9 ) += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 10 OF 72 ***
    // Wavefunction(s) for diagram number 10
    // (none)
    // Amplitude(s) for diagram number 10
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[2], w_fp[10], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) += 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 2 ) -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 8 ) -= 1. / 6. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 10 ) += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    J_ACCESS::kernelAccessIcol( jamps, 0 ) -= 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 5 ) += 1. / 6. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 12 OF 72 ***
    // Wavefunction(s) for diagram number 12
    // (none)
    // Amplitude(s) for diagram number 12
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[12], w_fp[10], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) -= 1. / 6. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 5 ) += 1. / 2. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 13 OF 72 ***
    // Wavefunction(s) for diagram number 13
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[4], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[6] );
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[12], COUPs[1], 1.0, 0., 0., w_fp[14] );
    // Amplitude(s) for diagram number 13
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[6], w_fp[14], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 4 ) -= 1. / 6. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 5 ) += 1. / 2. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 14 OF 72 ***
    // Wavefunction(s) for diagram number 14
    // (none)
    // Amplitude(s) for diagram number 14
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[6], w_fp[9], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 4 ) -= 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 5 ) += 1. / 6. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 15 OF 72 ***
    // Wavefunction(s) for diagram number 15
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[5], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[15] );
    // Amplitude(s) for diagram number 15
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[4], w_fp[14], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 1 ) += 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 4 ) -= 1. / 6. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 16 OF 72 ***
    // Wavefunction(s) for diagram number 16
    // (none)
    // Amplitude(s) for diagram number 16
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[12], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 1 ) += 1. / 6. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 4 ) -= 1. / 2. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 17 OF 72 ***
    // Wavefunction(s) for diagram number 17
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[12], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[16] );
    // Amplitude(s) for diagram number 17
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[16], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) -= 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 1 ) += 1. / 6. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 18 OF 72 ***
    // Wavefunction(s) for diagram number 18
    // (none)
    // Amplitude(s) for diagram number 18
    VVV5_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[8], w_fp[9], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 4 ) -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 19 OF 72 ***
    // Wavefunction(s) for diagram number 19
    // (none)
    // Amplitude(s) for diagram number 19
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[16], w_fp[10], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) -= 1. / 6. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 1 ) += 1. / 2. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 20 OF 72 ***
    // Wavefunction(s) for diagram number 20
    // (none)
    // Amplitude(s) for diagram number 20
    VVV5_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[10], w_fp[14], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 1 ) -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 5 ) += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    J_ACCESS::kernelAccessIcol( jamps, 6 ) += 1. / 6. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 8 ) -= 1. / 2. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 22 OF 72 ***
    // Wavefunction(s) for diagram number 22
    // (none)
    // Amplitude(s) for diagram number 22
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[14], w_fp[16], w_fp[10], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 6 ) += 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 8 ) -= 1. / 6. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 23 OF 72 ***
    // Wavefunction(s) for diagram number 23
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[14], w_fp[2], COUPs[1], 1.0, 0., 0., w_fp[12] );
    // Amplitude(s) for diagram number 23
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[6], w_fp[12], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 4 ) -= 1. / 6. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 10 ) += 1. / 2. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 24 OF 72 ***
    // Wavefunction(s) for diagram number 24
    // (none)
    // Amplitude(s) for diagram number 24
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[14], w_fp[6], w_fp[11], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 4 ) -= 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 10 ) += 1. / 6. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 25 OF 72 ***
    // Wavefunction(s) for diagram number 25
    // (none)
    // Amplitude(s) for diagram number 25
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[4], w_fp[12], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 4 ) -= 1. / 6. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 6 ) += 1. / 2. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 26 OF 72 ***
    // Wavefunction(s) for diagram number 26
    // (none)
    // Amplitude(s) for diagram number 26
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[2], w_fp[9], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 4 ) -= 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 6 ) += 1. / 6. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 27 OF 72 ***
    // Wavefunction(s) for diagram number 27
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[14], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[17] );
    // Amplitude(s) for diagram number 27
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[17], w_fp[4], w_fp[11], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 8 ) -= 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 10 ) += 1. / 6. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 28 OF 72 ***
    // Wavefunction(s) for diagram number 28
    // (none)
    // Amplitude(s) for diagram number 28
    VVV5_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[11], w_fp[9], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 4 ) += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 8 ) -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 29 OF 72 ***
    // Wavefunction(s) for diagram number 29
    // (none)
    // Amplitude(s) for diagram number 29
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[17], w_fp[2], w_fp[10], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 8 ) -= 1. / 6. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 10 ) += 1. / 2. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 30 OF 72 ***
    // Wavefunction(s) for diagram number 30
    // (none)
    // Amplitude(s) for diagram number 30
    VVV5_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[10], w_fp[12], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 6 ) -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 10 ) += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 31 OF 72 ***
    // Wavefunction(s) for diagram number 31
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[4], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[12] );
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[12], COUPs[1], 1.0, 0., 0., w_fp[17] );
    // Amplitude(s) for diagram number 31
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[16], w_fp[17], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 6 ) += 1. / 6. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 7 ) -= 1. / 2. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 32 OF 72 ***
    // Wavefunction(s) for diagram number 32
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[5], w_fp[12], COUPs[1], 1.0, 0., 0., w_fp[9] );
    // Amplitude(s) for diagram number 32
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[16], w_fp[9], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 6 ) += 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 7 ) -= 1. / 6. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 33 OF 72 ***
    // Wavefunction(s) for diagram number 33
    // (none)
    // Amplitude(s) for diagram number 33
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[2], w_fp[9], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 2 ) += 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 7 ) -= 1. / 6. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 34 OF 72 ***
    // Wavefunction(s) for diagram number 34
    // (none)
    // Amplitude(s) for diagram number 34
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[12], w_fp[11], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 2 ) += 1. / 6. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 7 ) -= 1. / 2. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 35 OF 72 ***
    // Wavefunction(s) for diagram number 35
    // (none)
    // Amplitude(s) for diagram number 35
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[2], w_fp[17], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 3 ) -= 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 6 ) += 1. / 6. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 36 OF 72 ***
    // Wavefunction(s) for diagram number 36
    // (none)
    // Amplitude(s) for diagram number 36
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[12], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 3 ) -= 1. / 6. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 6 ) += 1. / 2. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 37 OF 72 ***
    // Wavefunction(s) for diagram number 37
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[12], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[14] );
    // Amplitude(s) for diagram number 37
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[14], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 2 ) += 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 3 ) -= 1. / 6. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 38 OF 72 ***
    // Wavefunction(s) for diagram number 38
    // (none)
    // Amplitude(s) for diagram number 38
    VVV5_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[7], w_fp[9], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 2 ) -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 6 ) += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 39 OF 72 ***
    // Wavefunction(s) for diagram number 39
    // (none)
    // Amplitude(s) for diagram number 39
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[14], w_fp[11], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 2 ) += 1. / 6. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 3 ) -= 1. / 2. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 40 OF 72 ***
    // Wavefunction(s) for diagram number 40
    // (none)
    // Amplitude(s) for diagram number 40
    VVV5_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[11], w_fp[17], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 3 ) += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 7 ) -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 41 OF 72 ***
    // Wavefunction(s) for diagram number 41
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[5], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[17] );
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[17], w_fp[4], COUPs[1], 1.0, 0., 0., w_fp[14] );
    // Amplitude(s) for diagram number 41
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[16], w_fp[14], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 7 ) -= 1. / 6. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 9 ) += 1. / 2. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 42 OF 72 ***
    // Wavefunction(s) for diagram number 42
    // (none)
    // Amplitude(s) for diagram number 42
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[17], w_fp[16], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 7 ) -= 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 9 ) += 1. / 6. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 43 OF 72 ***
    // Wavefunction(s) for diagram number 43
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[17], w_fp[2], COUPs[1], 1.0, 0., 0., w_fp[9] );
    // Amplitude(s) for diagram number 43
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[4], w_fp[9], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 5 ) += 1. / 6. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 7 ) -= 1. / 2. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 44 OF 72 ***
    // Wavefunction(s) for diagram number 44
    // (none)
    // Amplitude(s) for diagram number 44
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[2], w_fp[14], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 5 ) += 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 7 ) -= 1. / 6. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 45 OF 72 ***
    // Wavefunction(s) for diagram number 45
    // (none)
    // Amplitude(s) for diagram number 45
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[6], w_fp[9], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 5 ) += 1. / 6. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) -= 1. / 2. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 46 OF 72 ***
    // Wavefunction(s) for diagram number 46
    // (none)
    // Amplitude(s) for diagram number 46
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[17], w_fp[6], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 5 ) += 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) -= 1. / 6. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 47 OF 72 ***
    // Wavefunction(s) for diagram number 47
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[17], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[12] );
    // Amplitude(s) for diagram number 47
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[4], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 9 ) += 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) -= 1. / 6. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 48 OF 72 ***
    // Wavefunction(s) for diagram number 48
    // (none)
    // Amplitude(s) for diagram number 48
    VVV5_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[7], w_fp[14], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 5 ) -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 9 ) += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 49 OF 72 ***
    // Wavefunction(s) for diagram number 49
    // (none)
    // Amplitude(s) for diagram number 49
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[2], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 9 ) += 1. / 6. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) -= 1. / 2. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 50 OF 72 ***
    // Wavefunction(s) for diagram number 50
    // (none)
    // Amplitude(s) for diagram number 50
    VVV5_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[8], w_fp[9], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 7 ) += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 51 OF 72 ***
    // Wavefunction(s) for diagram number 51
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[16], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[9] );
    // Amplitude(s) for diagram number 51
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[9], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 8 ) -= 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 9 ) += 1. / 6. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 52 OF 72 ***
    // Wavefunction(s) for diagram number 52
    VVV5P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[8], COUPs[0], 1.0, 0., 0., w_fp[12] );
    // Amplitude(s) for diagram number 52
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[16], w_fp[12], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 7 ) -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 8 ) += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 53 OF 72 ***
    // Wavefunction(s) for diagram number 53
    // (none)
    // Amplitude(s) for diagram number 53
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[9], w_fp[10], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 8 ) -= 1. / 6. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 9 ) += 1. / 2. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 54 OF 72 ***
    // Wavefunction(s) for diagram number 54
    VVV5P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[10], COUPs[0], 1.0, 0., 0., w_fp[9] );
    // Amplitude(s) for diagram number 54
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[16], w_fp[9], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 6 ) += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 9 ) -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 55 OF 72 ***
    // Wavefunction(s) for diagram number 55
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[13], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[16] );
    // Amplitude(s) for diagram number 55
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[4], w_fp[11], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) -= 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 2 ) += 1. / 6. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 56 OF 72 ***
    // Wavefunction(s) for diagram number 56
    VVV5P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[11], COUPs[0], 1.0, 0., 0., w_fp[14] );
    // Amplitude(s) for diagram number 56
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[4], w_fp[14], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 7 ) += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 57 OF 72 ***
    // Wavefunction(s) for diagram number 57
    // (none)
    // Amplitude(s) for diagram number 57
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[2], w_fp[10], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) -= 1. / 6. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 2 ) += 1. / 2. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 58 OF 72 ***
    // Wavefunction(s) for diagram number 58
    // (none)
    // Amplitude(s) for diagram number 58
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[2], w_fp[9], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 2 ) += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 5 ) -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 59 OF 72 ***
    // Wavefunction(s) for diagram number 59
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[6], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[13] );
    // Amplitude(s) for diagram number 59
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[13], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 10 ) += 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) -= 1. / 6. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 60 OF 72 ***
    // Wavefunction(s) for diagram number 60
    VVV5P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[7], COUPs[0], 1.0, 0., 0., w_fp[16] );
    // Amplitude(s) for diagram number 60
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[6], w_fp[16], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 5 ) += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 10 ) -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 61 OF 72 ***
    // Wavefunction(s) for diagram number 61
    // (none)
    // Amplitude(s) for diagram number 61
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[13], w_fp[11], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 10 ) += 1. / 6. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) -= 1. / 2. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 62 OF 72 ***
    // Wavefunction(s) for diagram number 62
    // (none)
    // Amplitude(s) for diagram number 62
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[6], w_fp[14], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 4 ) -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 63 OF 72 ***
    // Wavefunction(s) for diagram number 63
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[15], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[6] );
    // Amplitude(s) for diagram number 63
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[6], w_fp[4], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 1 ) += 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 3 ) -= 1. / 6. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 64 OF 72 ***
    // Wavefunction(s) for diagram number 64
    // (none)
    // Amplitude(s) for diagram number 64
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[4], w_fp[16], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 1 ) += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 6 ) -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 65 OF 72 ***
    // Wavefunction(s) for diagram number 65
    // (none)
    // Amplitude(s) for diagram number 65
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[6], w_fp[2], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 1 ) += 1. / 6. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 3 ) -= 1. / 2. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 66 OF 72 ***
    // Wavefunction(s) for diagram number 66
    // (none)
    // Amplitude(s) for diagram number 66
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[2], w_fp[12], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 3 ) -= 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 4 ) += 1. / 2. * cxtype( 0, 1 ) * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 67 OF 72 ***
    // Wavefunction(s) for diagram number 67
    // (none)
    // Amplitude(s) for diagram number 67
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[7], w_fp[10], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 1 ) += 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 2 ) -= 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 9 ) -= 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 10 ) += 1. / 2. * amp_sv[0];
    VVVV9_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[7], w_fp[10], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 1 ) += 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 5 ) -= 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 6 ) -= 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 10 ) += 1. / 2. * amp_sv[0];
    VVVV10_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[7], w_fp[10], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 2 ) += 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 5 ) -= 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 6 ) -= 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 9 ) += 1. / 2. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 68 OF 72 ***
    // Wavefunction(s) for diagram number 68
    // (none)
    // Amplitude(s) for diagram number 68
    VVV5_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[10], w_fp[16], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 1 ) += 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 5 ) -= 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 6 ) -= 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 10 ) += 1. / 2. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 69 OF 72 ***
    // Wavefunction(s) for diagram number 69
    // (none)
    // Amplitude(s) for diagram number 69
    VVV5_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[7], w_fp[9], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 2 ) += 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 5 ) -= 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 6 ) -= 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 9 ) += 1. / 2. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 70 OF 72 ***
    // Wavefunction(s) for diagram number 70
    // (none)
    // Amplitude(s) for diagram number 70
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[11], w_fp[8], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) -= 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 3 ) += 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 8 ) += 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) -= 1. / 2. * amp_sv[0];
    VVVV9_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[11], w_fp[8], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) -= 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 4 ) += 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 7 ) += 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) -= 1. / 2. * amp_sv[0];
    VVVV10_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[11], w_fp[8], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 3 ) -= 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 4 ) += 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 7 ) += 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 8 ) -= 1. / 2. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 71 OF 72 ***
    // Wavefunction(s) for diagram number 71
    // (none)
    // Amplitude(s) for diagram number 71
    VVV5_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[8], w_fp[14], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 0 ) -= 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 4 ) += 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 7 ) += 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 11 ) -= 1. / 2. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

  __global__ void
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
    // *** DIAGRAM 72 OF 72 ***
    // Wavefunction(s) for diagram number 72
    // (none)
    // Amplitude(s) for diagram number 72
    VVV5_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[11], w_fp[12], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    J_ACCESS::kernelAccessIcol( jamps, 3 ) -= 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 4 ) += 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 7 ) += 1. / 2. * amp_sv[0];
    J_ACCESS::kernelAccessIcol( jamps, 8 ) -= 1. / 2. * amp_sv[0];
  }
  
  //--------------------------------------------------------------------------

/* clang-format on */
