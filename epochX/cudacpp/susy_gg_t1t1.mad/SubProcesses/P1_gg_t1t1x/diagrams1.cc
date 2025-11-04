// Copyright (C) 2020-2025 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Sep 2025) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2025) for the MG5aMC CUDACPP plugin.

#include "GpuRuntime.h"
#include "HelAmps_MSSM_SLHA2.h"
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
  diagramgroup1( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
                 const fptype* cIPD,             // input: GPU __device__ or GPU host address of cIPD
#ifndef MGONGPU_RDC_DIAGRAMS
                 const short* cHelFlat,          // input: GPU __device__ or GPU host address of cHel
#else
                 const short (*cHel)[CPPProcess::npar], // input: GPU __device__ or GPU host address of cHel
#endif
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
#ifndef MGONGPU_RDC_DIAGRAMS
    // *** RETRIEVE WAVEFUNCTIONS FROM PREVIOUS DIAGRAM GROUPS ***
    // (none)
#endif
#endif

#ifndef MGONGPU_RDC_DIAGRAMS
    // Reinterpret the flat array pointer for helicities as a multidimensional array pointer
    constexpr int npar = CPPProcess::npar;
    const short (*cHel)[npar] = reinterpret_cast<const short(*)[npar]>( cHelFlat );
#endif

    // *** DIAGRAM 1 OF 6 ***
    // Wavefunction(s) for diagram number 1
    vxxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][0], -1, w_fp[0], 0 );
    vxxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][1], -1, w_fp[1], 1 );
    sxxxxx<M_ACCESS, W_ACCESS>( momenta, +1, w_fp[2], 2 );
    sxxxxx<M_ACCESS, W_ACCESS>( momenta, +1, w_fp[3], 3 );
    // Amplitude(s) for diagram number 1
    VVSS1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[3], w_fp[2], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[1] += amp_sv[0];
    VVSS1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[3], w_fp[2], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[0] += amp_sv[0];

    // *** DIAGRAM 2 OF 6 ***
    // Wavefunction(s) for diagram number 2
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], COUPs[1], 1.0, 0., 0., w_fp[4] );
    // Amplitude(s) for diagram number 2
    VSS1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[4], w_fp[3], w_fp[2], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 2 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3 OF 6 ***
    // Wavefunction(s) for diagram number 3
    VSS1_2<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[2], COUPs[2], 1.0, cIPD[0], cIPD[1], w_fp[4] );
    // Amplitude(s) for diagram number 3
    VSS1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[3], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 3 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += amp_sv[0];

    // *** DIAGRAM 4 OF 6 ***
    // Wavefunction(s) for diagram number 4
    VSS1_3<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[2], COUPs[3], -1.0, cIPD[2], cIPD[3], w_fp[4] );
    // Amplitude(s) for diagram number 4
    VSS1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[3], w_fp[4], COUPs[3], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 4 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += amp_sv[0];

    // *** DIAGRAM 5 OF 6 ***
    // Wavefunction(s) for diagram number 5
    VSS1_3<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[3], COUPs[2], 1.0, cIPD[0], cIPD[1], w_fp[4] );
    // Amplitude(s) for diagram number 5
    VSS1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[4], w_fp[2], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 5 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[1] += amp_sv[0];

    // *** DIAGRAM 6 OF 6 ***
    // Wavefunction(s) for diagram number 6
    VSS1_3<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[3], COUPs[3], 1.0, cIPD[2], cIPD[3], w_fp[4] );
    // Amplitude(s) for diagram number 6
    VSS1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[2], w_fp[4], COUPs[3], -1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 6 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[1] += amp_sv[0];

#if defined MGONGPUCPP_GPUIMPL and not defined MGONGPU_RDC_DIAGRAMS
    // *** STORE JAMPS ***
    // In CUDA (DCDIAG=0), copy the local jamp to the output global-memory jamp
    //printf( "diagramgroup1: nGoodHel=%d\n", nGoodHel );
    constexpr int ihel0 = 0; // allJamps buffer points to a specific helicity _within a super-buffer for nGoodHel helicities_
    for( int icol = 0; icol < ncolor; icol++ )
      J_ACCESS::kernelAccessIcolIhelNhel( jamps, icol, ihel0, nGoodHel ) = jamp_sv[icol]; // set jamps
#else
    // In C++, copy the local jamp to the output array passed as function argument
    // In CUDA (DCDIAG=1), copy the local jamp to the output array passed as function argument
    for( int icol = 0; icol < ncolor; icol++ )
      jamps[icol] = jamp_sv[icol]; // set jamps
#endif

#ifdef MGONGPUCPP_GPUIMPL
#ifndef MGONGPU_RDC_DIAGRAMS
    // *** STORE WAVEFUNCTIONS FOR NEXT DIAGRAM GROUPS ***
    // (none)
#endif
#endif
  }

  //--------------------------------------------------------------------------
}
