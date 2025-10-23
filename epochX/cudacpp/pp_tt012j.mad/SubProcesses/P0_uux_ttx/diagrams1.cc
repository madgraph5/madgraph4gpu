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

    // *** DIAGRAM 1 OF 1 ***
    // Wavefunction(s) for diagram number 1
    ixxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][0], +1, w_fp[0], 0 );
    oxxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][1], -1, w_fp[1], 1 );
    oxxxxx<M_ACCESS, W_ACCESS>( momenta, cIPD[0], cHel[ihel][2], +1, w_fp[2], 2 );
    ixxxxx<M_ACCESS, W_ACCESS>( momenta, cIPD[0], cHel[ihel][3], -1, w_fp[3], 3 );
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], COUPs[1], 1.0, 0., 0., w_fp[4] );
    // Amplitude(s) for diagram number 1
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[2], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 1 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += 1. / 6. * amp_sv[0];
    jamp_sv[1] -= 1. / 2. * amp_sv[0];

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
