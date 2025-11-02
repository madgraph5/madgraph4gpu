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

#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup19( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
#ifdef MGONGPUCPP_GPUIMPL
                  fptype* jamps,                  // output jamps[ncolor*2*nevt]
                  const int nGoodHel,             // input: number of good helicities
                  const fptype* couplings,        // input: dependent couplings[nevt*ndcoup*2] for all events
#else
                  cxtype_sv* jamps,               // output jamps[ncolor*2*neppV]
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
#ifndef MGONGPU_RDC_DIAGRAMS
    // *** RETRIEVE WAVEFUNCTIONS FROM PREVIOUS DIAGRAM GROUPS ***
    //for( int iwf = 0; iwf < nwf; iwf++ ) retrieveWf( wfs, w_cx, nevt, iwf );
    retrieveWf( wfs, w_cx, nevt, 1 );
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 10 );
    retrieveWf( wfs, w_cx, nevt, 45 );
    retrieveWf( wfs, w_cx, nevt, 66 );
    retrieveWf( wfs, w_cx, nevt, 68 );
    retrieveWf( wfs, w_cx, nevt, 72 );
    retrieveWf( wfs, w_cx, nevt, 73 );
    retrieveWf( wfs, w_cx, nevt, 78 );
    retrieveWf( wfs, w_cx, nevt, 79 );
    retrieveWf( wfs, w_cx, nevt, 82 );
    retrieveWf( wfs, w_cx, nevt, 86 );
    retrieveWf( wfs, w_cx, nevt, 87 );
    retrieveWf( wfs, w_cx, nevt, 91 );
    retrieveWf( wfs, w_cx, nevt, 95 );
    retrieveWf( wfs, w_cx, nevt, 97 );
    retrieveWf( wfs, w_cx, nevt, 100 );
    retrieveWf( wfs, w_cx, nevt, 107 );
    retrieveWf( wfs, w_cx, nevt, 108 );
    retrieveWf( wfs, w_cx, nevt, 111 );
    retrieveWf( wfs, w_cx, nevt, 113 );
    retrieveWf( wfs, w_cx, nevt, 115 );
    retrieveWf( wfs, w_cx, nevt, 116 );
    retrieveWf( wfs, w_cx, nevt, 123 );
    retrieveWf( wfs, w_cx, nevt, 125 );
    retrieveWf( wfs, w_cx, nevt, 127 );
    retrieveWf( wfs, w_cx, nevt, 133 );
    retrieveWf( wfs, w_cx, nevt, 134 );
    retrieveWf( wfs, w_cx, nevt, 135 );
    retrieveWf( wfs, w_cx, nevt, 174 );
    retrieveWf( wfs, w_cx, nevt, 175 );
    retrieveWf( wfs, w_cx, nevt, 177 );
    retrieveWf( wfs, w_cx, nevt, 184 );
    retrieveWf( wfs, w_cx, nevt, 196 );
    retrieveWf( wfs, w_cx, nevt, 198 );
    retrieveWf( wfs, w_cx, nevt, 199 );
    retrieveWf( wfs, w_cx, nevt, 202 );
    retrieveWf( wfs, w_cx, nevt, 206 );
    retrieveWf( wfs, w_cx, nevt, 216 );
    retrieveWf( wfs, w_cx, nevt, 221 );
    retrieveWf( wfs, w_cx, nevt, 243 );
    retrieveWf( wfs, w_cx, nevt, 247 );
    retrieveWf( wfs, w_cx, nevt, 248 );
    retrieveWf( wfs, w_cx, nevt, 249 );
    retrieveWf( wfs, w_cx, nevt, 259 );
    retrieveWf( wfs, w_cx, nevt, 282 );
    retrieveWf( wfs, w_cx, nevt, 305 );
    retrieveWf( wfs, w_cx, nevt, 325 );
    retrieveWf( wfs, w_cx, nevt, 332 );
    retrieveWf( wfs, w_cx, nevt, 333 );
    retrieveWf( wfs, w_cx, nevt, 334 );
    retrieveWf( wfs, w_cx, nevt, 344 );
    retrieveWf( wfs, w_cx, nevt, 345 );
#endif
#endif

    // *** DIAGRAM 1801 OF 15495 ***
    // Wavefunction(s) for diagram number 1801
    // (none)
    // Amplitude(s) for diagram number 1801
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[282], w_fp[333], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[42] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1802 OF 15495 ***
    // Wavefunction(s) for diagram number 1802
    // (none)
    // Amplitude(s) for diagram number 1802
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[325], w_fp[4], w_fp[5], w_fp[282], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[42] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[325], w_fp[4], w_fp[5], w_fp[282], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[42] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[325], w_fp[4], w_fp[5], w_fp[282], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1803 OF 15495 ***
    // Wavefunction(s) for diagram number 1803
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[174], w_fp[325], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[348] );
    // Amplitude(s) for diagram number 1803
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[348], w_fp[45], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1804 OF 15495 ***
    // Wavefunction(s) for diagram number 1804
    // (none)
    // Amplitude(s) for diagram number 1804
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[174], w_fp[45], w_fp[333], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[28] += amp_sv[0];
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[42] -= amp_sv[0];

    // *** DIAGRAM 1805 OF 15495 ***
    // Wavefunction(s) for diagram number 1805
    // (none)
    // Amplitude(s) for diagram number 1805
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[175], w_fp[45], w_fp[325], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[28] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[42] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1806 OF 15495 ***
    // Wavefunction(s) for diagram number 1806
    // (none)
    // Amplitude(s) for diagram number 1806
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[348], w_fp[95], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1807 OF 15495 ***
    // Wavefunction(s) for diagram number 1807
    // (none)
    // Amplitude(s) for diagram number 1807
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[174], w_fp[95], w_fp[332], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[52] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[66] -= amp_sv[0];

    // *** DIAGRAM 1808 OF 15495 ***
    // Wavefunction(s) for diagram number 1808
    // (none)
    // Amplitude(s) for diagram number 1808
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[202], w_fp[95], w_fp[325], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1809 OF 15495 ***
    // Wavefunction(s) for diagram number 1809
    // (none)
    // Amplitude(s) for diagram number 1809
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[175], w_fp[10], w_fp[332], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[96] -= amp_sv[0];

    // *** DIAGRAM 1810 OF 15495 ***
    // Wavefunction(s) for diagram number 1810
    // (none)
    // Amplitude(s) for diagram number 1810
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[202], w_fp[10], w_fp[333], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];

    // *** DIAGRAM 1811 OF 15495 ***
    // Wavefunction(s) for diagram number 1811
    // (none)
    // Amplitude(s) for diagram number 1811
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[174], w_fp[345], w_fp[66], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];

    // *** DIAGRAM 1812 OF 15495 ***
    // Wavefunction(s) for diagram number 1812
    // (none)
    // Amplitude(s) for diagram number 1812
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[325], w_fp[66], w_fp[282], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1813 OF 15495 ***
    // Wavefunction(s) for diagram number 1813
    // (none)
    // Amplitude(s) for diagram number 1813
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[174], w_fp[91], w_fp[325], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];

    // *** DIAGRAM 1814 OF 15495 ***
    // Wavefunction(s) for diagram number 1814
    // (none)
    // Amplitude(s) for diagram number 1814
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[221], w_fp[345], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[101] += amp_sv[0];

    // *** DIAGRAM 1815 OF 15495 ***
    // Wavefunction(s) for diagram number 1815
    // (none)
    // Amplitude(s) for diagram number 1815
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[345], w_fp[68], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1816 OF 15495 ***
    // Wavefunction(s) for diagram number 1816
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[325], w_fp[66], COUPs[0], 1.0, 0., 0., w_fp[349] );
    // Amplitude(s) for diagram number 1816
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[259], w_fp[349], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];

    // *** DIAGRAM 1817 OF 15495 ***
    // Wavefunction(s) for diagram number 1817
    // (none)
    // Amplitude(s) for diagram number 1817
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[259], w_fp[334], w_fp[66], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[33] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];

    // *** DIAGRAM 1818 OF 15495 ***
    // Wavefunction(s) for diagram number 1818
    // (none)
    // Amplitude(s) for diagram number 1818
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[259], w_fp[325], w_fp[68], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[32] += amp_sv[0];
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[101] += amp_sv[0];

    // *** DIAGRAM 1819 OF 15495 ***
    // Wavefunction(s) for diagram number 1819
    // (none)
    // Amplitude(s) for diagram number 1819
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[325], w_fp[66], w_fp[6], w_fp[259], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[325], w_fp[66], w_fp[6], w_fp[259], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[33] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[325], w_fp[66], w_fp[6], w_fp[259], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[33] += amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[87] += amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];

    // *** DIAGRAM 1820 OF 15495 ***
    // Wavefunction(s) for diagram number 1820
    // (none)
    // Amplitude(s) for diagram number 1820
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[344], w_fp[91], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] += amp_sv[0];
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[57] += amp_sv[0];

    // *** DIAGRAM 1821 OF 15495 ***
    // Wavefunction(s) for diagram number 1821
    // (none)
    // Amplitude(s) for diagram number 1821
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[91], w_fp[334], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[33] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1822 OF 15495 ***
    // Wavefunction(s) for diagram number 1822
    // (none)
    // Amplitude(s) for diagram number 1822
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[344], w_fp[177], w_fp[66], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[80] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[87] += amp_sv[0];

    // *** DIAGRAM 1823 OF 15495 ***
    // Wavefunction(s) for diagram number 1823
    // (none)
    // Amplitude(s) for diagram number 1823
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[177], w_fp[349], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[76] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[77] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[90] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1824 OF 15495 ***
    // Wavefunction(s) for diagram number 1824
    // (none)
    // Amplitude(s) for diagram number 1824
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[221], w_fp[177], w_fp[325], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[76] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];

    // *** DIAGRAM 1825 OF 15495 ***
    // Wavefunction(s) for diagram number 1825
    // (none)
    // Amplitude(s) for diagram number 1825
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[344], w_fp[10], w_fp[68], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[33] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1826 OF 15495 ***
    // Wavefunction(s) for diagram number 1826
    // (none)
    // Amplitude(s) for diagram number 1826
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[221], w_fp[10], w_fp[334], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[76] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[77] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[90] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1827 OF 15495 ***
    // Wavefunction(s) for diagram number 1827
    // (none)
    // Amplitude(s) for diagram number 1827
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[206], w_fp[345], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[99] += amp_sv[0];

    // *** DIAGRAM 1828 OF 15495 ***
    // Wavefunction(s) for diagram number 1828
    // (none)
    // Amplitude(s) for diagram number 1828
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[345], w_fp[87], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1829 OF 15495 ***
    // Wavefunction(s) for diagram number 1829
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[325], w_fp[86], COUPs[0], 1.0, 0., 0., w_fp[350] );
    // Amplitude(s) for diagram number 1829
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[259], w_fp[350], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[62] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[100] -= amp_sv[0];

    // *** DIAGRAM 1830 OF 15495 ***
    // Wavefunction(s) for diagram number 1830
    // (none)
    // Amplitude(s) for diagram number 1830
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[259], w_fp[333], w_fp[86], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[39] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];

    // *** DIAGRAM 1831 OF 15495 ***
    // Wavefunction(s) for diagram number 1831
    // (none)
    // Amplitude(s) for diagram number 1831
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[259], w_fp[325], w_fp[87], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[38] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[62] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[99] += amp_sv[0];
    jamp_sv[100] -= amp_sv[0];

    // *** DIAGRAM 1832 OF 15495 ***
    // Wavefunction(s) for diagram number 1832
    // (none)
    // Amplitude(s) for diagram number 1832
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[325], w_fp[86], w_fp[5], w_fp[259], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[62] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[325], w_fp[86], w_fp[5], w_fp[259], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[39] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[325], w_fp[86], w_fp[5], w_fp[259], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[39] += amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[100] += amp_sv[0];

    // *** DIAGRAM 1833 OF 15495 ***
    // Wavefunction(s) for diagram number 1833
    // (none)
    // Amplitude(s) for diagram number 1833
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[344], w_fp[111], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[38] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[81] += amp_sv[0];

    // *** DIAGRAM 1834 OF 15495 ***
    // Wavefunction(s) for diagram number 1834
    // (none)
    // Amplitude(s) for diagram number 1834
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[111], w_fp[333], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[39] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[40] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1835 OF 15495 ***
    // Wavefunction(s) for diagram number 1835
    // (none)
    // Amplitude(s) for diagram number 1835
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[344], w_fp[95], w_fp[86], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[56] += amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[63] += amp_sv[0];

    // *** DIAGRAM 1836 OF 15495 ***
    // Wavefunction(s) for diagram number 1836
    // (none)
    // Amplitude(s) for diagram number 1836
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[95], w_fp[350], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[52] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1837 OF 15495 ***
    // Wavefunction(s) for diagram number 1837
    // (none)
    // Amplitude(s) for diagram number 1837
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[206], w_fp[95], w_fp[325], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[52] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[67] += amp_sv[0];

    // *** DIAGRAM 1838 OF 15495 ***
    // Wavefunction(s) for diagram number 1838
    // (none)
    // Amplitude(s) for diagram number 1838
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[344], w_fp[10], w_fp[87], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[38] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[39] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1839 OF 15495 ***
    // Wavefunction(s) for diagram number 1839
    // (none)
    // Amplitude(s) for diagram number 1839
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[206], w_fp[10], w_fp[333], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1840 OF 15495 ***
    // Wavefunction(s) for diagram number 1840
    // (none)
    // Amplitude(s) for diagram number 1840
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[184], w_fp[345], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[97] += amp_sv[0];

    // *** DIAGRAM 1841 OF 15495 ***
    // Wavefunction(s) for diagram number 1841
    // (none)
    // Amplitude(s) for diagram number 1841
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[345], w_fp[115], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1842 OF 15495 ***
    // Wavefunction(s) for diagram number 1842
    // (none)
    // Amplitude(s) for diagram number 1842
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[259], w_fp[332], w_fp[113], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[19] += amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[85] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];

    // *** DIAGRAM 1843 OF 15495 ***
    // Wavefunction(s) for diagram number 1843
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[325], w_fp[113], COUPs[0], 1.0, 0., 0., w_fp[351] );
    // Amplitude(s) for diagram number 1843
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[259], w_fp[351], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[33] += amp_sv[0];
    jamp_sv[38] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[85] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[99] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];

    // *** DIAGRAM 1844 OF 15495 ***
    // Wavefunction(s) for diagram number 1844
    // (none)
    // Amplitude(s) for diagram number 1844
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[259], w_fp[325], w_fp[115], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[19] += amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[32] += amp_sv[0];
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[39] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[101] += amp_sv[0];

    // *** DIAGRAM 1845 OF 15495 ***
    // Wavefunction(s) for diagram number 1845
    // (none)
    // Amplitude(s) for diagram number 1845
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[325], w_fp[4], w_fp[113], w_fp[259], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[19] += amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[85] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[325], w_fp[4], w_fp[113], w_fp[259], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[33] += amp_sv[0];
    jamp_sv[38] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[85] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[99] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[325], w_fp[4], w_fp[113], w_fp[259], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[33] += amp_sv[0];
    jamp_sv[38] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[62] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[87] += amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[99] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];

    // *** DIAGRAM 1846 OF 15495 ***
    // Wavefunction(s) for diagram number 1846
    // (none)
    // Amplitude(s) for diagram number 1846
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[344], w_fp[45], w_fp[113], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] += amp_sv[0];
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[39] += amp_sv[0];

    // *** DIAGRAM 1847 OF 15495 ***
    // Wavefunction(s) for diagram number 1847
    // (none)
    // Amplitude(s) for diagram number 1847
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[45], w_fp[351], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[28] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[33] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[39] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[42] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[43] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1848 OF 15495 ***
    // Wavefunction(s) for diagram number 1848
    // (none)
    // Amplitude(s) for diagram number 1848
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[184], w_fp[45], w_fp[325], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[28] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[43] += amp_sv[0];

    // *** DIAGRAM 1849 OF 15495 ***
    // Wavefunction(s) for diagram number 1849
    // (none)
    // Amplitude(s) for diagram number 1849
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[344], w_fp[127], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[62] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[87] += amp_sv[0];

    // *** DIAGRAM 1850 OF 15495 ***
    // Wavefunction(s) for diagram number 1850
    // (none)
    // Amplitude(s) for diagram number 1850
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[127], w_fp[332], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1851 OF 15495 ***
    // Wavefunction(s) for diagram number 1851
    // (none)
    // Amplitude(s) for diagram number 1851
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[344], w_fp[10], w_fp[115], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[33] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[39] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1852 OF 15495 ***
    // Wavefunction(s) for diagram number 1852
    // (none)
    // Amplitude(s) for diagram number 1852
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[184], w_fp[10], w_fp[332], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[42] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[43] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1853 OF 15495 ***
    // Wavefunction(s) for diagram number 1853
    // (none)
    // Amplitude(s) for diagram number 1853
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[345], w_fp[133], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[345], w_fp[134], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[345], w_fp[135], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1854 OF 15495 ***
    // Wavefunction(s) for diagram number 1854
    // (none)
    // Amplitude(s) for diagram number 1854
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[325], w_fp[133], w_fp[259], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[19] += amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[32] += amp_sv[0];
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[39] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[101] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[325], w_fp[134], w_fp[259], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[39] += amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[325], w_fp[135], w_fp[259], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[33] += amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[87] += amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];

    // *** DIAGRAM 1855 OF 15495 ***
    // Wavefunction(s) for diagram number 1855
    // (none)
    // Amplitude(s) for diagram number 1855
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[247], w_fp[325], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[33] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[39] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[248], w_fp[325], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[38] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[39] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[249], w_fp[325], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[33] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1856 OF 15495 ***
    // Wavefunction(s) for diagram number 1856
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[10], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[345] );
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[345], w_fp[6], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[352] );
    // Amplitude(s) for diagram number 1856
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[216], w_fp[352], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] -= amp_sv[0];

    // *** DIAGRAM 1857 OF 15495 ***
    // Wavefunction(s) for diagram number 1857
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[345], w_fp[7], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[353] );
    // Amplitude(s) for diagram number 1857
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[216], w_fp[353], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[23] -= amp_sv[0];

    // *** DIAGRAM 1858 OF 15495 ***
    // Wavefunction(s) for diagram number 1858
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[345], w_fp[5], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[354] );
    // Amplitude(s) for diagram number 1858
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[198], w_fp[354], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] -= amp_sv[0];

    // *** DIAGRAM 1859 OF 15495 ***
    // Wavefunction(s) for diagram number 1859
    // (none)
    // Amplitude(s) for diagram number 1859
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[198], w_fp[353], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] -= amp_sv[0];

    // *** DIAGRAM 1860 OF 15495 ***
    // Wavefunction(s) for diagram number 1860
    // (none)
    // Amplitude(s) for diagram number 1860
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[199], w_fp[354], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= amp_sv[0];

    // *** DIAGRAM 1861 OF 15495 ***
    // Wavefunction(s) for diagram number 1861
    // (none)
    // Amplitude(s) for diagram number 1861
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[199], w_fp[352], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] -= amp_sv[0];

    // *** DIAGRAM 1862 OF 15495 ***
    // Wavefunction(s) for diagram number 1862
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[196], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[355] );
    // Amplitude(s) for diagram number 1862
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[355], w_fp[79], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= amp_sv[0];

    // *** DIAGRAM 1863 OF 15495 ***
    // Wavefunction(s) for diagram number 1863
    // (none)
    // Amplitude(s) for diagram number 1863
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[355], w_fp[78], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[70] -= amp_sv[0];

    // *** DIAGRAM 1864 OF 15495 ***
    // Wavefunction(s) for diagram number 1864
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[95], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[356] );
    // Amplitude(s) for diagram number 1864
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[198], w_fp[356], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[53] -= amp_sv[0];

    // *** DIAGRAM 1865 OF 15495 ***
    // Wavefunction(s) for diagram number 1865
    // (none)
    // Amplitude(s) for diagram number 1865
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[198], w_fp[78], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[67] -= amp_sv[0];

    // *** DIAGRAM 1866 OF 15495 ***
    // Wavefunction(s) for diagram number 1866
    // (none)
    // Amplitude(s) for diagram number 1866
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[199], w_fp[356], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] -= amp_sv[0];

    // *** DIAGRAM 1867 OF 15495 ***
    // Wavefunction(s) for diagram number 1867
    // (none)
    // Amplitude(s) for diagram number 1867
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[199], w_fp[79], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[61] -= amp_sv[0];

    // *** DIAGRAM 1868 OF 15495 ***
    // Wavefunction(s) for diagram number 1868
    // (none)
    // Amplitude(s) for diagram number 1868
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[355], w_fp[72], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[88] -= amp_sv[0];

    // *** DIAGRAM 1869 OF 15495 ***
    // Wavefunction(s) for diagram number 1869
    // (none)
    // Amplitude(s) for diagram number 1869
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[355], w_fp[73], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[94] -= amp_sv[0];

    // *** DIAGRAM 1870 OF 15495 ***
    // Wavefunction(s) for diagram number 1870
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[177], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[357] );
    // Amplitude(s) for diagram number 1870
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[216], w_fp[357], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[77] -= amp_sv[0];

    // *** DIAGRAM 1871 OF 15495 ***
    // Wavefunction(s) for diagram number 1871
    // (none)
    // Amplitude(s) for diagram number 1871
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[216], w_fp[73], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[91] -= amp_sv[0];

    // *** DIAGRAM 1872 OF 15495 ***
    // Wavefunction(s) for diagram number 1872
    // (none)
    // Amplitude(s) for diagram number 1872
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[199], w_fp[357], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[75] -= amp_sv[0];

    // *** DIAGRAM 1873 OF 15495 ***
    // Wavefunction(s) for diagram number 1873
    // (none)
    // Amplitude(s) for diagram number 1873
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[199], w_fp[72], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[85] -= amp_sv[0];

    // *** DIAGRAM 1874 OF 15495 ***
    // Wavefunction(s) for diagram number 1874
    // (none)
    // Amplitude(s) for diagram number 1874
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[355], w_fp[107], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[112] -= amp_sv[0];

    // *** DIAGRAM 1875 OF 15495 ***
    // Wavefunction(s) for diagram number 1875
    // (none)
    // Amplitude(s) for diagram number 1875
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[355], w_fp[82], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[118] -= amp_sv[0];

    // *** DIAGRAM 1876 OF 15495 ***
    // Wavefunction(s) for diagram number 1876
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[108], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[358] );
    // Amplitude(s) for diagram number 1876
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[216], w_fp[358], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[101] -= amp_sv[0];

    // *** DIAGRAM 1877 OF 15495 ***
    // Wavefunction(s) for diagram number 1877
    // (none)
    // Amplitude(s) for diagram number 1877
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[216], w_fp[82], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[115] -= amp_sv[0];

    // *** DIAGRAM 1878 OF 15495 ***
    // Wavefunction(s) for diagram number 1878
    // (none)
    // Amplitude(s) for diagram number 1878
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[198], w_fp[358], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] -= amp_sv[0];

    // *** DIAGRAM 1879 OF 15495 ***
    // Wavefunction(s) for diagram number 1879
    // (none)
    // Amplitude(s) for diagram number 1879
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[198], w_fp[107], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[109] -= amp_sv[0];

    // *** DIAGRAM 1880 OF 15495 ***
    // Wavefunction(s) for diagram number 1880
    // (none)
    // Amplitude(s) for diagram number 1880
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[243], w_fp[345], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1881 OF 15495 ***
    // Wavefunction(s) for diagram number 1881
    // (none)
    // Amplitude(s) for diagram number 1881
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[199], w_fp[345], w_fp[113], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1882 OF 15495 ***
    // Wavefunction(s) for diagram number 1882
    // (none)
    // Amplitude(s) for diagram number 1882
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[345], w_fp[116], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];

    // *** DIAGRAM 1883 OF 15495 ***
    // Wavefunction(s) for diagram number 1883
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[113], COUPs[0], 1.0, 0., 0., w_fp[359] );
    // Amplitude(s) for diagram number 1883
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[305], w_fp[359], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1884 OF 15495 ***
    // Wavefunction(s) for diagram number 1884
    // (none)
    // Amplitude(s) for diagram number 1884
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[305], w_fp[1], w_fp[116], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1885 OF 15495 ***
    // Wavefunction(s) for diagram number 1885
    // (none)
    // Amplitude(s) for diagram number 1885
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[113], w_fp[7], w_fp[305], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[113], w_fp[7], w_fp[305], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[64] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[113], w_fp[7], w_fp[305], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[64] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1886 OF 15495 ***
    // Wavefunction(s) for diagram number 1886
    // (none)
    // Amplitude(s) for diagram number 1886
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[355], w_fp[127], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1887 OF 15495 ***
    // Wavefunction(s) for diagram number 1887
    // (none)
    // Amplitude(s) for diagram number 1887
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[199], w_fp[127], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[61] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1888 OF 15495 ***
    // Wavefunction(s) for diagram number 1888
    // (none)
    // Amplitude(s) for diagram number 1888
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[355], w_fp[108], w_fp[113], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1889 OF 15495 ***
    // Wavefunction(s) for diagram number 1889
    // (none)
    // Amplitude(s) for diagram number 1889
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[108], w_fp[359], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];

    // *** DIAGRAM 1890 OF 15495 ***
    // Wavefunction(s) for diagram number 1890
    // (none)
    // Amplitude(s) for diagram number 1890
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[243], w_fp[108], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1891 OF 15495 ***
    // Wavefunction(s) for diagram number 1891
    // (none)
    // Amplitude(s) for diagram number 1891
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[355], w_fp[10], w_fp[116], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] += amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];

    // *** DIAGRAM 1892 OF 15495 ***
    // Wavefunction(s) for diagram number 1892
    // (none)
    // Amplitude(s) for diagram number 1892
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[199], w_fp[10], w_fp[359], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[61] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];

    // *** DIAGRAM 1893 OF 15495 ***
    // Wavefunction(s) for diagram number 1893
    // (none)
    // Amplitude(s) for diagram number 1893
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[123], w_fp[345], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1894 OF 15495 ***
    // Wavefunction(s) for diagram number 1894
    // (none)
    // Amplitude(s) for diagram number 1894
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[198], w_fp[345], w_fp[100], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1895 OF 15495 ***
    // Wavefunction(s) for diagram number 1895
    // (none)
    // Amplitude(s) for diagram number 1895
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[345], w_fp[125], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];

    // *** DIAGRAM 1896 OF 15495 ***
    // Wavefunction(s) for diagram number 1896
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[100], COUPs[0], 1.0, 0., 0., w_fp[360] );
    // Amplitude(s) for diagram number 1896
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[305], w_fp[360], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[77] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1897 OF 15495 ***
    // Wavefunction(s) for diagram number 1897
    // (none)
    // Amplitude(s) for diagram number 1897
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[305], w_fp[1], w_fp[125], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1898 OF 15495 ***
    // Wavefunction(s) for diagram number 1898
    // (none)
    // Amplitude(s) for diagram number 1898
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[100], w_fp[6], w_fp[305], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[77] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[100], w_fp[6], w_fp[305], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[77] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[100], w_fp[6], w_fp[305], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1899 OF 15495 ***
    // Wavefunction(s) for diagram number 1899
    // (none)
    // Amplitude(s) for diagram number 1899
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[355], w_fp[97], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1900 OF 15495 ***
    // Wavefunction(s) for diagram number 1900
    // (none)
    // Amplitude(s) for diagram number 1900
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[198], w_fp[97], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE JAMPS ***
    // In CUDA, copy the local jamp to the output global-memory jamp
    constexpr int ihel0 = 0; // the allJamps buffer already points to a specific helicity _within a super-buffer for nGoodHel helicities_
    for( int icol = 0; icol < ncolor; icol++ )
      J_ACCESS::kernelAccessIcolIhelNhel( jamps, icol, ihel0, nGoodHel ) += jamp_sv[icol]; // update jamps
#else
    // In C++, copy the local jamp to the output array passed as function argument
    for( int icol = 0; icol < ncolor; icol++ )
      jamps[icol] += jamp_sv[icol]; // update jamps
#endif

#ifdef MGONGPUCPP_GPUIMPL
#ifndef MGONGPU_RDC_DIAGRAMS
    // *** STORE WAVEFUNCTIONS FOR NEXT DIAGRAM GROUPS ***
    //for( int iwf = 0; iwf < nwf; iwf++ ) storeWf( wfs, w_cx, nevt, iwf );
    storeWf( wfs, w_cx, nevt, 345 );
    storeWf( wfs, w_cx, nevt, 348 );
    storeWf( wfs, w_cx, nevt, 349 );
    storeWf( wfs, w_cx, nevt, 350 );
    storeWf( wfs, w_cx, nevt, 351 );
    storeWf( wfs, w_cx, nevt, 352 );
    storeWf( wfs, w_cx, nevt, 353 );
    storeWf( wfs, w_cx, nevt, 354 );
    storeWf( wfs, w_cx, nevt, 355 );
    storeWf( wfs, w_cx, nevt, 356 );
    storeWf( wfs, w_cx, nevt, 357 );
    storeWf( wfs, w_cx, nevt, 358 );
    storeWf( wfs, w_cx, nevt, 359 );
    storeWf( wfs, w_cx, nevt, 360 );
#endif
#endif
  }

  //--------------------------------------------------------------------------
}
