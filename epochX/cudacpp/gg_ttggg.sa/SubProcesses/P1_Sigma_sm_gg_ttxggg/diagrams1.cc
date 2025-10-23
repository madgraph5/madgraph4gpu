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

    // *** DIAGRAM 1 OF 1240 ***
    // Wavefunction(s) for diagram number 1
    vxxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][0], -1, w_fp[0], 0 );
    vxxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][1], -1, w_fp[1], 1 );
    oxxxxx<M_ACCESS, W_ACCESS>( momenta, cIPD[0], cHel[ihel][2], +1, w_fp[2], 2 );
    ixxxxx<M_ACCESS, W_ACCESS>( momenta, cIPD[0], cHel[ihel][3], -1, w_fp[3], 3 );
    vxxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][4], +1, w_fp[4], 4 );
    vxxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][5], +1, w_fp[5], 5 );
    vxxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][6], +1, w_fp[6], 6 );
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], COUPs[0], 1.0, 0., 0., w_fp[7] );
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[2], COUPs[1], 1.0, 0., 0., w_fp[8] );
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[7], w_fp[4], COUPs[0], 1.0, 0., 0., w_fp[9] );
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[8], w_fp[5], COUPs[0], 1.0, 0., 0., w_fp[10] );
    // Amplitude(s) for diagram number 1
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[10], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[25] += amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];

    // *** DIAGRAM 2 OF 1240 ***
    // Wavefunction(s) for diagram number 2
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[8], w_fp[6], COUPs[0], 1.0, 0., 0., w_fp[11] );
    // Amplitude(s) for diagram number 2
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[11], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];

    // *** DIAGRAM 3 OF 1240 ***
    // Wavefunction(s) for diagram number 3
    // (none)
    // Amplitude(s) for diagram number 3
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[5], w_fp[6], w_fp[9], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[25] += amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[5], w_fp[6], w_fp[9], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[5], w_fp[6], w_fp[9], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[1] += amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[90] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];

    // *** DIAGRAM 4 OF 1240 ***
    // Wavefunction(s) for diagram number 4
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[7], w_fp[5], COUPs[0], 1.0, 0., 0., w_fp[12] );
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[8], w_fp[4], COUPs[0], 1.0, 0., 0., w_fp[13] );
    // Amplitude(s) for diagram number 4
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[13], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[27] += amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[62] += amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[73] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];

    // *** DIAGRAM 5 OF 1240 ***
    // Wavefunction(s) for diagram number 5
    // (none)
    // Amplitude(s) for diagram number 5
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[11], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[62] += amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];

    // *** DIAGRAM 6 OF 1240 ***
    // Wavefunction(s) for diagram number 6
    // (none)
    // Amplitude(s) for diagram number 6
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[4], w_fp[6], w_fp[12], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[27] += amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[62] += amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[73] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[4], w_fp[6], w_fp[12], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[62] += amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[4], w_fp[6], w_fp[12], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[3] += amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];

    // *** DIAGRAM 7 OF 1240 ***
    // Wavefunction(s) for diagram number 7
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[7], w_fp[6], COUPs[0], 1.0, 0., 0., w_fp[14] );
    // Amplitude(s) for diagram number 7
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[14], w_fp[13], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[62] += amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[73] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];

    // *** DIAGRAM 8 OF 1240 ***
    // Wavefunction(s) for diagram number 8
    // (none)
    // Amplitude(s) for diagram number 8
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[14], w_fp[10], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[73] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];

    // *** DIAGRAM 9 OF 1240 ***
    // Wavefunction(s) for diagram number 9
    // (none)
    // Amplitude(s) for diagram number 9
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[4], w_fp[5], w_fp[14], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[62] += amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[73] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[4], w_fp[5], w_fp[14], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[73] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[4], w_fp[5], w_fp[14], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[60] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];

    // *** DIAGRAM 10 OF 1240 ***
    // Wavefunction(s) for diagram number 10
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[7], w_fp[4], w_fp[5], COUPs[2], 1.0, 0., 0., w_fp[15] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[7], w_fp[4], w_fp[5], COUPs[2], 1.0, 0., 0., w_fp[16] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[7], w_fp[4], w_fp[5], COUPs[2], 1.0, 0., 0., w_fp[17] );
    // Amplitude(s) for diagram number 10
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[6], w_fp[15], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[62] += amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[6], w_fp[16], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[62] += amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[6], w_fp[17], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];

    // *** DIAGRAM 11 OF 1240 ***
    // Wavefunction(s) for diagram number 11
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[7], w_fp[4], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[18] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[7], w_fp[4], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[19] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[7], w_fp[4], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[20] );
    // Amplitude(s) for diagram number 11
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[5], w_fp[18], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[73] += amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[5], w_fp[19], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[73] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[5], w_fp[20], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[25] += amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];

    // *** DIAGRAM 12 OF 1240 ***
    // Wavefunction(s) for diagram number 12
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[7], w_fp[5], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[21] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[7], w_fp[5], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[22] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[7], w_fp[5], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[23] );
    // Amplitude(s) for diagram number 12
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[4], w_fp[21], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[4], w_fp[22], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[62] += amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[73] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[4], w_fp[23], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[27] += amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[62] += amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[73] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];

    // *** DIAGRAM 13 OF 1240 ***
    // Wavefunction(s) for diagram number 13
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[4], w_fp[5], COUPs[0], 1.0, 0., 0., w_fp[24] );
    // Amplitude(s) for diagram number 13
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[7], w_fp[8], w_fp[24], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[7], w_fp[8], w_fp[24], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[62] += amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[7], w_fp[8], w_fp[24], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[62] += amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];

    // *** DIAGRAM 14 OF 1240 ***
    // Wavefunction(s) for diagram number 14
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[7], w_fp[8], COUPs[0], 1.0, 0., 0., w_fp[25] );
    // Amplitude(s) for diagram number 14
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[24], w_fp[6], w_fp[25], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];

    // *** DIAGRAM 15 OF 1240 ***
    // Wavefunction(s) for diagram number 15
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[7], w_fp[24], COUPs[0], 1.0, 0., 0., w_fp[26] );
    // Amplitude(s) for diagram number 15
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[6], w_fp[26], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[62] += amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];

    // *** DIAGRAM 16 OF 1240 ***
    // Wavefunction(s) for diagram number 16
    // (none)
    // Amplitude(s) for diagram number 16
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[24], w_fp[14], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[62] += amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];

    // *** DIAGRAM 17 OF 1240 ***
    // Wavefunction(s) for diagram number 17
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[4], w_fp[6], COUPs[0], 1.0, 0., 0., w_fp[27] );
    // Amplitude(s) for diagram number 17
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[7], w_fp[8], w_fp[27], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[3] += amp_sv[0];
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[7], w_fp[8], w_fp[27], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[73] += amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[7], w_fp[8], w_fp[27], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[27] += amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[73] += amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];

    // *** DIAGRAM 18 OF 1240 ***
    // Wavefunction(s) for diagram number 18
    // (none)
    // Amplitude(s) for diagram number 18
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[27], w_fp[5], w_fp[25], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[3] += amp_sv[0];
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];

    // *** DIAGRAM 19 OF 1240 ***
    // Wavefunction(s) for diagram number 19
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[7], w_fp[27], COUPs[0], 1.0, 0., 0., w_fp[28] );
    // Amplitude(s) for diagram number 19
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[5], w_fp[28], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[73] += amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];

    // *** DIAGRAM 20 OF 1240 ***
    // Wavefunction(s) for diagram number 20
    // (none)
    // Amplitude(s) for diagram number 20
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[27], w_fp[12], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[27] += amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[73] += amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];

    // *** DIAGRAM 21 OF 1240 ***
    // Wavefunction(s) for diagram number 21
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[5], w_fp[6], COUPs[0], 1.0, 0., 0., w_fp[29] );
    // Amplitude(s) for diagram number 21
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[7], w_fp[8], w_fp[4], w_fp[29], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[25] += amp_sv[0];
    jamp_sv[27] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[7], w_fp[8], w_fp[4], w_fp[29], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[25] += amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[7], w_fp[8], w_fp[4], w_fp[29], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];

    // *** DIAGRAM 22 OF 1240 ***
    // Wavefunction(s) for diagram number 22
    // (none)
    // Amplitude(s) for diagram number 22
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[4], w_fp[29], w_fp[25], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[25] += amp_sv[0];
    jamp_sv[27] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];

    // *** DIAGRAM 23 OF 1240 ***
    // Wavefunction(s) for diagram number 23
    // (none)
    // Amplitude(s) for diagram number 23
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[29], w_fp[9], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[25] += amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];

    // *** DIAGRAM 24 OF 1240 ***
    // Wavefunction(s) for diagram number 24
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[7], w_fp[29], COUPs[0], 1.0, 0., 0., w_fp[25] );
    // Amplitude(s) for diagram number 24
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[4], w_fp[25], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];

    // *** DIAGRAM 25 OF 1240 ***
    // Wavefunction(s) for diagram number 25
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[4], w_fp[5], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[30] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[4], w_fp[5], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[31] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[4], w_fp[5], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[32] );
    // Amplitude(s) for diagram number 25
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[7], w_fp[8], w_fp[30], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[25] += amp_sv[0];
    jamp_sv[27] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[7], w_fp[8], w_fp[31], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[2] += amp_sv[0];
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[4] += amp_sv[0];
    jamp_sv[25] += amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[27] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[7], w_fp[8], w_fp[32], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[2] += amp_sv[0];
    jamp_sv[4] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];

    // *** DIAGRAM 26 OF 1240 ***
    // Wavefunction(s) for diagram number 26
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[4], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[33] );
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[7], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[34] );
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[33], w_fp[5], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[35] );
    // Amplitude(s) for diagram number 26
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[34], w_fp[35], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 27 OF 1240 ***
    // Wavefunction(s) for diagram number 27
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[33], w_fp[6], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[36] );
    // Amplitude(s) for diagram number 27
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[34], w_fp[36], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 28 OF 1240 ***
    // Wavefunction(s) for diagram number 28
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[33], COUPs[1], 1.0, 0., 0., w_fp[37] );
    // Amplitude(s) for diagram number 28
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[37], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 29 OF 1240 ***
    // Wavefunction(s) for diagram number 29
    // (none)
    // Amplitude(s) for diagram number 29
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[36], w_fp[12], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];

    // *** DIAGRAM 30 OF 1240 ***
    // Wavefunction(s) for diagram number 30
    // (none)
    // Amplitude(s) for diagram number 30
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[14], w_fp[37], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 31 OF 1240 ***
    // Wavefunction(s) for diagram number 31
    // (none)
    // Amplitude(s) for diagram number 31
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[35], w_fp[14], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[60] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];

    // *** DIAGRAM 32 OF 1240 ***
    // Wavefunction(s) for diagram number 32
    // (none)
    // Amplitude(s) for diagram number 32
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[33], w_fp[21], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[33], w_fp[22], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[33], w_fp[23], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 33 OF 1240 ***
    // Wavefunction(s) for diagram number 33
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[5], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[38] );
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[33], w_fp[7], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[39] );
    // Amplitude(s) for diagram number 33
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[39], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 34 OF 1240 ***
    // Wavefunction(s) for diagram number 34
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[38], w_fp[7], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[40] );
    // Amplitude(s) for diagram number 34
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[40], w_fp[33], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 35 OF 1240 ***
    // Wavefunction(s) for diagram number 35
    // (none)
    // Amplitude(s) for diagram number 35
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[33], w_fp[14], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];

    // *** DIAGRAM 36 OF 1240 ***
    // Wavefunction(s) for diagram number 36
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[6], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[41] );
    // Amplitude(s) for diagram number 36
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[39], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 37 OF 1240 ***
    // Wavefunction(s) for diagram number 37
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[41], w_fp[7], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[42] );
    // Amplitude(s) for diagram number 37
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[42], w_fp[33], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 38 OF 1240 ***
    // Wavefunction(s) for diagram number 38
    // (none)
    // Amplitude(s) for diagram number 38
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[33], w_fp[12], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[62] += amp_sv[0];

    // *** DIAGRAM 39 OF 1240 ***
    // Wavefunction(s) for diagram number 39
    // (none)
    // Amplitude(s) for diagram number 39
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[39], w_fp[29], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] += amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];

    // *** DIAGRAM 40 OF 1240 ***
    // Wavefunction(s) for diagram number 40
    // (none)
    // Amplitude(s) for diagram number 40
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[34], w_fp[33], w_fp[29], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];

    // *** DIAGRAM 41 OF 1240 ***
    // Wavefunction(s) for diagram number 41
    // (none)
    // Amplitude(s) for diagram number 41
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[33], w_fp[25], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 42 OF 1240 ***
    // Wavefunction(s) for diagram number 42
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[5], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[39] );
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[39], w_fp[4], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[43] );
    // Amplitude(s) for diagram number 42
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[34], w_fp[43], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[88] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[89] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 43 OF 1240 ***
    // Wavefunction(s) for diagram number 43
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[39], w_fp[6], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[44] );
    // Amplitude(s) for diagram number 43
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[34], w_fp[44], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[95] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 44 OF 1240 ***
    // Wavefunction(s) for diagram number 44
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[39], COUPs[1], 1.0, 0., 0., w_fp[45] );
    // Amplitude(s) for diagram number 44
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[45], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[72] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[90] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[95] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 45 OF 1240 ***
    // Wavefunction(s) for diagram number 45
    // (none)
    // Amplitude(s) for diagram number 45
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[44], w_fp[9], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[90] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];

    // *** DIAGRAM 46 OF 1240 ***
    // Wavefunction(s) for diagram number 46
    // (none)
    // Amplitude(s) for diagram number 46
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[14], w_fp[45], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[73] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[89] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[90] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 47 OF 1240 ***
    // Wavefunction(s) for diagram number 47
    // (none)
    // Amplitude(s) for diagram number 47
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[43], w_fp[14], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[84] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];

    // *** DIAGRAM 48 OF 1240 ***
    // Wavefunction(s) for diagram number 48
    // (none)
    // Amplitude(s) for diagram number 48
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[39], w_fp[18], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[72] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[73] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[89] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[95] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[39], w_fp[19], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[73] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[89] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[90] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[39], w_fp[20], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[72] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[90] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[95] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 49 OF 1240 ***
    // Wavefunction(s) for diagram number 49
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[4], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[46] );
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[39], w_fp[7], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[47] );
    // Amplitude(s) for diagram number 49
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[47], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[73] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 50 OF 1240 ***
    // Wavefunction(s) for diagram number 50
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[46], w_fp[7], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[48] );
    // Amplitude(s) for diagram number 50
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[48], w_fp[39], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[90] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 51 OF 1240 ***
    // Wavefunction(s) for diagram number 51
    // (none)
    // Amplitude(s) for diagram number 51
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[39], w_fp[14], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[73] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];

    // *** DIAGRAM 52 OF 1240 ***
    // Wavefunction(s) for diagram number 52
    // (none)
    // Amplitude(s) for diagram number 52
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[47], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[72] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 53 OF 1240 ***
    // Wavefunction(s) for diagram number 53
    // (none)
    // Amplitude(s) for diagram number 53
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[42], w_fp[39], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[84] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 54 OF 1240 ***
    // Wavefunction(s) for diagram number 54
    // (none)
    // Amplitude(s) for diagram number 54
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[39], w_fp[9], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[72] += amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];

    // *** DIAGRAM 55 OF 1240 ***
    // Wavefunction(s) for diagram number 55
    // (none)
    // Amplitude(s) for diagram number 55
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[47], w_fp[27], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[72] += amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];

    // *** DIAGRAM 56 OF 1240 ***
    // Wavefunction(s) for diagram number 56
    // (none)
    // Amplitude(s) for diagram number 56
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[34], w_fp[39], w_fp[27], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[88] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];

    // *** DIAGRAM 57 OF 1240 ***
    // Wavefunction(s) for diagram number 57
    // (none)
    // Amplitude(s) for diagram number 57
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[39], w_fp[28], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[72] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[73] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[89] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[95] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 58 OF 1240 ***
    // Wavefunction(s) for diagram number 58
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[6], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[47] );
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[47], w_fp[4], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[49] );
    // Amplitude(s) for diagram number 58
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[34], w_fp[49], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 59 OF 1240 ***
    // Wavefunction(s) for diagram number 59
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[47], w_fp[5], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[50] );
    // Amplitude(s) for diagram number 59
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[34], w_fp[50], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 60 OF 1240 ***
    // Wavefunction(s) for diagram number 60
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[47], COUPs[1], 1.0, 0., 0., w_fp[51] );
    // Amplitude(s) for diagram number 60
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[51], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 61 OF 1240 ***
    // Wavefunction(s) for diagram number 61
    // (none)
    // Amplitude(s) for diagram number 61
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[50], w_fp[9], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[114] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];

    // *** DIAGRAM 62 OF 1240 ***
    // Wavefunction(s) for diagram number 62
    // (none)
    // Amplitude(s) for diagram number 62
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[51], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 63 OF 1240 ***
    // Wavefunction(s) for diagram number 63
    // (none)
    // Amplitude(s) for diagram number 63
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[49], w_fp[12], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[108] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];

    // *** DIAGRAM 64 OF 1240 ***
    // Wavefunction(s) for diagram number 64
    // (none)
    // Amplitude(s) for diagram number 64
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[47], w_fp[15], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[47], w_fp[16], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[47], w_fp[17], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 65 OF 1240 ***
    // Wavefunction(s) for diagram number 65
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[47], w_fp[7], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[52] );
    // Amplitude(s) for diagram number 65
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[52], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 66 OF 1240 ***
    // Wavefunction(s) for diagram number 66
    // (none)
    // Amplitude(s) for diagram number 66
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[48], w_fp[47], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[114] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 67 OF 1240 ***
    // Wavefunction(s) for diagram number 67
    // (none)
    // Amplitude(s) for diagram number 67
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[47], w_fp[12], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];

    // *** DIAGRAM 68 OF 1240 ***
    // Wavefunction(s) for diagram number 68
    // (none)
    // Amplitude(s) for diagram number 68
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[52], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 69 OF 1240 ***
    // Wavefunction(s) for diagram number 69
    // (none)
    // Amplitude(s) for diagram number 69
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[40], w_fp[47], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[108] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 70 OF 1240 ***
    // Wavefunction(s) for diagram number 70
    // (none)
    // Amplitude(s) for diagram number 70
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[47], w_fp[9], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];

    // *** DIAGRAM 71 OF 1240 ***
    // Wavefunction(s) for diagram number 71
    // (none)
    // Amplitude(s) for diagram number 71
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[52], w_fp[24], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];

    // *** DIAGRAM 72 OF 1240 ***
    // Wavefunction(s) for diagram number 72
    // (none)
    // Amplitude(s) for diagram number 72
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[34], w_fp[47], w_fp[24], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[112] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];

    // *** DIAGRAM 73 OF 1240 ***
    // Wavefunction(s) for diagram number 73
    // (none)
    // Amplitude(s) for diagram number 73
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[47], w_fp[26], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 74 OF 1240 ***
    // Wavefunction(s) for diagram number 74
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[7], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[52] );
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[46], w_fp[5], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[7] );
    // Amplitude(s) for diagram number 74
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[7], w_fp[52], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 75 OF 1240 ***
    // Wavefunction(s) for diagram number 75
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[46], w_fp[6], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[53] );
    // Amplitude(s) for diagram number 75
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[53], w_fp[52], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 76 OF 1240 ***
    // Wavefunction(s) for diagram number 76
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[46], w_fp[2], COUPs[1], 1.0, 0., 0., w_fp[54] );
    // Amplitude(s) for diagram number 76
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[54], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[73] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 77 OF 1240 ***
    // Wavefunction(s) for diagram number 77
    // (none)
    // Amplitude(s) for diagram number 77
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[53], w_fp[2], w_fp[12], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];

    // *** DIAGRAM 78 OF 1240 ***
    // Wavefunction(s) for diagram number 78
    // (none)
    // Amplitude(s) for diagram number 78
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[14], w_fp[54], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[73] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[90] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 79 OF 1240 ***
    // Wavefunction(s) for diagram number 79
    // (none)
    // Amplitude(s) for diagram number 79
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[7], w_fp[2], w_fp[14], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];

    // *** DIAGRAM 80 OF 1240 ***
    // Wavefunction(s) for diagram number 80
    // (none)
    // Amplitude(s) for diagram number 80
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[2], w_fp[21], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[90] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[2], w_fp[22], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[73] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[90] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[2], w_fp[23], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[73] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 81 OF 1240 ***
    // Wavefunction(s) for diagram number 81
    // (none)
    // Amplitude(s) for diagram number 81
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[52], w_fp[29], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];

    // *** DIAGRAM 82 OF 1240 ***
    // Wavefunction(s) for diagram number 82
    // (none)
    // Amplitude(s) for diagram number 82
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[48], w_fp[2], w_fp[29], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[90] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];

    // *** DIAGRAM 83 OF 1240 ***
    // Wavefunction(s) for diagram number 83
    // (none)
    // Amplitude(s) for diagram number 83
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[2], w_fp[25], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[90] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 84 OF 1240 ***
    // Wavefunction(s) for diagram number 84
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[38], w_fp[4], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[25] );
    // Amplitude(s) for diagram number 84
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[25], w_fp[52], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 85 OF 1240 ***
    // Wavefunction(s) for diagram number 85
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[38], w_fp[6], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[48] );
    // Amplitude(s) for diagram number 85
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[48], w_fp[52], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 86 OF 1240 ***
    // Wavefunction(s) for diagram number 86
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[38], w_fp[2], COUPs[1], 1.0, 0., 0., w_fp[23] );
    // Amplitude(s) for diagram number 86
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[23], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 87 OF 1240 ***
    // Wavefunction(s) for diagram number 87
    // (none)
    // Amplitude(s) for diagram number 87
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[48], w_fp[2], w_fp[9], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];

    // *** DIAGRAM 88 OF 1240 ***
    // Wavefunction(s) for diagram number 88
    // (none)
    // Amplitude(s) for diagram number 88
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[14], w_fp[23], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 89 OF 1240 ***
    // Wavefunction(s) for diagram number 89
    // (none)
    // Amplitude(s) for diagram number 89
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[25], w_fp[2], w_fp[14], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];

    // *** DIAGRAM 90 OF 1240 ***
    // Wavefunction(s) for diagram number 90
    // (none)
    // Amplitude(s) for diagram number 90
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[2], w_fp[18], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[2], w_fp[19], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[2], w_fp[20], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 91 OF 1240 ***
    // Wavefunction(s) for diagram number 91
    // (none)
    // Amplitude(s) for diagram number 91
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[52], w_fp[27], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];

    // *** DIAGRAM 92 OF 1240 ***
    // Wavefunction(s) for diagram number 92
    // (none)
    // Amplitude(s) for diagram number 92
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[40], w_fp[2], w_fp[27], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];

    // *** DIAGRAM 93 OF 1240 ***
    // Wavefunction(s) for diagram number 93
    // (none)
    // Amplitude(s) for diagram number 93
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[2], w_fp[28], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 94 OF 1240 ***
    // Wavefunction(s) for diagram number 94
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[41], w_fp[4], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[28] );
    // Amplitude(s) for diagram number 94
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[28], w_fp[52], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 95 OF 1240 ***
    // Wavefunction(s) for diagram number 95
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[41], w_fp[5], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[40] );
    // Amplitude(s) for diagram number 95
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[40], w_fp[52], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[24] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 96 OF 1240 ***
    // Wavefunction(s) for diagram number 96
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[41], w_fp[2], COUPs[1], 1.0, 0., 0., w_fp[20] );
    // Amplitude(s) for diagram number 96
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[20], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[24] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[48] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[72] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 97 OF 1240 ***
    // Wavefunction(s) for diagram number 97
    // (none)
    // Amplitude(s) for diagram number 97
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[40], w_fp[2], w_fp[9], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];

    // *** DIAGRAM 98 OF 1240 ***
    // Wavefunction(s) for diagram number 98
    // (none)
    // Amplitude(s) for diagram number 98
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[20], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[48] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[72] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 99 OF 1240 ***
    // Wavefunction(s) for diagram number 99
    // (none)
    // Amplitude(s) for diagram number 99
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[28], w_fp[2], w_fp[12], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[78] += amp_sv[0];

    // *** DIAGRAM 100 OF 1240 ***
    // Wavefunction(s) for diagram number 100
    // (none)
    // Amplitude(s) for diagram number 100
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[2], w_fp[15], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[24] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[2], w_fp[16], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[48] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[72] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[2], w_fp[17], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[24] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[48] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[72] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] -= cxtype( 0, 1 ) * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 26 );
    storeWf( wfs, w_cx, nevt, 27 );
    storeWf( wfs, w_cx, nevt, 28 );
    storeWf( wfs, w_cx, nevt, 29 );
    storeWf( wfs, w_cx, nevt, 30 );
    storeWf( wfs, w_cx, nevt, 31 );
    storeWf( wfs, w_cx, nevt, 32 );
    storeWf( wfs, w_cx, nevt, 33 );
    storeWf( wfs, w_cx, nevt, 34 );
    storeWf( wfs, w_cx, nevt, 35 );
    storeWf( wfs, w_cx, nevt, 36 );
    storeWf( wfs, w_cx, nevt, 37 );
    storeWf( wfs, w_cx, nevt, 38 );
    storeWf( wfs, w_cx, nevt, 39 );
    storeWf( wfs, w_cx, nevt, 40 );
    storeWf( wfs, w_cx, nevt, 41 );
    storeWf( wfs, w_cx, nevt, 42 );
    storeWf( wfs, w_cx, nevt, 43 );
    storeWf( wfs, w_cx, nevt, 44 );
    storeWf( wfs, w_cx, nevt, 45 );
    storeWf( wfs, w_cx, nevt, 46 );
    storeWf( wfs, w_cx, nevt, 47 );
    storeWf( wfs, w_cx, nevt, 48 );
    storeWf( wfs, w_cx, nevt, 49 );
    storeWf( wfs, w_cx, nevt, 50 );
    storeWf( wfs, w_cx, nevt, 51 );
    storeWf( wfs, w_cx, nevt, 52 );
    storeWf( wfs, w_cx, nevt, 53 );
    storeWf( wfs, w_cx, nevt, 54 );
#endif
  }

  //--------------------------------------------------------------------------
}
