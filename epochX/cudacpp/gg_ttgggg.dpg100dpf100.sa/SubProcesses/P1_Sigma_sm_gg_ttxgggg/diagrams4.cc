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
  diagramgroup4( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 8 );
    retrieveWf( wfs, w_cx, nevt, 10 );
    retrieveWf( wfs, w_cx, nevt, 12 );
    retrieveWf( wfs, w_cx, nevt, 13 );
    retrieveWf( wfs, w_cx, nevt, 26 );
    retrieveWf( wfs, w_cx, nevt, 37 );
    retrieveWf( wfs, w_cx, nevt, 39 );
    retrieveWf( wfs, w_cx, nevt, 40 );
    retrieveWf( wfs, w_cx, nevt, 44 );
    retrieveWf( wfs, w_cx, nevt, 46 );
    retrieveWf( wfs, w_cx, nevt, 47 );
    retrieveWf( wfs, w_cx, nevt, 51 );
    retrieveWf( wfs, w_cx, nevt, 52 );
    retrieveWf( wfs, w_cx, nevt, 53 );
    retrieveWf( wfs, w_cx, nevt, 54 );
    retrieveWf( wfs, w_cx, nevt, 55 );
    retrieveWf( wfs, w_cx, nevt, 56 );
    retrieveWf( wfs, w_cx, nevt, 63 );
    retrieveWf( wfs, w_cx, nevt, 64 );
    retrieveWf( wfs, w_cx, nevt, 65 );
    retrieveWf( wfs, w_cx, nevt, 67 );
    retrieveWf( wfs, w_cx, nevt, 84 );
    retrieveWf( wfs, w_cx, nevt, 85 );
    retrieveWf( wfs, w_cx, nevt, 100 );
    retrieveWf( wfs, w_cx, nevt, 101 );
    retrieveWf( wfs, w_cx, nevt, 113 );
    retrieveWf( wfs, w_cx, nevt, 114 );
    retrieveWf( wfs, w_cx, nevt, 116 );
    retrieveWf( wfs, w_cx, nevt, 118 );
    retrieveWf( wfs, w_cx, nevt, 119 );
    retrieveWf( wfs, w_cx, nevt, 120 );
    retrieveWf( wfs, w_cx, nevt, 121 );
    retrieveWf( wfs, w_cx, nevt, 122 );
    retrieveWf( wfs, w_cx, nevt, 123 );
    retrieveWf( wfs, w_cx, nevt, 125 );
    retrieveWf( wfs, w_cx, nevt, 127 );
    retrieveWf( wfs, w_cx, nevt, 128 );
    retrieveWf( wfs, w_cx, nevt, 129 );
    retrieveWf( wfs, w_cx, nevt, 131 );
    retrieveWf( wfs, w_cx, nevt, 151 );
    retrieveWf( wfs, w_cx, nevt, 152 );
    retrieveWf( wfs, w_cx, nevt, 153 );
    retrieveWf( wfs, w_cx, nevt, 154 );
    retrieveWf( wfs, w_cx, nevt, 155 );
    retrieveWf( wfs, w_cx, nevt, 156 );
    retrieveWf( wfs, w_cx, nevt, 157 );
    retrieveWf( wfs, w_cx, nevt, 158 );
    retrieveWf( wfs, w_cx, nevt, 159 );
    retrieveWf( wfs, w_cx, nevt, 160 );
    retrieveWf( wfs, w_cx, nevt, 161 );
    retrieveWf( wfs, w_cx, nevt, 163 );
    retrieveWf( wfs, w_cx, nevt, 164 );
    retrieveWf( wfs, w_cx, nevt, 166 );
    retrieveWf( wfs, w_cx, nevt, 167 );
    retrieveWf( wfs, w_cx, nevt, 169 );
    retrieveWf( wfs, w_cx, nevt, 174 );
    retrieveWf( wfs, w_cx, nevt, 176 );
    retrieveWf( wfs, w_cx, nevt, 177 );
    retrieveWf( wfs, w_cx, nevt, 179 );
#endif
#endif

    // *** DIAGRAM 301 OF 15495 ***
    // Wavefunction(s) for diagram number 301
    // (none)
    // Amplitude(s) for diagram number 301
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[179], w_fp[156], w_fp[114], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[266] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[300] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 302 OF 15495 ***
    // Wavefunction(s) for diagram number 302
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[113], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[184] );
    // Amplitude(s) for diagram number 302
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[184], w_fp[169], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[244] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];

    // *** DIAGRAM 303 OF 15495 ***
    // Wavefunction(s) for diagram number 303
    // (none)
    // Amplitude(s) for diagram number 303
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[169], w_fp[116], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[266] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[268] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[269] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 304 OF 15495 ***
    // Wavefunction(s) for diagram number 304
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[156], w_fp[113], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[185] );
    // Amplitude(s) for diagram number 304
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[157], w_fp[185], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[304] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];

    // *** DIAGRAM 305 OF 15495 ***
    // Wavefunction(s) for diagram number 305
    // (none)
    // Amplitude(s) for diagram number 305
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[157], w_fp[163], w_fp[113], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[352] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];

    // *** DIAGRAM 306 OF 15495 ***
    // Wavefunction(s) for diagram number 306
    // (none)
    // Amplitude(s) for diagram number 306
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[157], w_fp[156], w_fp[116], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[304] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[305] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[328] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 307 OF 15495 ***
    // Wavefunction(s) for diagram number 307
    // (none)
    // Amplitude(s) for diagram number 307
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[114], w_fp[164], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[264] -= amp_sv[0];
    jamp_sv[266] += amp_sv[0];
    jamp_sv[300] -= amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[324] += amp_sv[0];
    jamp_sv[326] -= amp_sv[0];
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[337] += amp_sv[0];
    jamp_sv[342] += amp_sv[0];
    jamp_sv[343] -= amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];

    // *** DIAGRAM 308 OF 15495 ***
    // Wavefunction(s) for diagram number 308
    // (none)
    // Amplitude(s) for diagram number 308
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[163], w_fp[114], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[336] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[343] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 309 OF 15495 ***
    // Wavefunction(s) for diagram number 309
    // (none)
    // Amplitude(s) for diagram number 309
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[44], w_fp[164], w_fp[113], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[244] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[300] -= amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[324] += amp_sv[0];
    jamp_sv[326] -= amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[337] += amp_sv[0];
    jamp_sv[342] += amp_sv[0];
    jamp_sv[343] -= amp_sv[0];

    // *** DIAGRAM 310 OF 15495 ***
    // Wavefunction(s) for diagram number 310
    // (none)
    // Amplitude(s) for diagram number 310
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[185], w_fp[44], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[300] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[304] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[305] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[328] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 311 OF 15495 ***
    // Wavefunction(s) for diagram number 311
    // (none)
    // Amplitude(s) for diagram number 311
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[184], w_fp[156], w_fp[44], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[244] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[268] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[269] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[336] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[343] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 312 OF 15495 ***
    // Wavefunction(s) for diagram number 312
    // (none)
    // Amplitude(s) for diagram number 312
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[164], w_fp[116], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[264] -= amp_sv[0];
    jamp_sv[266] += amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];

    // *** DIAGRAM 313 OF 15495 ***
    // Wavefunction(s) for diagram number 313
    // (none)
    // Amplitude(s) for diagram number 313
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[184], w_fp[163], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[336] += amp_sv[0];
    jamp_sv[337] -= amp_sv[0];
    jamp_sv[342] -= amp_sv[0];
    jamp_sv[343] += amp_sv[0];

    // *** DIAGRAM 314 OF 15495 ***
    // Wavefunction(s) for diagram number 314
    // (none)
    // Amplitude(s) for diagram number 314
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[156], w_fp[121], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[242] += amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[264] += amp_sv[0];
    jamp_sv[266] -= amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[353] += amp_sv[0];
    jamp_sv[358] += amp_sv[0];
    jamp_sv[359] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[156], w_fp[122], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[244] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[300] -= amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[324] += amp_sv[0];
    jamp_sv[326] -= amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[337] += amp_sv[0];
    jamp_sv[342] += amp_sv[0];
    jamp_sv[343] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[156], w_fp[123], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[264] -= amp_sv[0];
    jamp_sv[266] += amp_sv[0];
    jamp_sv[300] -= amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[324] += amp_sv[0];
    jamp_sv[326] -= amp_sv[0];
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[337] += amp_sv[0];
    jamp_sv[342] += amp_sv[0];
    jamp_sv[343] -= amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];

    // *** DIAGRAM 315 OF 15495 ***
    // Wavefunction(s) for diagram number 315
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[100], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[186] );
    // Amplitude(s) for diagram number 315
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[186], w_fp[169], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[242] += amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[266] -= amp_sv[0];
    jamp_sv[267] += amp_sv[0];

    // *** DIAGRAM 316 OF 15495 ***
    // Wavefunction(s) for diagram number 316
    // (none)
    // Amplitude(s) for diagram number 316
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[169], w_fp[125], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[241] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[266] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[267] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[268] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 317 OF 15495 ***
    // Wavefunction(s) for diagram number 317
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[156], w_fp[100], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[187] );
    // Amplitude(s) for diagram number 317
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[157], w_fp[187], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[310] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[353] += amp_sv[0];

    // *** DIAGRAM 318 OF 15495 ***
    // Wavefunction(s) for diagram number 318
    // (none)
    // Amplitude(s) for diagram number 318
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[157], w_fp[161], w_fp[100], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[328] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];

    // *** DIAGRAM 319 OF 15495 ***
    // Wavefunction(s) for diagram number 319
    // (none)
    // Amplitude(s) for diagram number 319
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[157], w_fp[156], w_fp[125], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[310] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[311] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[328] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 320 OF 15495 ***
    // Wavefunction(s) for diagram number 320
    // (none)
    // Amplitude(s) for diagram number 320
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[101], w_fp[164], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[241] += amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[265] -= amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[312] -= amp_sv[0];
    jamp_sv[313] += amp_sv[0];
    jamp_sv[318] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];

    // *** DIAGRAM 321 OF 15495 ***
    // Wavefunction(s) for diagram number 321
    // (none)
    // Amplitude(s) for diagram number 321
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[161], w_fp[101], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[312] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[313] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[328] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 322 OF 15495 ***
    // Wavefunction(s) for diagram number 322
    // (none)
    // Amplitude(s) for diagram number 322
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[37], w_fp[164], w_fp[100], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[242] += amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[266] -= amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[312] -= amp_sv[0];
    jamp_sv[313] += amp_sv[0];
    jamp_sv[318] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[353] += amp_sv[0];

    // *** DIAGRAM 323 OF 15495 ***
    // Wavefunction(s) for diagram number 323
    // (none)
    // Amplitude(s) for diagram number 323
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[187], w_fp[37], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[306] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[310] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[311] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[350] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 324 OF 15495 ***
    // Wavefunction(s) for diagram number 324
    // (none)
    // Amplitude(s) for diagram number 324
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[186], w_fp[156], w_fp[37], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[242] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[266] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[267] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[312] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[313] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 325 OF 15495 ***
    // Wavefunction(s) for diagram number 325
    // (none)
    // Amplitude(s) for diagram number 325
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[164], w_fp[125], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[241] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[243] += amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[265] -= amp_sv[0];
    jamp_sv[266] += amp_sv[0];
    jamp_sv[267] -= amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];

    // *** DIAGRAM 326 OF 15495 ***
    // Wavefunction(s) for diagram number 326
    // (none)
    // Amplitude(s) for diagram number 326
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[186], w_fp[161], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[312] += amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[318] -= amp_sv[0];
    jamp_sv[319] += amp_sv[0];

    // *** DIAGRAM 327 OF 15495 ***
    // Wavefunction(s) for diagram number 327
    // (none)
    // Amplitude(s) for diagram number 327
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[156], w_fp[127], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[242] += amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[266] -= amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[353] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[156], w_fp[128], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[242] += amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[266] -= amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[312] -= amp_sv[0];
    jamp_sv[313] += amp_sv[0];
    jamp_sv[318] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[353] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[156], w_fp[129], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[241] += amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[265] -= amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[312] -= amp_sv[0];
    jamp_sv[313] += amp_sv[0];
    jamp_sv[318] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];

    // *** DIAGRAM 328 OF 15495 ***
    // Wavefunction(s) for diagram number 328
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[84], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[188] );
    // Amplitude(s) for diagram number 328
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[188], w_fp[169], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] += amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[264] -= amp_sv[0];
    jamp_sv[265] += amp_sv[0];

    // *** DIAGRAM 329 OF 15495 ***
    // Wavefunction(s) for diagram number 329
    // (none)
    // Amplitude(s) for diagram number 329
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[169], w_fp[131], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[241] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[267] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[269] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 330 OF 15495 ***
    // Wavefunction(s) for diagram number 330
    // (none)
    // Amplitude(s) for diagram number 330
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[157], w_fp[158], w_fp[84], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[304] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];

    // *** DIAGRAM 331 OF 15495 ***
    // Wavefunction(s) for diagram number 331
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[156], w_fp[84], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[189] );
    // Amplitude(s) for diagram number 331
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[157], w_fp[189], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[334] += amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];

    // *** DIAGRAM 332 OF 15495 ***
    // Wavefunction(s) for diagram number 332
    // (none)
    // Amplitude(s) for diagram number 332
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[157], w_fp[156], w_fp[131], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[304] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[305] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[310] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[311] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 333 OF 15495 ***
    // Wavefunction(s) for diagram number 333
    // (none)
    // Amplitude(s) for diagram number 333
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[26], w_fp[164], w_fp[84], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] += amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[264] -= amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[294] += amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[354] += amp_sv[0];
    jamp_sv[356] -= amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];

    // *** DIAGRAM 334 OF 15495 ***
    // Wavefunction(s) for diagram number 334
    // (none)
    // Amplitude(s) for diagram number 334
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[189], w_fp[26], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[330] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[354] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 335 OF 15495 ***
    // Wavefunction(s) for diagram number 335
    // (none)
    // Amplitude(s) for diagram number 335
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[188], w_fp[156], w_fp[26], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[241] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[294] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[295] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 336 OF 15495 ***
    // Wavefunction(s) for diagram number 336
    // (none)
    // Amplitude(s) for diagram number 336
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[85], w_fp[164], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[243] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[267] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[294] += amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[354] += amp_sv[0];
    jamp_sv[356] -= amp_sv[0];

    // *** DIAGRAM 337 OF 15495 ***
    // Wavefunction(s) for diagram number 337
    // (none)
    // Amplitude(s) for diagram number 337
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[158], w_fp[85], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[288] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[294] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[295] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[304] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[305] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[310] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[311] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 338 OF 15495 ***
    // Wavefunction(s) for diagram number 338
    // (none)
    // Amplitude(s) for diagram number 338
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[164], w_fp[131], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] += amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[264] -= amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];

    // *** DIAGRAM 339 OF 15495 ***
    // Wavefunction(s) for diagram number 339
    // (none)
    // Amplitude(s) for diagram number 339
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[188], w_fp[158], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[288] += amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[294] -= amp_sv[0];
    jamp_sv[295] += amp_sv[0];

    // *** DIAGRAM 340 OF 15495 ***
    // Wavefunction(s) for diagram number 340
    // (none)
    // Amplitude(s) for diagram number 340
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[156], w_fp[120], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[241] += amp_sv[0];
    jamp_sv[243] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[264] += amp_sv[0];
    jamp_sv[265] -= amp_sv[0];
    jamp_sv[267] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[358] += amp_sv[0];
    jamp_sv[359] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[156], w_fp[119], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[243] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[267] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[294] += amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[354] += amp_sv[0];
    jamp_sv[356] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[156], w_fp[118], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] += amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[264] -= amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[294] += amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[354] += amp_sv[0];
    jamp_sv[356] -= amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];

    // *** DIAGRAM 341 OF 15495 ***
    // Wavefunction(s) for diagram number 341
    // (none)
    // Amplitude(s) for diagram number 341
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[169], w_fp[151], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[241] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[267] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[269] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[169], w_fp[152], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[241] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[266] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[267] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[268] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[169], w_fp[153], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[266] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[268] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[269] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 342 OF 15495 ***
    // Wavefunction(s) for diagram number 342
    // (none)
    // Amplitude(s) for diagram number 342
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[157], w_fp[156], w_fp[151], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[304] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[305] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[310] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[311] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[157], w_fp[156], w_fp[152], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[310] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[311] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[328] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[157], w_fp[156], w_fp[153], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[304] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[305] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[328] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 343 OF 15495 ***
    // Wavefunction(s) for diagram number 343
    // (none)
    // Amplitude(s) for diagram number 343
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[156], w_fp[67], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[241] += amp_sv[0];
    jamp_sv[243] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[264] += amp_sv[0];
    jamp_sv[265] -= amp_sv[0];
    jamp_sv[267] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[358] += amp_sv[0];
    jamp_sv[359] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[156], w_fp[154], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[241] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[243] += amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[265] -= amp_sv[0];
    jamp_sv[266] += amp_sv[0];
    jamp_sv[267] -= amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[156], w_fp[155], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[264] -= amp_sv[0];
    jamp_sv[266] += amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];

    // *** DIAGRAM 344 OF 15495 ***
    // Wavefunction(s) for diagram number 344
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[5], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[169] );
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[169], w_fp[4], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[190] );
    // Amplitude(s) for diagram number 344
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[159], w_fp[190], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[430] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 345 OF 15495 ***
    // Wavefunction(s) for diagram number 345
    // (none)
    // Amplitude(s) for diagram number 345
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[160], w_fp[190], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[424] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 346 OF 15495 ***
    // Wavefunction(s) for diagram number 346
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[169], w_fp[6], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[191] );
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[157], w_fp[4], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[192] );
    // Amplitude(s) for diagram number 346
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[192], w_fp[191], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[454] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 347 OF 15495 ***
    // Wavefunction(s) for diagram number 347
    // (none)
    // Amplitude(s) for diagram number 347
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[160], w_fp[191], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[448] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 348 OF 15495 ***
    // Wavefunction(s) for diagram number 348
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[169], w_fp[7], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[193] );
    // Amplitude(s) for diagram number 348
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[192], w_fp[193], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[478] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 349 OF 15495 ***
    // Wavefunction(s) for diagram number 349
    // (none)
    // Amplitude(s) for diagram number 349
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[159], w_fp[193], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[472] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 350 OF 15495 ***
    // Wavefunction(s) for diagram number 350
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[169], COUPs[1], 1.0, 0., 0., w_fp[194] );
    // Amplitude(s) for diagram number 350
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[10], w_fp[194], w_fp[6], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[414] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[10], w_fp[194], w_fp[6], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[414] += amp_sv[0];
    jamp_sv[432] -= amp_sv[0];
    jamp_sv[438] += amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[10], w_fp[194], w_fp[6], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] += amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[415] += amp_sv[0];
    jamp_sv[432] -= amp_sv[0];
    jamp_sv[438] += amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];

    // *** DIAGRAM 351 OF 15495 ***
    // Wavefunction(s) for diagram number 351
    // (none)
    // Amplitude(s) for diagram number 351
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[194], w_fp[7], w_fp[12], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[414] += amp_sv[0];
    jamp_sv[432] -= amp_sv[0];
    jamp_sv[438] += amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];

    // *** DIAGRAM 352 OF 15495 ***
    // Wavefunction(s) for diagram number 352
    // (none)
    // Amplitude(s) for diagram number 352
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[194], w_fp[6], w_fp[13], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] += amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[415] += amp_sv[0];
    jamp_sv[432] -= amp_sv[0];
    jamp_sv[438] += amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];

    // *** DIAGRAM 353 OF 15495 ***
    // Wavefunction(s) for diagram number 353
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[10], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[195] );
    // Amplitude(s) for diagram number 353
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[195], w_fp[191], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[450] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];

    // *** DIAGRAM 354 OF 15495 ***
    // Wavefunction(s) for diagram number 354
    // (none)
    // Amplitude(s) for diagram number 354
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[191], w_fp[13], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[432] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[438] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[454] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 355 OF 15495 ***
    // Wavefunction(s) for diagram number 355
    // (none)
    // Amplitude(s) for diagram number 355
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[195], w_fp[193], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[474] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];

    // *** DIAGRAM 356 OF 15495 ***
    // Wavefunction(s) for diagram number 356
    // (none)
    // Amplitude(s) for diagram number 356
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[193], w_fp[12], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[456] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 357 OF 15495 ***
    // Wavefunction(s) for diagram number 357
    // (none)
    // Amplitude(s) for diagram number 357
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[37], w_fp[194], w_fp[4], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[362] += amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];
    jamp_sv[432] -= amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[438] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[473] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[37], w_fp[194], w_fp[4], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[362] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[414] += amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[432] -= amp_sv[0];
    jamp_sv[438] += amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[473] += amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[37], w_fp[194], w_fp[4], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[363] += amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[414] += amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];

    // *** DIAGRAM 358 OF 15495 ***
    // Wavefunction(s) for diagram number 358
    // (none)
    // Amplitude(s) for diagram number 358
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[194], w_fp[7], w_fp[40], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[362] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[414] += amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[432] -= amp_sv[0];
    jamp_sv[438] += amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[473] += amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];

    // *** DIAGRAM 359 OF 15495 ***
    // Wavefunction(s) for diagram number 359
    // (none)
    // Amplitude(s) for diagram number 359
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[194], w_fp[4], w_fp[39], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[363] += amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[414] += amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];

    // *** DIAGRAM 360 OF 15495 ***
    // Wavefunction(s) for diagram number 360
    // (none)
    // Amplitude(s) for diagram number 360
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[166], w_fp[190], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[426] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];

    // *** DIAGRAM 361 OF 15495 ***
    // Wavefunction(s) for diagram number 361
    // (none)
    // Amplitude(s) for diagram number 361
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[190], w_fp[39], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[408] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 362 OF 15495 ***
    // Wavefunction(s) for diagram number 362
    // (none)
    // Amplitude(s) for diagram number 362
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[166], w_fp[193], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[468] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[473] += amp_sv[0];

    // *** DIAGRAM 363 OF 15495 ***
    // Wavefunction(s) for diagram number 363
    // (none)
    // Amplitude(s) for diagram number 363
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[193], w_fp[40], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[457] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 364 OF 15495 ***
    // Wavefunction(s) for diagram number 364
    // (none)
    // Amplitude(s) for diagram number 364
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[44], w_fp[194], w_fp[4], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[457] += amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[44], w_fp[194], w_fp[4], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[415] += amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[44], w_fp[194], w_fp[4], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[365] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[415] += amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];

    // *** DIAGRAM 365 OF 15495 ***
    // Wavefunction(s) for diagram number 365
    // (none)
    // Amplitude(s) for diagram number 365
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[194], w_fp[6], w_fp[47], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[415] += amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[462] += amp_sv[0];

    // *** DIAGRAM 366 OF 15495 ***
    // Wavefunction(s) for diagram number 366
    // (none)
    // Amplitude(s) for diagram number 366
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[194], w_fp[4], w_fp[46], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[365] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[415] += amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];

    // *** DIAGRAM 367 OF 15495 ***
    // Wavefunction(s) for diagram number 367
    // (none)
    // Amplitude(s) for diagram number 367
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[167], w_fp[190], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[420] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];

    // *** DIAGRAM 368 OF 15495 ***
    // Wavefunction(s) for diagram number 368
    // (none)
    // Amplitude(s) for diagram number 368
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[190], w_fp[46], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[424] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 369 OF 15495 ***
    // Wavefunction(s) for diagram number 369
    // (none)
    // Amplitude(s) for diagram number 369
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[167], w_fp[191], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[444] += amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];

    // *** DIAGRAM 370 OF 15495 ***
    // Wavefunction(s) for diagram number 370
    // (none)
    // Amplitude(s) for diagram number 370
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[191], w_fp[47], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[433] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 371 OF 15495 ***
    // Wavefunction(s) for diagram number 371
    // (none)
    // Amplitude(s) for diagram number 371
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[51], w_fp[194], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[457] += amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[194], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[432] += amp_sv[0];
    jamp_sv[438] -= amp_sv[0];
    jamp_sv[457] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[53], w_fp[194], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] -= amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[432] += amp_sv[0];
    jamp_sv[438] -= amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[456] += amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];

    // *** DIAGRAM 372 OF 15495 ***
    // Wavefunction(s) for diagram number 372
    // (none)
    // Amplitude(s) for diagram number 372
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[193], w_fp[51], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[456] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[193], w_fp[52], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[457] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[193], w_fp[53], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[456] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 373 OF 15495 ***
    // Wavefunction(s) for diagram number 373
    // (none)
    // Amplitude(s) for diagram number 373
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[54], w_fp[194], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] += amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[432] -= amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[438] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[55], w_fp[194], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[456] += amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[56], w_fp[194], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[432] += amp_sv[0];
    jamp_sv[438] -= amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[456] += amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];

    // *** DIAGRAM 374 OF 15495 ***
    // Wavefunction(s) for diagram number 374
    // (none)
    // Amplitude(s) for diagram number 374
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[191], w_fp[54], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[432] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[433] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[438] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[454] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[191], w_fp[55], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[433] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[191], w_fp[56], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[432] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[438] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[454] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 375 OF 15495 ***
    // Wavefunction(s) for diagram number 375
    // (none)
    // Amplitude(s) for diagram number 375
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[63], w_fp[194], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[363] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[414] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[64], w_fp[194], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[457] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[65], w_fp[194], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[457] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];

    // *** DIAGRAM 376 OF 15495 ***
    // Wavefunction(s) for diagram number 376
    // (none)
    // Amplitude(s) for diagram number 376
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[190], w_fp[63], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[408] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[424] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[190], w_fp[64], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[409] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[424] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[190], w_fp[65], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 377 OF 15495 ***
    // Wavefunction(s) for diagram number 377
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[4], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[196] );
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[169], w_fp[8], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[197] );
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[196], w_fp[6], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[198] );
    // Amplitude(s) for diagram number 377
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[198], w_fp[197], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[365] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 378 OF 15495 ***
    // Wavefunction(s) for diagram number 378
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[196], w_fp[7], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[199] );
    // Amplitude(s) for diagram number 378
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[199], w_fp[197], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 379 OF 15495 ***
    // Wavefunction(s) for diagram number 379
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[196], w_fp[8], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[200] );
    // Amplitude(s) for diagram number 379
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[200], w_fp[191], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[450] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 380 OF 15495 ***
    // Wavefunction(s) for diagram number 380
    // (none)
    // Amplitude(s) for diagram number 380
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[200], w_fp[193], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[474] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 381 OF 15495 ***
    // Wavefunction(s) for diagram number 381
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[196], w_fp[169], COUPs[1], 1.0, 0., 0., w_fp[201] );
    // Amplitude(s) for diagram number 381
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[37], w_fp[201], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[433] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 382 OF 15495 ***
    // Wavefunction(s) for diagram number 382
    // (none)
    // Amplitude(s) for diagram number 382
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[193], w_fp[37], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[457] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];

    // *** DIAGRAM 383 OF 15495 ***
    // Wavefunction(s) for diagram number 383
    // (none)
    // Amplitude(s) for diagram number 383
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[199], w_fp[169], w_fp[37], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[363] += amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];

    // *** DIAGRAM 384 OF 15495 ***
    // Wavefunction(s) for diagram number 384
    // (none)
    // Amplitude(s) for diagram number 384
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[44], w_fp[201], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[365] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[433] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 385 OF 15495 ***
    // Wavefunction(s) for diagram number 385
    // (none)
    // Amplitude(s) for diagram number 385
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[191], w_fp[44], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[433] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];

    // *** DIAGRAM 386 OF 15495 ***
    // Wavefunction(s) for diagram number 386
    // (none)
    // Amplitude(s) for diagram number 386
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[198], w_fp[169], w_fp[44], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[365] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];

    // *** DIAGRAM 387 OF 15495 ***
    // Wavefunction(s) for diagram number 387
    // (none)
    // Amplitude(s) for diagram number 387
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[199], w_fp[191], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[433] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 388 OF 15495 ***
    // Wavefunction(s) for diagram number 388
    // (none)
    // Amplitude(s) for diagram number 388
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[198], w_fp[193], w_fp[8], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[457] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 389 OF 15495 ***
    // Wavefunction(s) for diagram number 389
    // (none)
    // Amplitude(s) for diagram number 389
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[169], w_fp[63], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[169], w_fp[64], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[365] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[433] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[169], w_fp[65], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[433] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 390 OF 15495 ***
    // Wavefunction(s) for diagram number 390
    // (none)
    // Amplitude(s) for diagram number 390
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[197], w_fp[84], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[363] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];

    // *** DIAGRAM 391 OF 15495 ***
    // Wavefunction(s) for diagram number 391
    // (none)
    // Amplitude(s) for diagram number 391
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[200], w_fp[169], w_fp[84], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[450] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];

    // *** DIAGRAM 392 OF 15495 ***
    // Wavefunction(s) for diagram number 392
    // (none)
    // Amplitude(s) for diagram number 392
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[169], w_fp[85], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 393 OF 15495 ***
    // Wavefunction(s) for diagram number 393
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[174], w_fp[4], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[202] );
    // Amplitude(s) for diagram number 393
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[202], w_fp[197], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 394 OF 15495 ***
    // Wavefunction(s) for diagram number 394
    // (none)
    // Amplitude(s) for diagram number 394
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[176], w_fp[197], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 395 OF 15495 ***
    // Wavefunction(s) for diagram number 395
    // (none)
    // Amplitude(s) for diagram number 395
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[177], w_fp[190], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[426] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 396 OF 15495 ***
    // Wavefunction(s) for diagram number 396
    // (none)
    // Amplitude(s) for diagram number 396
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[177], w_fp[193], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[468] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 397 OF 15495 ***
    // Wavefunction(s) for diagram number 397
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[174], w_fp[169], COUPs[1], 1.0, 0., 0., w_fp[203] );
    // Amplitude(s) for diagram number 397
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[10], w_fp[203], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[456] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 398 OF 15495 ***
    // Wavefunction(s) for diagram number 398
    // (none)
    // Amplitude(s) for diagram number 398
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[174], w_fp[193], w_fp[10], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[456] += amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];

    // *** DIAGRAM 399 OF 15495 ***
    // Wavefunction(s) for diagram number 399
    // (none)
    // Amplitude(s) for diagram number 399
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[176], w_fp[169], w_fp[10], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] += amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[415] += amp_sv[0];

    // *** DIAGRAM 400 OF 15495 ***
    // Wavefunction(s) for diagram number 400
    // (none)
    // Amplitude(s) for diagram number 400
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[44], w_fp[203], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[456] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] += cxtype( 0, 1 ) * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 169 );
    storeWf( wfs, w_cx, nevt, 184 );
    storeWf( wfs, w_cx, nevt, 185 );
    storeWf( wfs, w_cx, nevt, 186 );
    storeWf( wfs, w_cx, nevt, 187 );
    storeWf( wfs, w_cx, nevt, 188 );
    storeWf( wfs, w_cx, nevt, 189 );
    storeWf( wfs, w_cx, nevt, 190 );
    storeWf( wfs, w_cx, nevt, 191 );
    storeWf( wfs, w_cx, nevt, 192 );
    storeWf( wfs, w_cx, nevt, 193 );
    storeWf( wfs, w_cx, nevt, 194 );
    storeWf( wfs, w_cx, nevt, 195 );
    storeWf( wfs, w_cx, nevt, 196 );
    storeWf( wfs, w_cx, nevt, 197 );
    storeWf( wfs, w_cx, nevt, 198 );
    storeWf( wfs, w_cx, nevt, 199 );
    storeWf( wfs, w_cx, nevt, 200 );
    storeWf( wfs, w_cx, nevt, 201 );
    storeWf( wfs, w_cx, nevt, 202 );
    storeWf( wfs, w_cx, nevt, 203 );
#endif
#endif
  }

  //--------------------------------------------------------------------------
}
