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
  diagramgroup54( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 9 );
    retrieveWf( wfs, w_cx, nevt, 16 );
    retrieveWf( wfs, w_cx, nevt, 27 );
    retrieveWf( wfs, w_cx, nevt, 84 );
    retrieveWf( wfs, w_cx, nevt, 102 );
    retrieveWf( wfs, w_cx, nevt, 128 );
    retrieveWf( wfs, w_cx, nevt, 144 );
    retrieveWf( wfs, w_cx, nevt, 150 );
    retrieveWf( wfs, w_cx, nevt, 156 );
    retrieveWf( wfs, w_cx, nevt, 163 );
    retrieveWf( wfs, w_cx, nevt, 164 );
    retrieveWf( wfs, w_cx, nevt, 174 );
    retrieveWf( wfs, w_cx, nevt, 176 );
    retrieveWf( wfs, w_cx, nevt, 179 );
    retrieveWf( wfs, w_cx, nevt, 181 );
    retrieveWf( wfs, w_cx, nevt, 188 );
    retrieveWf( wfs, w_cx, nevt, 225 );
    retrieveWf( wfs, w_cx, nevt, 228 );
    retrieveWf( wfs, w_cx, nevt, 258 );
    retrieveWf( wfs, w_cx, nevt, 263 );
    retrieveWf( wfs, w_cx, nevt, 265 );
    retrieveWf( wfs, w_cx, nevt, 275 );
    retrieveWf( wfs, w_cx, nevt, 285 );
    retrieveWf( wfs, w_cx, nevt, 288 );
    retrieveWf( wfs, w_cx, nevt, 304 );
    retrieveWf( wfs, w_cx, nevt, 312 );
    retrieveWf( wfs, w_cx, nevt, 314 );
    retrieveWf( wfs, w_cx, nevt, 318 );
    retrieveWf( wfs, w_cx, nevt, 319 );
    retrieveWf( wfs, w_cx, nevt, 320 );
    retrieveWf( wfs, w_cx, nevt, 324 );
    retrieveWf( wfs, w_cx, nevt, 437 );
    retrieveWf( wfs, w_cx, nevt, 438 );
    retrieveWf( wfs, w_cx, nevt, 444 );
    retrieveWf( wfs, w_cx, nevt, 447 );
    retrieveWf( wfs, w_cx, nevt, 449 );
    retrieveWf( wfs, w_cx, nevt, 451 );
    retrieveWf( wfs, w_cx, nevt, 452 );
    retrieveWf( wfs, w_cx, nevt, 454 );
    retrieveWf( wfs, w_cx, nevt, 474 );
    retrieveWf( wfs, w_cx, nevt, 479 );
    retrieveWf( wfs, w_cx, nevt, 480 );
    retrieveWf( wfs, w_cx, nevt, 481 );
    retrieveWf( wfs, w_cx, nevt, 483 );
    retrieveWf( wfs, w_cx, nevt, 488 );
    retrieveWf( wfs, w_cx, nevt, 490 );
    retrieveWf( wfs, w_cx, nevt, 496 );
    retrieveWf( wfs, w_cx, nevt, 502 );
    retrieveWf( wfs, w_cx, nevt, 509 );
    retrieveWf( wfs, w_cx, nevt, 510 );
    retrieveWf( wfs, w_cx, nevt, 515 );
    retrieveWf( wfs, w_cx, nevt, 520 );
    retrieveWf( wfs, w_cx, nevt, 521 );
    retrieveWf( wfs, w_cx, nevt, 528 );
    retrieveWf( wfs, w_cx, nevt, 530 );
    retrieveWf( wfs, w_cx, nevt, 532 );
    retrieveWf( wfs, w_cx, nevt, 536 );
    retrieveWf( wfs, w_cx, nevt, 539 );
    retrieveWf( wfs, w_cx, nevt, 545 );
    retrieveWf( wfs, w_cx, nevt, 547 );
    retrieveWf( wfs, w_cx, nevt, 551 );
    retrieveWf( wfs, w_cx, nevt, 576 );
    retrieveWf( wfs, w_cx, nevt, 577 );
    retrieveWf( wfs, w_cx, nevt, 578 );
    retrieveWf( wfs, w_cx, nevt, 579 );
    retrieveWf( wfs, w_cx, nevt, 581 );
    retrieveWf( wfs, w_cx, nevt, 585 );
    retrieveWf( wfs, w_cx, nevt, 586 );
    retrieveWf( wfs, w_cx, nevt, 588 );
    retrieveWf( wfs, w_cx, nevt, 589 );
    retrieveWf( wfs, w_cx, nevt, 590 );
#endif
#endif

    // *** DIAGRAM 5301 OF 15495 ***
    // Wavefunction(s) for diagram number 5301
    // (none)
    // Amplitude(s) for diagram number 5301
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[258], w_fp[150], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[151] -= amp_sv[0];
    jamp_sv[162] -= amp_sv[0];
    jamp_sv[164] += amp_sv[0];
    jamp_sv[265] -= amp_sv[0];
    jamp_sv[271] += amp_sv[0];
    jamp_sv[282] += amp_sv[0];
    jamp_sv[284] -= amp_sv[0];
    jamp_sv[378] += amp_sv[0];
    jamp_sv[380] -= amp_sv[0];
    jamp_sv[612] += amp_sv[0];
    jamp_sv[614] -= amp_sv[0];
    jamp_sv[672] -= amp_sv[0];
    jamp_sv[674] += amp_sv[0];

    // *** DIAGRAM 5302 OF 15495 ***
    // Wavefunction(s) for diagram number 5302
    // (none)
    // Amplitude(s) for diagram number 5302
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[176], w_fp[483], w_fp[479], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[145] += amp_sv[0];
    jamp_sv[151] -= amp_sv[0];
    jamp_sv[265] -= amp_sv[0];
    jamp_sv[271] += amp_sv[0];

    // *** DIAGRAM 5303 OF 15495 ***
    // Wavefunction(s) for diagram number 5303
    // (none)
    // Amplitude(s) for diagram number 5303
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[265], w_fp[150], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[162] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[282] -= amp_sv[0];
    jamp_sv[284] += amp_sv[0];
    jamp_sv[361] += amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[380] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[654] += amp_sv[0];
    jamp_sv[656] -= amp_sv[0];

    // *** DIAGRAM 5304 OF 15495 ***
    // Wavefunction(s) for diagram number 5304
    // (none)
    // Amplitude(s) for diagram number 5304
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[174], w_fp[2], w_fp[490], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[162] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[282] -= amp_sv[0];
    jamp_sv[284] += amp_sv[0];
    jamp_sv[361] += amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[380] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[654] += amp_sv[0];
    jamp_sv[656] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[174], w_fp[2], w_fp[585], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[151] += amp_sv[0];
    jamp_sv[162] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[271] -= amp_sv[0];
    jamp_sv[282] -= amp_sv[0];
    jamp_sv[284] += amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[380] += amp_sv[0];
    jamp_sv[612] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[174], w_fp[2], w_fp[586], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[151] += amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[271] -= amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[612] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[632] -= amp_sv[0];
    jamp_sv[654] -= amp_sv[0];
    jamp_sv[656] += amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    jamp_sv[674] -= amp_sv[0];

    // *** DIAGRAM 5305 OF 15495 ***
    // Wavefunction(s) for diagram number 5305
    // (none)
    // Amplitude(s) for diagram number 5305
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[521], w_fp[144], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[150] += amp_sv[0];
    jamp_sv[264] += amp_sv[0];
    jamp_sv[270] -= amp_sv[0];
    jamp_sv[360] -= amp_sv[0];
    jamp_sv[366] += amp_sv[0];
    jamp_sv[492] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[510] += amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[534] -= amp_sv[0];
    jamp_sv[536] += amp_sv[0];
    jamp_sv[552] += amp_sv[0];
    jamp_sv[554] -= amp_sv[0];

    // *** DIAGRAM 5306 OF 15495 ***
    // Wavefunction(s) for diagram number 5306
    // (none)
    // Amplitude(s) for diagram number 5306
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[181], w_fp[2], w_fp[521], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5307 OF 15495 ***
    // Wavefunction(s) for diagram number 5307
    // (none)
    // Amplitude(s) for diagram number 5307
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[285], w_fp[532], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[60] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[372] -= amp_sv[0];
    jamp_sv[374] += amp_sv[0];

    // *** DIAGRAM 5308 OF 15495 ***
    // Wavefunction(s) for diagram number 5308
    // (none)
    // Amplitude(s) for diagram number 5308
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[179], w_fp[532], w_fp[263], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5309 OF 15495 ***
    // Wavefunction(s) for diagram number 5309
    // (none)
    // Amplitude(s) for diagram number 5309
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[181], w_fp[532], w_fp[258], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[360] -= amp_sv[0];
    jamp_sv[366] += amp_sv[0];

    // *** DIAGRAM 5310 OF 15495 ***
    // Wavefunction(s) for diagram number 5310
    // (none)
    // Amplitude(s) for diagram number 5310
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[536], w_fp[483], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[156] += amp_sv[0];
    jamp_sv[158] -= amp_sv[0];
    jamp_sv[276] -= amp_sv[0];
    jamp_sv[278] += amp_sv[0];

    // *** DIAGRAM 5311 OF 15495 ***
    // Wavefunction(s) for diagram number 5311
    // (none)
    // Amplitude(s) for diagram number 5311
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[536], w_fp[2], w_fp[263], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[156] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5312 OF 15495 ***
    // Wavefunction(s) for diagram number 5312
    // (none)
    // Amplitude(s) for diagram number 5312
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[179], w_fp[483], w_fp[547], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5313 OF 15495 ***
    // Wavefunction(s) for diagram number 5313
    // (none)
    // Amplitude(s) for diagram number 5313
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[285], w_fp[2], w_fp[547], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[60] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5314 OF 15495 ***
    // Wavefunction(s) for diagram number 5314
    // (none)
    // Amplitude(s) for diagram number 5314
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[547], w_fp[258], w_fp[144], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[62] += amp_sv[0];
    jamp_sv[144] += amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[156] -= amp_sv[0];
    jamp_sv[158] += amp_sv[0];
    jamp_sv[264] -= amp_sv[0];
    jamp_sv[270] += amp_sv[0];
    jamp_sv[276] += amp_sv[0];
    jamp_sv[278] -= amp_sv[0];
    jamp_sv[372] += amp_sv[0];
    jamp_sv[374] -= amp_sv[0];
    jamp_sv[492] += amp_sv[0];
    jamp_sv[494] -= amp_sv[0];
    jamp_sv[552] -= amp_sv[0];
    jamp_sv[554] += amp_sv[0];

    // *** DIAGRAM 5315 OF 15495 ***
    // Wavefunction(s) for diagram number 5315
    // (none)
    // Amplitude(s) for diagram number 5315
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[181], w_fp[483], w_fp[479], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[144] += amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[264] -= amp_sv[0];
    jamp_sv[270] += amp_sv[0];

    // *** DIAGRAM 5316 OF 15495 ***
    // Wavefunction(s) for diagram number 5316
    // (none)
    // Amplitude(s) for diagram number 5316
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[263], w_fp[144], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[60] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[156] += amp_sv[0];
    jamp_sv[158] -= amp_sv[0];
    jamp_sv[276] -= amp_sv[0];
    jamp_sv[278] += amp_sv[0];
    jamp_sv[360] += amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[372] -= amp_sv[0];
    jamp_sv[374] += amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[534] += amp_sv[0];
    jamp_sv[536] -= amp_sv[0];

    // *** DIAGRAM 5317 OF 15495 ***
    // Wavefunction(s) for diagram number 5317
    // (none)
    // Amplitude(s) for diagram number 5317
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[179], w_fp[2], w_fp[510], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[60] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[156] += amp_sv[0];
    jamp_sv[158] -= amp_sv[0];
    jamp_sv[276] -= amp_sv[0];
    jamp_sv[278] += amp_sv[0];
    jamp_sv[360] += amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[372] -= amp_sv[0];
    jamp_sv[374] += amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[534] += amp_sv[0];
    jamp_sv[536] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[179], w_fp[2], w_fp[481], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[60] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[150] += amp_sv[0];
    jamp_sv[156] += amp_sv[0];
    jamp_sv[158] -= amp_sv[0];
    jamp_sv[264] += amp_sv[0];
    jamp_sv[270] -= amp_sv[0];
    jamp_sv[276] -= amp_sv[0];
    jamp_sv[278] += amp_sv[0];
    jamp_sv[372] -= amp_sv[0];
    jamp_sv[374] += amp_sv[0];
    jamp_sv[492] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[552] += amp_sv[0];
    jamp_sv[554] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[179], w_fp[2], w_fp[480], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[150] += amp_sv[0];
    jamp_sv[264] += amp_sv[0];
    jamp_sv[270] -= amp_sv[0];
    jamp_sv[360] -= amp_sv[0];
    jamp_sv[366] += amp_sv[0];
    jamp_sv[492] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[510] += amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[534] -= amp_sv[0];
    jamp_sv[536] += amp_sv[0];
    jamp_sv[552] += amp_sv[0];
    jamp_sv[554] -= amp_sv[0];

    // *** DIAGRAM 5318 OF 15495 ***
    // Wavefunction(s) for diagram number 5318
    // (none)
    // Amplitude(s) for diagram number 5318
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[128], w_fp[521], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[580] -= amp_sv[0];
    jamp_sv[581] += amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    jamp_sv[591] += amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[595] -= amp_sv[0];
    jamp_sv[700] += amp_sv[0];
    jamp_sv[701] -= amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[710] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];
    jamp_sv[714] -= amp_sv[0];
    jamp_sv[715] += amp_sv[0];

    // *** DIAGRAM 5319 OF 15495 ***
    // Wavefunction(s) for diagram number 5319
    // (none)
    // Amplitude(s) for diagram number 5319
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[188], w_fp[2], w_fp[521], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[144] += amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[151] += amp_sv[0];
    jamp_sv[264] -= amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[270] += amp_sv[0];
    jamp_sv[271] -= amp_sv[0];
    jamp_sv[360] += amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[367] += amp_sv[0];

    // *** DIAGRAM 5320 OF 15495 ***
    // Wavefunction(s) for diagram number 5320
    // (none)
    // Amplitude(s) for diagram number 5320
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[275], w_fp[532], w_fp[84], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5321 OF 15495 ***
    // Wavefunction(s) for diagram number 5321
    // (none)
    // Amplitude(s) for diagram number 5321
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[532], w_fp[288], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    jamp_sv[360] += amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[376] -= amp_sv[0];
    jamp_sv[377] += amp_sv[0];
    jamp_sv[382] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];

    // *** DIAGRAM 5322 OF 15495 ***
    // Wavefunction(s) for diagram number 5322
    // (none)
    // Amplitude(s) for diagram number 5322
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[188], w_fp[532], w_fp[258], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5323 OF 15495 ***
    // Wavefunction(s) for diagram number 5323
    // (none)
    // Amplitude(s) for diagram number 5323
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[530], w_fp[483], w_fp[84], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5324 OF 15495 ***
    // Wavefunction(s) for diagram number 5324
    // (none)
    // Amplitude(s) for diagram number 5324
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[530], w_fp[2], w_fp[288], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[161] += amp_sv[0];
    jamp_sv[166] += amp_sv[0];
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[286] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    jamp_sv[591] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[710] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];

    // *** DIAGRAM 5325 OF 15495 ***
    // Wavefunction(s) for diagram number 5325
    // (none)
    // Amplitude(s) for diagram number 5325
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[530], w_fp[128], w_fp[258], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5326 OF 15495 ***
    // Wavefunction(s) for diagram number 5326
    // (none)
    // Amplitude(s) for diagram number 5326
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[483], w_fp[581], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[150] += amp_sv[0];
    jamp_sv[151] -= amp_sv[0];
    jamp_sv[160] += amp_sv[0];
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[167] += amp_sv[0];
    jamp_sv[264] += amp_sv[0];
    jamp_sv[265] -= amp_sv[0];
    jamp_sv[270] -= amp_sv[0];
    jamp_sv[271] += amp_sv[0];
    jamp_sv[280] -= amp_sv[0];
    jamp_sv[281] += amp_sv[0];
    jamp_sv[286] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];

    // *** DIAGRAM 5327 OF 15495 ***
    // Wavefunction(s) for diagram number 5327
    // (none)
    // Amplitude(s) for diagram number 5327
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[275], w_fp[2], w_fp[581], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[382] -= amp_sv[0];
    jamp_sv[383] += amp_sv[0];
    jamp_sv[580] += amp_sv[0];
    jamp_sv[581] -= amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[595] += amp_sv[0];
    jamp_sv[700] -= amp_sv[0];
    jamp_sv[701] += amp_sv[0];
    jamp_sv[714] += amp_sv[0];
    jamp_sv[715] -= amp_sv[0];

    // *** DIAGRAM 5328 OF 15495 ***
    // Wavefunction(s) for diagram number 5328
    // (none)
    // Amplitude(s) for diagram number 5328
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[188], w_fp[483], w_fp[479], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5329 OF 15495 ***
    // Wavefunction(s) for diagram number 5329
    // (none)
    // Amplitude(s) for diagram number 5329
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[275], w_fp[128], w_fp[479], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5330 OF 15495 ***
    // Wavefunction(s) for diagram number 5330
    // (none)
    // Amplitude(s) for diagram number 5330
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[324], w_fp[589], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[371] += amp_sv[0];

    // *** DIAGRAM 5331 OF 15495 ***
    // Wavefunction(s) for diagram number 5331
    // (none)
    // Amplitude(s) for diagram number 5331
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[324], w_fp[590], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[68] += amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[380] -= amp_sv[0];
    jamp_sv[381] += amp_sv[0];

    // *** DIAGRAM 5332 OF 15495 ***
    // Wavefunction(s) for diagram number 5332
    // (none)
    // Amplitude(s) for diagram number 5332
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[312], w_fp[7], w_fp[454], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[60] += amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[69] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[362] += amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[372] -= amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[380] += amp_sv[0];
    jamp_sv[381] -= amp_sv[0];
    jamp_sv[382] += amp_sv[0];

    // *** DIAGRAM 5333 OF 15495 ***
    // Wavefunction(s) for diagram number 5333
    // (none)
    // Amplitude(s) for diagram number 5333
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[590], w_fp[312], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[67] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5334 OF 15495 ***
    // Wavefunction(s) for diagram number 5334
    // (none)
    // Amplitude(s) for diagram number 5334
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[314], w_fp[4], w_fp[454], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[363] += amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[373] -= amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[382] += amp_sv[0];

    // *** DIAGRAM 5335 OF 15495 ***
    // Wavefunction(s) for diagram number 5335
    // (none)
    // Amplitude(s) for diagram number 5335
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[589], w_fp[314], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5336 OF 15495 ***
    // Wavefunction(s) for diagram number 5336
    // (none)
    // Amplitude(s) for diagram number 5336
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[532], w_fp[318], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[51] += amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[60] += amp_sv[0];
    jamp_sv[61] -= amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[69] += amp_sv[0];
    jamp_sv[362] += amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[371] += amp_sv[0];
    jamp_sv[372] -= amp_sv[0];
    jamp_sv[373] += amp_sv[0];
    jamp_sv[380] += amp_sv[0];
    jamp_sv[381] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[532], w_fp[319], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[61] -= amp_sv[0];
    jamp_sv[67] -= amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[366] += amp_sv[0];
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[371] += amp_sv[0];
    jamp_sv[373] += amp_sv[0];
    jamp_sv[379] += amp_sv[0];
    jamp_sv[382] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[532], w_fp[320], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[50] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[67] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[366] += amp_sv[0];
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[372] += amp_sv[0];
    jamp_sv[379] += amp_sv[0];
    jamp_sv[380] -= amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[382] -= amp_sv[0];

    // *** DIAGRAM 5337 OF 15495 ***
    // Wavefunction(s) for diagram number 5337
    // (none)
    // Amplitude(s) for diagram number 5337
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[576], w_fp[496], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[212] += amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[525] += amp_sv[0];

    // *** DIAGRAM 5338 OF 15495 ***
    // Wavefunction(s) for diagram number 5338
    // (none)
    // Amplitude(s) for diagram number 5338
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[449], w_fp[496], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[202] += amp_sv[0];
    jamp_sv[203] -= amp_sv[0];
    jamp_sv[514] -= amp_sv[0];
    jamp_sv[515] += amp_sv[0];

    // *** DIAGRAM 5339 OF 15495 ***
    // Wavefunction(s) for diagram number 5339
    // (none)
    // Amplitude(s) for diagram number 5339
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[312], w_fp[7], w_fp[520], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[203] += amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[645] -= amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[668] += amp_sv[0];
    jamp_sv[669] -= amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];

    // *** DIAGRAM 5340 OF 15495 ***
    // Wavefunction(s) for diagram number 5340
    // (none)
    // Amplitude(s) for diagram number 5340
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[449], w_fp[2], w_fp[312], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5341 OF 15495 ***
    // Wavefunction(s) for diagram number 5341
    // (none)
    // Amplitude(s) for diagram number 5341
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[314], w_fp[4], w_fp[520], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[213] += amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[346] -= amp_sv[0];
    jamp_sv[347] += amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[645] -= amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];

    // *** DIAGRAM 5342 OF 15495 ***
    // Wavefunction(s) for diagram number 5342
    // (none)
    // Amplitude(s) for diagram number 5342
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[576], w_fp[2], w_fp[314], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5343 OF 15495 ***
    // Wavefunction(s) for diagram number 5343
    // (none)
    // Amplitude(s) for diagram number 5343
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[530], w_fp[2], w_fp[318], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[203] += amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[346] += amp_sv[0];
    jamp_sv[347] -= amp_sv[0];
    jamp_sv[356] -= amp_sv[0];
    jamp_sv[357] += amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[668] += amp_sv[0];
    jamp_sv[669] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[530], w_fp[2], w_fp[319], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[212] += amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[280] -= amp_sv[0];
    jamp_sv[281] += amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[323] -= amp_sv[0];
    jamp_sv[346] += amp_sv[0];
    jamp_sv[347] -= amp_sv[0];
    jamp_sv[356] -= amp_sv[0];
    jamp_sv[357] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[704] += amp_sv[0];
    jamp_sv[705] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[530], w_fp[2], w_fp[320], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[202] += amp_sv[0];
    jamp_sv[203] -= amp_sv[0];
    jamp_sv[280] -= amp_sv[0];
    jamp_sv[281] += amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[323] -= amp_sv[0];
    jamp_sv[514] -= amp_sv[0];
    jamp_sv[515] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[658] += amp_sv[0];
    jamp_sv[659] -= amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    jamp_sv[669] += amp_sv[0];
    jamp_sv[704] += amp_sv[0];
    jamp_sv[705] -= amp_sv[0];

    // *** DIAGRAM 5344 OF 15495 ***
    // Wavefunction(s) for diagram number 5344
    // (none)
    // Amplitude(s) for diagram number 5344
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[496], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5345 OF 15495 ***
    // Wavefunction(s) for diagram number 5345
    // (none)
    // Amplitude(s) for diagram number 5345
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[496], w_fp[539], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[200] -= amp_sv[0];
    jamp_sv[204] += amp_sv[0];
    jamp_sv[211] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[213] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[516] -= amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[526] += amp_sv[0];

    // *** DIAGRAM 5346 OF 15495 ***
    // Wavefunction(s) for diagram number 5346
    // (none)
    // Amplitude(s) for diagram number 5346
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[324], w_fp[588], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5347 OF 15495 ***
    // Wavefunction(s) for diagram number 5347
    // (none)
    // Amplitude(s) for diagram number 5347
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[324], w_fp[2], w_fp[539], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[293] += amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[615] -= amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[675] += amp_sv[0];

    // *** DIAGRAM 5348 OF 15495 ***
    // Wavefunction(s) for diagram number 5348
    // (none)
    // Amplitude(s) for diagram number 5348
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[588], w_fp[314], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[246] += amp_sv[0];
    jamp_sv[248] -= amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[290] += amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[366] += amp_sv[0];
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[371] += amp_sv[0];

    // *** DIAGRAM 5349 OF 15495 ***
    // Wavefunction(s) for diagram number 5349
    // (none)
    // Amplitude(s) for diagram number 5349
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[2], w_fp[314], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[526] -= amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[646] -= amp_sv[0];
    jamp_sv[703] -= amp_sv[0];
    jamp_sv[704] += amp_sv[0];
    jamp_sv[705] -= amp_sv[0];
    jamp_sv[706] += amp_sv[0];

    // *** DIAGRAM 5350 OF 15495 ***
    // Wavefunction(s) for diagram number 5350
    // (none)
    // Amplitude(s) for diagram number 5350
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[437], w_fp[496], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5351 OF 15495 ***
    // Wavefunction(s) for diagram number 5351
    // (none)
    // Amplitude(s) for diagram number 5351
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[496], w_fp[578], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[200] -= amp_sv[0];
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[203] += amp_sv[0];
    jamp_sv[205] += amp_sv[0];
    jamp_sv[211] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[517] -= amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[526] += amp_sv[0];

    // *** DIAGRAM 5352 OF 15495 ***
    // Wavefunction(s) for diagram number 5352
    // (none)
    // Amplitude(s) for diagram number 5352
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[324], w_fp[438], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5353 OF 15495 ***
    // Wavefunction(s) for diagram number 5353
    // (none)
    // Amplitude(s) for diagram number 5353
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[324], w_fp[2], w_fp[578], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[69] += amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[293] += amp_sv[0];
    jamp_sv[338] -= amp_sv[0];
    jamp_sv[339] += amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[349] -= amp_sv[0];
    jamp_sv[380] += amp_sv[0];
    jamp_sv[381] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[615] -= amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[675] += amp_sv[0];

    // *** DIAGRAM 5354 OF 15495 ***
    // Wavefunction(s) for diagram number 5354
    // (none)
    // Amplitude(s) for diagram number 5354
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[438], w_fp[312], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[67] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[379] += amp_sv[0];
    jamp_sv[380] -= amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[382] -= amp_sv[0];
    jamp_sv[613] += amp_sv[0];
    jamp_sv[614] -= amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[616] -= amp_sv[0];
    jamp_sv[673] -= amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    jamp_sv[676] += amp_sv[0];

    // *** DIAGRAM 5355 OF 15495 ***
    // Wavefunction(s) for diagram number 5355
    // (none)
    // Amplitude(s) for diagram number 5355
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[437], w_fp[2], w_fp[312], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[198] -= amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[203] -= amp_sv[0];
    jamp_sv[276] += amp_sv[0];
    jamp_sv[278] -= amp_sv[0];
    jamp_sv[280] -= amp_sv[0];
    jamp_sv[281] += amp_sv[0];
    jamp_sv[318] -= amp_sv[0];
    jamp_sv[320] += amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[323] -= amp_sv[0];
    jamp_sv[510] += amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[514] -= amp_sv[0];
    jamp_sv[515] += amp_sv[0];

    // *** DIAGRAM 5356 OF 15495 ***
    // Wavefunction(s) for diagram number 5356
    // (none)
    // Amplitude(s) for diagram number 5356
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[496], w_fp[545], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[195] += amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[203] -= amp_sv[0];
    jamp_sv[204] += amp_sv[0];
    jamp_sv[205] -= amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[213] += amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    jamp_sv[514] -= amp_sv[0];
    jamp_sv[515] += amp_sv[0];
    jamp_sv[516] -= amp_sv[0];
    jamp_sv[517] += amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[496], w_fp[528], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[195] += amp_sv[0];
    jamp_sv[198] -= amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[203] -= amp_sv[0];
    jamp_sv[205] -= amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    jamp_sv[510] += amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[514] -= amp_sv[0];
    jamp_sv[515] += amp_sv[0];
    jamp_sv[517] += amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[526] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[496], w_fp[447], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[194] += amp_sv[0];
    jamp_sv[198] -= amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[204] -= amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[510] += amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[516] += amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[526] -= amp_sv[0];

    // *** DIAGRAM 5357 OF 15495 ***
    // Wavefunction(s) for diagram number 5357
    // (none)
    // Amplitude(s) for diagram number 5357
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[324], w_fp[2], w_fp[545], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[338] += amp_sv[0];
    jamp_sv[339] -= amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[349] += amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[380] -= amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[324], w_fp[2], w_fp[528], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[68] += amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[338] += amp_sv[0];
    jamp_sv[339] -= amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[349] += amp_sv[0];
    jamp_sv[380] -= amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[614] -= amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[324], w_fp[2], w_fp[447], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[371] += amp_sv[0];
    jamp_sv[614] -= amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[650] += amp_sv[0];
    jamp_sv[651] -= amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[661] += amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[675] -= amp_sv[0];

    // *** DIAGRAM 5358 OF 15495 ***
    // Wavefunction(s) for diagram number 5358
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[479], w_fp[304], COUPs[0], 1.0, 0., 0., w_fp[521] );
    // Amplitude(s) for diagram number 5358
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[521], w_fp[27], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5359 OF 15495 ***
    // Wavefunction(s) for diagram number 5359
    // (none)
    // Amplitude(s) for diagram number 5359
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[521], w_fp[16], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5360 OF 15495 ***
    // Wavefunction(s) for diagram number 5360
    // (none)
    // Amplitude(s) for diagram number 5360
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[4], w_fp[7], w_fp[521], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[4], w_fp[7], w_fp[521], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[4], w_fp[7], w_fp[521], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5361 OF 15495 ***
    // Wavefunction(s) for diagram number 5361
    // (none)
    // Amplitude(s) for diagram number 5361
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[474], w_fp[312], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5362 OF 15495 ***
    // Wavefunction(s) for diagram number 5362
    // (none)
    // Amplitude(s) for diagram number 5362
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[474], w_fp[314], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5363 OF 15495 ***
    // Wavefunction(s) for diagram number 5363
    // (none)
    // Amplitude(s) for diagram number 5363
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[304], w_fp[4], w_fp[7], w_fp[474], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[304], w_fp[4], w_fp[7], w_fp[474], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[304], w_fp[4], w_fp[7], w_fp[474], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5364 OF 15495 ***
    // Wavefunction(s) for diagram number 5364
    // (none)
    // Amplitude(s) for diagram number 5364
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[577], w_fp[551], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5365 OF 15495 ***
    // Wavefunction(s) for diagram number 5365
    // (none)
    // Amplitude(s) for diagram number 5365
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[577], w_fp[314], w_fp[9], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5366 OF 15495 ***
    // Wavefunction(s) for diagram number 5366
    // (none)
    // Amplitude(s) for diagram number 5366
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[577], w_fp[304], w_fp[16], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5367 OF 15495 ***
    // Wavefunction(s) for diagram number 5367
    // (none)
    // Amplitude(s) for diagram number 5367
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[304], w_fp[9], w_fp[7], w_fp[577], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[304], w_fp[9], w_fp[7], w_fp[577], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[304], w_fp[9], w_fp[7], w_fp[577], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5368 OF 15495 ***
    // Wavefunction(s) for diagram number 5368
    // (none)
    // Amplitude(s) for diagram number 5368
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[551], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5369 OF 15495 ***
    // Wavefunction(s) for diagram number 5369
    // (none)
    // Amplitude(s) for diagram number 5369
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[312], w_fp[9], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5370 OF 15495 ***
    // Wavefunction(s) for diagram number 5370
    // (none)
    // Amplitude(s) for diagram number 5370
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[304], w_fp[27], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[67] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5371 OF 15495 ***
    // Wavefunction(s) for diagram number 5371
    // (none)
    // Amplitude(s) for diagram number 5371
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[304], w_fp[9], w_fp[4], w_fp[488], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[304], w_fp[9], w_fp[4], w_fp[488], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[304], w_fp[9], w_fp[4], w_fp[488], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5372 OF 15495 ***
    // Wavefunction(s) for diagram number 5372
    // (none)
    // Amplitude(s) for diagram number 5372
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[312], w_fp[16], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5373 OF 15495 ***
    // Wavefunction(s) for diagram number 5373
    // (none)
    // Amplitude(s) for diagram number 5373
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[314], w_fp[27], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5374 OF 15495 ***
    // Wavefunction(s) for diagram number 5374
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[479], w_fp[304], w_fp[4], COUPs[2], 1.0, 0., 0., w_fp[480] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[479], w_fp[304], w_fp[4], COUPs[2], 1.0, 0., 0., w_fp[481] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[479], w_fp[304], w_fp[4], COUPs[2], 1.0, 0., 0., w_fp[510] );
    // Amplitude(s) for diagram number 5374
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[7], w_fp[480], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[7], w_fp[481], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[7], w_fp[510], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5375 OF 15495 ***
    // Wavefunction(s) for diagram number 5375
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[479], w_fp[304], w_fp[7], COUPs[2], 1.0, 0., 0., w_fp[586] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[479], w_fp[304], w_fp[7], COUPs[2], 1.0, 0., 0., w_fp[585] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[479], w_fp[304], w_fp[7], COUPs[2], 1.0, 0., 0., w_fp[490] );
    // Amplitude(s) for diagram number 5375
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[4], w_fp[586], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[4], w_fp[585], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[4], w_fp[490], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5376 OF 15495 ***
    // Wavefunction(s) for diagram number 5376
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[479], w_fp[9], w_fp[4], COUPs[2], 1.0, 0., 0., w_fp[45] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[479], w_fp[9], w_fp[4], COUPs[2], 1.0, 0., 0., w_fp[526] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[479], w_fp[9], w_fp[4], COUPs[2], 1.0, 0., 0., w_fp[525] );
    // Amplitude(s) for diagram number 5376
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[304], w_fp[7], w_fp[45], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[304], w_fp[7], w_fp[526], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[304], w_fp[7], w_fp[525], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5377 OF 15495 ***
    // Wavefunction(s) for diagram number 5377
    // (none)
    // Amplitude(s) for diagram number 5377
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[304], w_fp[4], w_fp[452], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[304], w_fp[4], w_fp[444], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[304], w_fp[4], w_fp[509], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5378 OF 15495 ***
    // Wavefunction(s) for diagram number 5378
    // (none)
    // Amplitude(s) for diagram number 5378
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[304], w_fp[9], w_fp[545], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[304], w_fp[9], w_fp[528], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[304], w_fp[9], w_fp[447], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5379 OF 15495 ***
    // Wavefunction(s) for diagram number 5379
    // (none)
    // Amplitude(s) for diagram number 5379
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[304], w_fp[9], w_fp[102], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[304], w_fp[9], w_fp[102], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[304], w_fp[9], w_fp[102], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5380 OF 15495 ***
    // Wavefunction(s) for diagram number 5380
    // (none)
    // Amplitude(s) for diagram number 5380
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[102], w_fp[521], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5381 OF 15495 ***
    // Wavefunction(s) for diagram number 5381
    // (none)
    // Amplitude(s) for diagram number 5381
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[304], w_fp[102], w_fp[474], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5382 OF 15495 ***
    // Wavefunction(s) for diagram number 5382
    // (none)
    // Amplitude(s) for diagram number 5382
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[304], w_fp[9], w_fp[579], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5383 OF 15495 ***
    // Wavefunction(s) for diagram number 5383
    // (none)
    // Amplitude(s) for diagram number 5383
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[521], w_fp[164], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[246] += amp_sv[0];
    jamp_sv[248] -= amp_sv[0];
    jamp_sv[276] -= amp_sv[0];
    jamp_sv[278] += amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[290] += amp_sv[0];
    jamp_sv[318] += amp_sv[0];
    jamp_sv[320] -= amp_sv[0];
    jamp_sv[338] -= amp_sv[0];
    jamp_sv[339] += amp_sv[0];
    jamp_sv[346] += amp_sv[0];
    jamp_sv[347] -= amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[349] -= amp_sv[0];
    jamp_sv[356] -= amp_sv[0];
    jamp_sv[357] += amp_sv[0];

    // *** DIAGRAM 5384 OF 15495 ***
    // Wavefunction(s) for diagram number 5384
    // (none)
    // Amplitude(s) for diagram number 5384
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[163], w_fp[521], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5385 OF 15495 ***
    // Wavefunction(s) for diagram number 5385
    // (none)
    // Amplitude(s) for diagram number 5385
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[324], w_fp[515], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[250] += amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[293] += amp_sv[0];

    // *** DIAGRAM 5386 OF 15495 ***
    // Wavefunction(s) for diagram number 5386
    // (none)
    // Amplitude(s) for diagram number 5386
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[515], w_fp[314], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5387 OF 15495 ***
    // Wavefunction(s) for diagram number 5387
    // (none)
    // Amplitude(s) for diagram number 5387
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[530], w_fp[502], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[280] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[323] += amp_sv[0];

    // *** DIAGRAM 5388 OF 15495 ***
    // Wavefunction(s) for diagram number 5388
    // (none)
    // Amplitude(s) for diagram number 5388
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[530], w_fp[156], w_fp[314], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5389 OF 15495 ***
    // Wavefunction(s) for diagram number 5389
    // (none)
    // Amplitude(s) for diagram number 5389
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[530], w_fp[163], w_fp[304], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[346] += amp_sv[0];
    jamp_sv[347] -= amp_sv[0];
    jamp_sv[356] -= amp_sv[0];
    jamp_sv[357] += amp_sv[0];

    // *** DIAGRAM 5390 OF 15495 ***
    // Wavefunction(s) for diagram number 5390
    // (none)
    // Amplitude(s) for diagram number 5390
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[502], w_fp[488], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[276] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5391 OF 15495 ***
    // Wavefunction(s) for diagram number 5391
    // (none)
    // Amplitude(s) for diagram number 5391
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[324], w_fp[156], w_fp[488], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5392 OF 15495 ***
    // Wavefunction(s) for diagram number 5392
    // (none)
    // Amplitude(s) for diagram number 5392
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[304], w_fp[164], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[276] += amp_sv[0];
    jamp_sv[278] -= amp_sv[0];
    jamp_sv[280] -= amp_sv[0];
    jamp_sv[281] += amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[318] -= amp_sv[0];
    jamp_sv[320] += amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[323] -= amp_sv[0];
    jamp_sv[338] += amp_sv[0];
    jamp_sv[339] -= amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[349] += amp_sv[0];

    // *** DIAGRAM 5393 OF 15495 ***
    // Wavefunction(s) for diagram number 5393
    // (none)
    // Amplitude(s) for diagram number 5393
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[324], w_fp[163], w_fp[479], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[338] += amp_sv[0];
    jamp_sv[339] -= amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[349] += amp_sv[0];

    // *** DIAGRAM 5394 OF 15495 ***
    // Wavefunction(s) for diagram number 5394
    // (none)
    // Amplitude(s) for diagram number 5394
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[479], w_fp[314], w_fp[164], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[246] -= amp_sv[0];
    jamp_sv[248] += amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[288] += amp_sv[0];
    jamp_sv[290] -= amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[293] += amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[346] -= amp_sv[0];
    jamp_sv[347] += amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[357] -= amp_sv[0];

    // *** DIAGRAM 5395 OF 15495 ***
    // Wavefunction(s) for diagram number 5395
    // (none)
    // Amplitude(s) for diagram number 5395
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[156], w_fp[586], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[246] -= amp_sv[0];
    jamp_sv[248] += amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[288] += amp_sv[0];
    jamp_sv[290] -= amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[293] += amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[346] -= amp_sv[0];
    jamp_sv[347] += amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[156], w_fp[585], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[250] += amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[276] -= amp_sv[0];
    jamp_sv[278] += amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[293] += amp_sv[0];
    jamp_sv[318] += amp_sv[0];
    jamp_sv[320] -= amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[338] -= amp_sv[0];
    jamp_sv[339] += amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[349] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[156], w_fp[490], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[246] += amp_sv[0];
    jamp_sv[248] -= amp_sv[0];
    jamp_sv[276] -= amp_sv[0];
    jamp_sv[278] += amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[290] += amp_sv[0];
    jamp_sv[318] += amp_sv[0];
    jamp_sv[320] -= amp_sv[0];
    jamp_sv[338] -= amp_sv[0];
    jamp_sv[339] += amp_sv[0];
    jamp_sv[346] += amp_sv[0];
    jamp_sv[347] -= amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[349] -= amp_sv[0];
    jamp_sv[356] -= amp_sv[0];
    jamp_sv[357] += amp_sv[0];

    // *** DIAGRAM 5396 OF 15495 ***
    // Wavefunction(s) for diagram number 5396
    // (none)
    // Amplitude(s) for diagram number 5396
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[179], w_fp[156], w_fp[521], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5397 OF 15495 ***
    // Wavefunction(s) for diagram number 5397
    // (none)
    // Amplitude(s) for diagram number 5397
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[179], w_fp[515], w_fp[304], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[246] += amp_sv[0];
    jamp_sv[248] -= amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[290] += amp_sv[0];

    // *** DIAGRAM 5398 OF 15495 ***
    // Wavefunction(s) for diagram number 5398
    // (none)
    // Amplitude(s) for diagram number 5398
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[536], w_fp[156], w_fp[304], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[276] += amp_sv[0];
    jamp_sv[278] -= amp_sv[0];
    jamp_sv[318] -= amp_sv[0];
    jamp_sv[320] += amp_sv[0];

    // *** DIAGRAM 5399 OF 15495 ***
    // Wavefunction(s) for diagram number 5399
    // (none)
    // Amplitude(s) for diagram number 5399
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[521], w_fp[228], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[613] += amp_sv[0];
    jamp_sv[616] -= amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[658] += amp_sv[0];
    jamp_sv[659] -= amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    jamp_sv[669] += amp_sv[0];
    jamp_sv[673] -= amp_sv[0];
    jamp_sv[676] += amp_sv[0];
    jamp_sv[703] += amp_sv[0];
    jamp_sv[706] -= amp_sv[0];

    // *** DIAGRAM 5400 OF 15495 ***
    // Wavefunction(s) for diagram number 5400
    // (none)
    // Amplitude(s) for diagram number 5400
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[225], w_fp[521], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] += cxtype( 0, 1 ) * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 45 );
    storeWf( wfs, w_cx, nevt, 480 );
    storeWf( wfs, w_cx, nevt, 481 );
    storeWf( wfs, w_cx, nevt, 490 );
    storeWf( wfs, w_cx, nevt, 510 );
    storeWf( wfs, w_cx, nevt, 521 );
    storeWf( wfs, w_cx, nevt, 525 );
    storeWf( wfs, w_cx, nevt, 526 );
    storeWf( wfs, w_cx, nevt, 585 );
    storeWf( wfs, w_cx, nevt, 586 );
#endif
#endif
  }

  //--------------------------------------------------------------------------
}
