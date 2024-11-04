// Copyright (C) 2010 The MadGraph5_aMC@NLO development team and contributors.
// Created by: J. Alwall (Oct 2010) for the MG5aMC CPP backend.
//==========================================================================
// Copyright (C) 2020-2024 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Modified by: S. Roiser (Feb 2020) for the MG5aMC CUDACPP plugin.
// Further modified by: S. Hageboeck, O. Mattelaer, S. Roiser, J. Teig, A. Valassi, Z. Wettersten (2020-2024) for the MG5aMC CUDACPP plugin.
//==========================================================================
// This file has been automatically generated for CUDA/C++ standalone by
// MadGraph5_aMC@NLO v. 3.6.0, 2024-09-30
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "CPPProcess.h"

#include "mgOnGpuConfig.h"

#include "HelAmps_MSSM_SLHA2.h"
#include "MemoryAccessAmplitudes.h"
#include "MemoryAccessChannelIds.h"
#include "MemoryAccessCouplings.h"
#include "MemoryAccessCouplingsFixed.h"
#include "MemoryAccessGs.h"
#include "MemoryAccessMatrixElements.h"
#include "MemoryAccessMomenta.h"
#include "MemoryAccessWavefunctions.h"

#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
#include "MemoryAccessDenominators.h"
#include "MemoryAccessNumerators.h"
#include "coloramps.h"
#endif

#include <algorithm>
#include <array>
#include <cfenv>  // for feenableexcept, fegetexcept and FE_XXX
#include <cfloat> // for FLT_MIN
#include <cstring>
#include <iomanip>
#include <iostream>
#include <memory>

// Test ncu metrics for CUDA thread divergence
#undef MGONGPU_TEST_DIVERGENCE
//#define MGONGPU_TEST_DIVERGENCE 1

//--------------------------------------------------------------------------

// Enable FPE traps (see #701, #733, #831 - except on MacOS where feenableexcept is not defined #730)
// [NB1: Fortran default is -ffpe-trap=none, i.e. FPE traps are not enabled, https://gcc.gnu.org/onlinedocs/gfortran/Debugging-Options.html]
// [NB2: Fortran default is -ffpe-summary=invalid,zero,overflow,underflow,denormal, i.e. warn at the end on STOP]
inline void
fpeEnable()
{
  static bool first = true; // FIXME: quick and dirty hack to do this only once (can be removed when separate C++/CUDA builds are implemented)
  if( !first ) return;
  first = false;
#ifndef __APPLE__ // on MacOS feenableexcept is not defined #730
  //int fpes = fegetexcept();
  //std::cout << "fpeEnable: analyse fegetexcept()=" << fpes << std::endl;
  //std::cout << "fpeEnable:     FE_DIVBYZERO is" << ( ( fpes & FE_DIVBYZERO ) ? " " : " NOT " ) << "enabled" << std::endl;
  //std::cout << "fpeEnable:     FE_INEXACT is" << ( ( fpes & FE_INEXACT ) ? " " : " NOT " ) << "enabled" << std::endl;
  //std::cout << "fpeEnable:     FE_INVALID is" << ( ( fpes & FE_INVALID ) ? " " : " NOT " ) << "enabled" << std::endl;
  //std::cout << "fpeEnable:     FE_OVERFLOW is" << ( ( fpes & FE_OVERFLOW ) ? " " : " NOT " ) << "enabled" << std::endl;
  //std::cout << "fpeEnable:     FE_UNDERFLOW is" << ( ( fpes & FE_UNDERFLOW ) ? " " : " NOT " ) << "enabled" << std::endl;
  constexpr bool enableFPE = true; // this is hardcoded and no longer controlled by getenv( "CUDACPP_RUNTIME_ENABLEFPE" )
  if( enableFPE )
  {
    std::cout << "INFO: The following Floating Point Exceptions will cause SIGFPE program aborts: FE_DIVBYZERO, FE_INVALID, FE_OVERFLOW" << std::endl;
    feenableexcept( FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW ); // new strategy #831 (do not enable FE_UNDERFLOW)
    //fpes = fegetexcept();
    //std::cout << "fpeEnable: analyse fegetexcept()=" << fpes << std::endl;
    //std::cout << "fpeEnable:     FE_DIVBYZERO is" << ( ( fpes & FE_DIVBYZERO ) ? " " : " NOT " ) << "enabled" << std::endl;
    //std::cout << "fpeEnable:     FE_INEXACT is" << ( ( fpes & FE_INEXACT ) ? " " : " NOT " ) << "enabled" << std::endl;
    //std::cout << "fpeEnable:     FE_INVALID is" << ( ( fpes & FE_INVALID ) ? " " : " NOT " ) << "enabled" << std::endl;
    //std::cout << "fpeEnable:     FE_OVERFLOW is" << ( ( fpes & FE_OVERFLOW ) ? " " : " NOT " ) << "enabled" << std::endl;
    //std::cout << "fpeEnable:     FE_UNDERFLOW is" << ( ( fpes & FE_UNDERFLOW ) ? " " : " NOT " ) << "enabled" << std::endl;
  }
  else
  {
    //std::cout << "INFO: Do not enable SIGFPE traps for Floating Point Exceptions" << std::endl;
  }
#else
  //std::cout << "INFO: Keep default SIGFPE settings because feenableexcept is not available on MacOS" << std::endl;
#endif
}

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: g g > t1 t1~ WEIGHTED<=2 @1

#ifdef MGONGPUCPP_GPUIMPL
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
  constexpr int nw6 = CPPProcess::nw6;       // dimensions of each wavefunction (HELAS KEK 91-11): e.g. 6 for e+ e- -> mu+ mu- (fermions and vectors)
  constexpr int npar = CPPProcess::npar;     // #particles in total (external = initial + final): e.g. 4 for e+ e- -> mu+ mu-
  constexpr int ncomb = CPPProcess::ncomb;   // #helicity combinations: e.g. 16 for e+ e- -> mu+ mu- (2**4 = fermion spin up/down ** npar)
  constexpr int ncolor = CPPProcess::ncolor; // the number of leading colors

  // [NB: I am currently unable to get the right value of nwf in CPPProcess.h - will hardcode it in CPPProcess.cc instead (#644)]
  //using CPPProcess::nwf; // #wavefunctions = #external (npar) + #internal: e.g. 5 for e+ e- -> mu+ mu- (1 internal is gamma or Z)

  using Parameters_MSSM_SLHA2_dependentCouplings::ndcoup;   // #couplings that vary event by event (depend on running alphas QCD)
  using Parameters_MSSM_SLHA2_independentCouplings::nicoup; // #couplings that are fixed for all events (do not depend on running alphas QCD)

  // The number of SIMD vectors of events processed by calculate_wavefunction
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
  constexpr int nParity = 2;
#else
  constexpr int nParity = 1;
#endif

  // Physics parameters (masses, coupling, etc...)
  // For CUDA performance, hardcoded constexpr's would be better: fewer registers and a tiny throughput increase
  // However, physics parameters are user-defined through card files: use CUDA constant memory instead (issue #39)
  // [NB if hardcoded parameters are used, it's better to define them here to avoid silent shadowing (issue #263)]
  constexpr int nIPD = 4; // SM independent parameters used in this CPPProcess.cc (FIXME? rename as sm_IndepParam?)
  // Note: in the Python code generator, nIPD == nparam, while nIPC <= nicoup, because (see #823)
  // nIPC may vary from one P*/CPPProcess.cc to another, while nicoup is defined in src/Param.h and is common to all P*
  constexpr int nIPC = 0; // SM independent couplings used in this CPPProcess.cc (FIXME? rename as sm_IndepCoupl?)
  static_assert( nIPC <= nicoup );
  static_assert( nIPD >= 0 ); // Hack to avoid build warnings when nIPD==0 is unused
  static_assert( nIPC >= 0 ); // Hack to avoid build warnings when nIPC==0 is unused
#ifdef MGONGPU_HARDCODE_PARAM
  __device__ const fptype cIPD[nIPD] = { (fptype)Parameters_MSSM_SLHA2::mdl_Msu3, (fptype)Parameters_MSSM_SLHA2::mdl_Wsu3, (fptype)Parameters_MSSM_SLHA2::mdl_Msu6, (fptype)Parameters_MSSM_SLHA2::mdl_Wsu6 };
  __device__ const fptype* cIPC = nullptr; // unused as nIPC=0
#else
#ifdef MGONGPUCPP_GPUIMPL
  __device__ __constant__ fptype cIPD[nIPD];
  __device__ __constant__ fptype* cIPC = nullptr; // unused as nIPC=0
#else
  static fptype cIPD[nIPD];
  static fptype* cIPC = nullptr; // unused as nIPC=0
#endif
#endif

  // AV Jan 2024 (PR #625): this ugly #define was the only way I found to avoid creating arrays[nBsm] in CPPProcess.cc if nBsm is 0
  // The problem is that nBsm is determined when generating Parameters.h, which happens after CPPProcess.cc has already been generated
  // For simplicity, keep this code hardcoded also for SM processes (a nullptr is needed as in the case nBsm == 0)
#ifdef MGONGPUCPP_NBSMINDEPPARAM_GT_0
#ifdef MGONGPU_HARDCODE_PARAM
  __device__ const double* bsmIndepParam = Parameters_MSSM_SLHA2::mdl_bsmIndepParam;
#else
#ifdef MGONGPUCPP_GPUIMPL
  __device__ __constant__ double bsmIndepParam[Parameters_MSSM_SLHA2::nBsmIndepParam];
#else
  static double bsmIndepParam[Parameters_MSSM_SLHA2::nBsmIndepParam];
#endif
#endif
#else
#ifdef MGONGPU_HARDCODE_PARAM
  __device__ const double* bsmIndepParam = nullptr;
#else
#ifdef MGONGPUCPP_GPUIMPL
  __device__ __constant__ double* bsmIndepParam = nullptr;
#else
  static double* bsmIndepParam = nullptr;
#endif
#endif
#endif

  // Helicity combinations (and filtering of "good" helicity combinations)
#ifdef MGONGPUCPP_GPUIMPL
  __device__ __constant__ short cHel[ncomb][npar];
  __device__ __constant__ int dcNGoodHel;
  __device__ __constant__ int dcGoodHel[ncomb];
#else
  static short cHel[ncomb][npar];
#endif
  static int cNGoodHel;
  static int cGoodHel[ncomb];

  //--------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL
  class DeviceAccessJamp2
  {
  public:
    static __device__ inline fptype&
    kernelAccessIcol( fptype* buffer, const int icol )
    {
      const int nevt = gridDim.x * blockDim.x;
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x;
      return buffer[icol * nevt + ievt];
    }
    static __device__ inline const fptype&
    kernelAccessIcolConst( const fptype* buffer, const int icol )
    {
      const int nevt = gridDim.x * blockDim.x;
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x;
      return buffer[icol * nevt + ievt];
    }
  };
#endif

  //--------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL
  __device__ INLINE unsigned int
  gpu_channelId( const unsigned int* allChannelIds )
  {
    unsigned int channelId = 0; // disable multichannel single-diagram enhancement unless allChannelIds != nullptr
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    using CID_ACCESS = DeviceAccessChannelIds; // non-trivial access: buffer includes all events
    // SCALAR channelId for the current event (CUDA)
    if( allChannelIds != nullptr )
    {
      const unsigned int* channelIds = allChannelIds;                            // fix #899 (distinguish channelIds and allChannelIds)
      const uint_sv channelIds_sv = CID_ACCESS::kernelAccessConst( channelIds ); // fix #895 (compute this only once for all diagrams)
      // NB: channelIds_sv is a scalar in CUDA
      channelId = channelIds_sv;
      assert( channelId > 0 ); // SANITY CHECK: scalar channelId must be > 0 if multichannel is enabled (allChannelIds != nullptr)
    }
#endif
    return channelId;
  }
#endif

  //--------------------------------------------------------------------------

  // Evaluate |M|^2 for each subprocess
  // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum of |M|^2 over helicities for the given event(s)
  // (similarly, it also ADDS the numerator and denominator for a given ihel to their running sums over helicities)
  // In CUDA, this function computes the ME for a single event
  // ** NB1: NEW Nov2024! In CUDA this is now a kernel function (it used to be a device function)
  // ** NB2: NEW Nov2024! in CUDA this now takes a channelId array as input (it used to take a scalar channelId as input)
  // In C++, this function computes the ME for a single event "page" or SIMD vector (or for two in "mixed" precision mode, nParity=2)
  // *** NB: in C++, calculate_wavefunction accepts a SCALAR channelId because it is GUARANTEED that all events in a SIMD vector have the same channelId #898
  __global__ INLINE void /* clang-format off */
  calculate_wavefunctions( int ihel,
                           const fptype* allmomenta,          // input: momenta[nevt*npar*4]
                           const fptype* allcouplings,        // input: couplings[nevt*ndcoup*2]
                           fptype* allMEs,                    // output: allMEs[nevt], |M|^2 running_sum_over_helicities
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
#ifdef MGONGPUCPP_GPUIMPL
                           const unsigned int* allChannelIds, // input: multichannel channelIds[nevt] (1 to #diagrams); nullptr to disable SDE enhancement (fix #899/#911)
#else
                           const unsigned int channelId,      // input: multichannel SCALAR channelId (1 to #diagrams, 0 to disable SDE) for this event or SIMD vector
#endif
                           fptype* allNumerators,             // output: multichannel numerators[nevt], running_sum_over_helicities
                           fptype* allDenominators,           // output: multichannel denominators[nevt], running_sum_over_helicities
#endif
#ifdef MGONGPUCPP_GPUIMPL
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
                           fptype_sv* allJamp2s,              // output: jamp2[ncolor][nevt] for color choice (nullptr if disabled)
#endif
                           const int nevt                     // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
#else
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
                           fptype_sv* jamp2_sv,               // output: jamp2[nParity][ncolor][neppV] for color choice (nullptr if disabled)
#endif
                           const int ievt00                   // input: first event number in current C++ event page (for CUDA, ievt depends on threadid)
#endif
                           )
  //ALWAYS_INLINE // attributes are not permitted in a function definition
  {
#ifdef MGONGPUCPP_GPUIMPL
    using namespace mg5amcGpu;
    using M_ACCESS = DeviceAccessMomenta;         // non-trivial access: buffer includes all events
    using E_ACCESS = DeviceAccessMatrixElements;  // non-trivial access: buffer includes all events
    using W_ACCESS = DeviceAccessWavefunctions;   // TRIVIAL ACCESS (no kernel splitting yet): buffer for one event
    using A_ACCESS = DeviceAccessAmplitudes;      // TRIVIAL ACCESS (no kernel splitting yet): buffer for one event
    using CD_ACCESS = DeviceAccessCouplings;      // non-trivial access (dependent couplings): buffer includes all events
    using CI_ACCESS = DeviceAccessCouplingsFixed; // TRIVIAL access (independent couplings): buffer for one event
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    using NUM_ACCESS = DeviceAccessNumerators;    // non-trivial access: buffer includes all events
    using DEN_ACCESS = DeviceAccessDenominators;  // non-trivial access: buffer includes all events
#endif
#else
    using namespace mg5amcCpu;
    using M_ACCESS = HostAccessMomenta;         // non-trivial access: buffer includes all events
    using E_ACCESS = HostAccessMatrixElements;  // non-trivial access: buffer includes all events
    using W_ACCESS = HostAccessWavefunctions;   // TRIVIAL ACCESS (no kernel splitting yet): buffer for one event
    using A_ACCESS = HostAccessAmplitudes;      // TRIVIAL ACCESS (no kernel splitting yet): buffer for one event
    using CD_ACCESS = HostAccessCouplings;      // non-trivial access (dependent couplings): buffer includes all events
    using CI_ACCESS = HostAccessCouplingsFixed; // TRIVIAL access (independent couplings): buffer for one event
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    using NUM_ACCESS = HostAccessNumerators;    // non-trivial access: buffer includes all events
    using DEN_ACCESS = HostAccessDenominators;  // non-trivial access: buffer includes all events
#endif
#endif
    mgDebug( 0, __FUNCTION__ );
    //bool debug = true;
#ifndef MGONGPUCPP_GPUIMPL
    //debug = ( ievt00 >= 64 && ievt00 < 80 && ihel == 3 ); // example: debug #831
    //if( debug ) printf( "calculate_wavefunctions: ievt00=%d ihel=%2d\n", ievt00, ihel );
#else
    //const int ievt = blockDim.x * blockIdx.x + threadIdx.x;
    //debug = ( ievt == 0 );
    //if( debug ) printf( "calculate_wavefunctions: ievt=%6d ihel=%2d\n", ievt, ihel );
#endif /* clang-format on */

    // The variable nwf (which is specific to each P1 subdirectory, #644) is only used here
    // It is hardcoded here because various attempts to hardcode it in CPPProcess.h at generation time gave the wrong result...
    static const int nwf = 5; // #wavefunctions = #external (npar) + #internal: e.g. 5 for e+ e- -> mu+ mu- (1 internal is gamma or Z)

    // Local TEMPORARY variables for a subset of Feynman diagrams in the given CUDA event (ievt) or C++ event page (ipagV)
    // [NB these variables are reused several times (and re-initialised each time) within the same event or event page]
    // ** NB: in other words, amplitudes and wavefunctions still have TRIVIAL ACCESS: there is currently no need
    // ** NB: to have large memory structurs for wavefunctions/amplitudes in all events (no kernel splitting yet)!
    //MemoryBufferWavefunctions w_buffer[nwf]{ neppV };
    cxtype_sv w_sv[nwf][nw6]; // particle wavefunctions within Feynman diagrams (nw6 is often 6, the dimension of spin 1/2 or spin 1 particles)
    cxtype_sv amp_sv[1];      // invariant amplitude for one given Feynman diagram

    // Proof of concept for using fptype* in the interface
    fptype* w_fp[nwf];
    for( int iwf = 0; iwf < nwf; iwf++ ) w_fp[iwf] = reinterpret_cast<fptype*>( w_sv[iwf] );
    fptype* amp_fp;
    amp_fp = reinterpret_cast<fptype*>( amp_sv );

    // Local variables for the given CUDA event (ievt) or C++ event page (ipagV)
    // [jamp: sum (for one event or event page) of the invariant amplitudes for all Feynman diagrams in a given color combination]
    cxtype_sv jamp_sv[ncolor] = {}; // all zeros (NB: vector cxtype_v IS initialized to 0, but scalar cxtype is NOT, if "= {}" is missing!)

    // === Calculate wavefunctions and amplitudes for all diagrams in all processes         ===
    // === (for one event in CUDA, for one - or two in mixed mode - SIMD event pages in C++ ===
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
    // Mixed fptypes #537: float for color algebra and double elsewhere
    // Delay color algebra and ME updates (only on even pages)
    cxtype_sv jamp_sv_previous[ncolor] = {};
    fptype* MEs_previous = 0;
#endif
    for( int iParity = 0; iParity < nParity; ++iParity )
    { // START LOOP ON IPARITY
#ifndef MGONGPUCPP_GPUIMPL
      const int ievt0 = ievt00 + iParity * neppV;
#endif
      //constexpr size_t nxcoup = ndcoup + nicoup; // both dependent and independent couplings (BUG #823)
      constexpr size_t nxcoup = ndcoup + nIPC; // both dependent and independent couplings (FIX #823)
      const fptype* allCOUPs[nxcoup];
#ifdef __CUDACC__ // this must be __CUDACC__ (not MGONGPUCPP_GPUIMPL)
#pragma nv_diagnostic push
#pragma nv_diag_suppress 186 // e.g. <<warning #186-D: pointless comparison of unsigned integer with zero>>
#endif
      for( size_t idcoup = 0; idcoup < ndcoup; idcoup++ )
        allCOUPs[idcoup] = CD_ACCESS::idcoupAccessBufferConst( allcouplings, idcoup ); // dependent couplings, vary event-by-event
      //for( size_t iicoup = 0; iicoup < nicoup; iicoup++ )                             // BUG #823
      for( size_t iicoup = 0; iicoup < nIPC; iicoup++ )                                 // FIX #823
        allCOUPs[ndcoup + iicoup] = CI_ACCESS::iicoupAccessBufferConst( cIPC, iicoup ); // independent couplings, fixed for all events
#ifdef MGONGPUCPP_GPUIMPL
#ifdef __CUDACC__ // this must be __CUDACC__ (not MGONGPUCPP_GPUIMPL)
#pragma nv_diagnostic pop
#endif
      // CUDA kernels take input/output buffers with momenta/MEs for all events
      const fptype* momenta = allmomenta;
      const fptype* COUPs[nxcoup];
      for( size_t ixcoup = 0; ixcoup < nxcoup; ixcoup++ ) COUPs[ixcoup] = allCOUPs[ixcoup];
      fptype* MEs = allMEs;
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      fptype* numerators = allNumerators;
      fptype* denominators = allDenominators;
#endif
#else
      // C++ kernels take input/output buffers with momenta/MEs for one specific event (the first in the current event page)
      const fptype* momenta = M_ACCESS::ieventAccessRecordConst( allmomenta, ievt0 );
      const fptype* COUPs[nxcoup];
      for( size_t idcoup = 0; idcoup < ndcoup; idcoup++ )
        COUPs[idcoup] = CD_ACCESS::ieventAccessRecordConst( allCOUPs[idcoup], ievt0 ); // dependent couplings, vary event-by-event
      //for( size_t iicoup = 0; iicoup < nicoup; iicoup++ ) // BUG #823
      for( size_t iicoup = 0; iicoup < nIPC; iicoup++ )     // FIX #823
        COUPs[ndcoup + iicoup] = allCOUPs[ndcoup + iicoup]; // independent couplings, fixed for all events
      fptype* MEs = E_ACCESS::ieventAccessRecord( allMEs, ievt0 );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      fptype* numerators = NUM_ACCESS::ieventAccessRecord( allNumerators, ievt0 );
      fptype* denominators = DEN_ACCESS::ieventAccessRecord( allDenominators, ievt0 );
#endif
#endif

      // Reset color flows (reset jamp_sv) at the beginning of a new event or event page
      for( int i = 0; i < ncolor; i++ ) { jamp_sv[i] = cxzero_sv(); }

#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
#ifdef MGONGPUCPP_GPUIMPL
      // SCALAR channelId for the current event (CUDA)
      unsigned int channelId = gpu_channelId( allChannelIds );
#endif
      // Numerators and denominators for the current event (CUDA) or SIMD event page (C++)
      fptype_sv& numerators_sv = NUM_ACCESS::kernelAccess( numerators );
      fptype_sv& denominators_sv = DEN_ACCESS::kernelAccess( denominators );
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

#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // *** COLOR CHOICE BELOW ***
      // Store the leading color flows for choice of color
#ifndef MGONGPUCPP_GPUIMPL
      if( jamp2_sv ) // disable color choice if nullptr
        for( int icol = 0; icol < ncolor; icol++ )
          jamp2_sv[ncolor * iParity + icol] += cxabs2( jamp_sv[icol] ); // may underflow #831
#else
      assert( iParity == 0 ); // sanity check for J_ACCESS
      using J_ACCESS = DeviceAccessJamp2;
      if( allJamp2s ) // disable color choice if nullptr
        for( int icol = 0; icol < ncolor; icol++ )
          J_ACCESS::kernelAccessIcol( allJamp2s, icol ) += cxabs2( jamp_sv[icol] );
#endif
#endif

      // *** COLOR MATRIX BELOW ***
      // (This method used to be called CPPProcess::matrix_1_gg_t1t1x()?)

      // The color denominators (initialize all array elements, with ncolor=2)
      // [NB do keep 'static' for these constexpr arrays, see issue #283]
      static constexpr fptype2 denom[ncolor] = { 3, 3 }; // 1-D array[2]

      // The color matrix (initialize all array elements, with ncolor=2)
      // [NB do keep 'static' for these constexpr arrays, see issue #283]
      static constexpr fptype2 cf[ncolor][ncolor] = {
        { 16, -2 },
        { -2, 16 } }; // 2-D array[2][2]

#ifndef MGONGPUCPP_GPUIMPL
      // Pre-compute a constexpr triangular color matrix properly normalized #475
      struct TriangularNormalizedColorMatrix
      {
        // See https://stackoverflow.com/a/34465458
        __host__ __device__ constexpr TriangularNormalizedColorMatrix()
          : value()
        {
          for( int icol = 0; icol < ncolor; icol++ )
          {
            // Diagonal terms
            value[icol][icol] = cf[icol][icol] / denom[icol];
            // Off-diagonal terms
            for( int jcol = icol + 1; jcol < ncolor; jcol++ )
              value[icol][jcol] = 2 * cf[icol][jcol] / denom[icol];
          }
        }
        fptype2 value[ncolor][ncolor];
      };
      static constexpr auto cf2 = TriangularNormalizedColorMatrix();
#endif

#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
      if( iParity == 0 ) // NB: first page is 0! skip even pages, compute on odd pages
      {
        // Mixed fptypes: delay color algebra and ME updates to next (odd) ipagV
        for( int icol = 0; icol < ncolor; icol++ )
          jamp_sv_previous[icol] = jamp_sv[icol];
        MEs_previous = MEs;
        continue; // go to next iParity in the loop: skip color algebra and ME update on odd pages
      }
      fptype_sv deltaMEs_previous = { 0 };
#endif

      // Sum and square the color flows to get the matrix element
      // (compute |M|^2 by squaring |M|, taking into account colours)
      // Sum and square the color flows to get the matrix element
      // (compute |M|^2 by squaring |M|, taking into account colours)
      fptype_sv deltaMEs = { 0 }; // all zeros https://en.cppreference.com/w/c/language/array_initialization#Notes

      // Use the property that M is a real matrix (see #475):
      // we can rewrite the quadratic form (A-iB)(M)(A+iB) as AMA - iBMA + iBMA + BMB = AMA + BMB
      // In addition, on C++ use the property that M is symmetric (see #475),
      // and also use constexpr to compute "2*" and "/denom[icol]" once and for all at compile time:
      // we gain (not a factor 2...) in speed here as we only loop over the up diagonal part of the matrix.
      // Strangely, CUDA is slower instead, so keep the old implementation for the moment.
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
      fptype2_sv jampR_sv[ncolor] = { 0 };
      fptype2_sv jampI_sv[ncolor] = { 0 };
      for( int icol = 0; icol < ncolor; icol++ )
      {
        jampR_sv[icol] = fpvmerge( cxreal( jamp_sv_previous[icol] ), cxreal( jamp_sv[icol] ) );
        jampI_sv[icol] = fpvmerge( cximag( jamp_sv_previous[icol] ), cximag( jamp_sv[icol] ) );
      }
#endif
      for( int icol = 0; icol < ncolor; icol++ )
      {
        //if( debug ) printf( "calculate_wavefunctions... icol=%d\n", icol );
#ifndef MGONGPUCPP_GPUIMPL
        // === C++ START ===
        // Diagonal terms
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
        fptype2_sv& jampRi_sv = jampR_sv[icol];
        fptype2_sv& jampIi_sv = jampI_sv[icol];
#else
        fptype2_sv jampRi_sv = (fptype2_sv)( cxreal( jamp_sv[icol] ) );
        fptype2_sv jampIi_sv = (fptype2_sv)( cximag( jamp_sv[icol] ) );
#endif
        fptype2_sv ztempR_sv = cf2.value[icol][icol] * jampRi_sv;
        fptype2_sv ztempI_sv = cf2.value[icol][icol] * jampIi_sv;
        // Off-diagonal terms
        for( int jcol = icol + 1; jcol < ncolor; jcol++ )
        {
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
          fptype2_sv& jampRj_sv = jampR_sv[jcol];
          fptype2_sv& jampIj_sv = jampI_sv[jcol];
#else
          fptype2_sv jampRj_sv = (fptype2_sv)( cxreal( jamp_sv[jcol] ) );
          fptype2_sv jampIj_sv = (fptype2_sv)( cximag( jamp_sv[jcol] ) );
#endif
          ztempR_sv += cf2.value[icol][jcol] * jampRj_sv;
          ztempI_sv += cf2.value[icol][jcol] * jampIj_sv;
        }
        fptype2_sv deltaMEs2 = ( jampRi_sv * ztempR_sv + jampIi_sv * ztempI_sv ); // may underflow #831
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
        deltaMEs_previous += fpvsplit0( deltaMEs2 );
        deltaMEs += fpvsplit1( deltaMEs2 );
#else
        deltaMEs += deltaMEs2;
#endif
        // === C++ END ===
#else
        // === CUDA START ===
        fptype2_sv ztempR_sv = { 0 };
        fptype2_sv ztempI_sv = { 0 };
        for( int jcol = 0; jcol < ncolor; jcol++ )
        {
          fptype2_sv jampRj_sv = cxreal( jamp_sv[jcol] );
          fptype2_sv jampIj_sv = cximag( jamp_sv[jcol] );
          ztempR_sv += cf[icol][jcol] * jampRj_sv;
          ztempI_sv += cf[icol][jcol] * jampIj_sv;
        }
        deltaMEs += ( ztempR_sv * cxreal( jamp_sv[icol] ) + ztempI_sv * cximag( jamp_sv[icol] ) ) / denom[icol];
        // === CUDA END ===
#endif
      }

      // *** STORE THE RESULTS ***

      // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum of |M|^2 over helicities for the given event(s)
      fptype_sv& MEs_sv = E_ACCESS::kernelAccess( MEs );
      MEs_sv += deltaMEs; // fix #435
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
      fptype_sv& MEs_sv_previous = E_ACCESS::kernelAccess( MEs_previous );
      MEs_sv_previous += deltaMEs_previous;
#endif
      /*
#ifdef MGONGPUCPP_GPUIMPL
      if ( cNGoodHel > 0 ) printf( "calculate_wavefunctions: ievt=%6d ihel=%2d me_running=%f\n", blockDim.x * blockIdx.x + threadIdx.x, ihel, MEs_sv );
#else
#ifdef MGONGPU_CPPSIMD
      if( cNGoodHel > 0 )
        for( int ieppV = 0; ieppV < neppV; ieppV++ )
          printf( "calculate_wavefunctions: ievt=%6d ihel=%2d me_running=%f\n", ievt0 + ieppV, ihel, MEs_sv[ieppV] );
#else
      if ( cNGoodHel > 0 ) printf( "calculate_wavefunctions: ievt=%6d ihel=%2d me_running=%f\n", ievt0, ihel, MEs_sv );
#endif
#endif
      */
    } // END LOOP ON IPARITY
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  CPPProcess::CPPProcess( bool verbose,
                          bool debug )
    : m_verbose( verbose )
    , m_debug( debug )
#ifndef MGONGPU_HARDCODE_PARAM
    , m_pars( 0 )
#endif
    , m_masses()
  {
    // Helicities for the process [NB do keep 'static' for this constexpr array, see issue #283]
    // *** NB There is no automatic check yet that these are in the same order as Fortran! #569 ***
    static constexpr short tHel[ncomb][npar] = {
      { -1, -1, 0, 0 },
      { -1, 1, 0, 0 },
      { 1, -1, 0, 0 },
      { 1, 1, 0, 0 } };
#ifdef MGONGPUCPP_GPUIMPL
    gpuMemcpyToSymbol( cHel, tHel, ncomb * npar * sizeof( short ) );
#else
    memcpy( cHel, tHel, ncomb * npar * sizeof( short ) );
#endif
    fpeEnable(); // enable SIGFPE traps for Floating Point Exceptions
  }

  //--------------------------------------------------------------------------

  CPPProcess::~CPPProcess() {}

  //--------------------------------------------------------------------------

#ifndef MGONGPU_HARDCODE_PARAM
  // Initialize process (with parameters read from user cards)
  void
  CPPProcess::initProc( const std::string& param_card_name )
  {
    // Instantiate the model class and set parameters that stay fixed during run
    m_pars = Parameters_MSSM_SLHA2::getInstance();
    SLHAReader slha( param_card_name, m_verbose );
    m_pars->setIndependentParameters( slha );
    m_pars->setIndependentCouplings();
    //m_pars->setDependentParameters(); // now computed event-by-event (running alphas #373)
    //m_pars->setDependentCouplings(); // now computed event-by-event (running alphas #373)
    if( m_verbose )
    {
      m_pars->printIndependentParameters();
      m_pars->printIndependentCouplings();
      //m_pars->printDependentParameters(); // now computed event-by-event (running alphas #373)
      //m_pars->printDependentCouplings(); // now computed event-by-event (running alphas #373)
    }
    // Set external particle masses for this matrix element
    m_masses.push_back( m_pars->ZERO );
    m_masses.push_back( m_pars->ZERO );
    m_masses.push_back( m_pars->mdl_Msu3 );
    m_masses.push_back( m_pars->mdl_Msu3 );
    // Read physics parameters like masses and couplings from user configuration files (static: initialize once)
    // Then copy them to CUDA constant memory (issue #39) or its C++ emulation in file-scope static memory
    const fptype tIPD[nIPD] = { (fptype)m_pars->mdl_Msu3, (fptype)m_pars->mdl_Wsu3, (fptype)m_pars->mdl_Msu6, (fptype)m_pars->mdl_Wsu6 };
    //const cxtype tIPC[0] = { ... }; // nIPC=0
#ifdef MGONGPUCPP_GPUIMPL
    gpuMemcpyToSymbol( cIPD, tIPD, nIPD * sizeof( fptype ) );
    //gpuMemcpyToSymbol( cIPC, tIPC, 0 * sizeof( cxtype ) ); // nIPC=0
#ifdef MGONGPUCPP_NBSMINDEPPARAM_GT_0
    if( Parameters_MSSM_SLHA2::nBsmIndepParam > 0 )
      gpuMemcpyToSymbol( bsmIndepParam, m_pars->mdl_bsmIndepParam, Parameters_MSSM_SLHA2::nBsmIndepParam * sizeof( double ) );
#endif
#else
    memcpy( cIPD, tIPD, nIPD * sizeof( fptype ) );
    //memcpy( cIPC, tIPC, nIPC * sizeof( cxtype ) ); // nIPC=0
#ifdef MGONGPUCPP_NBSMINDEPPARAM_GT_0
    if( Parameters_MSSM_SLHA2::nBsmIndepParam > 0 )
      memcpy( bsmIndepParam, m_pars->mdl_bsmIndepParam, Parameters_MSSM_SLHA2::nBsmIndepParam * sizeof( double ) );
#endif
#endif
    //for ( int i=0; i<nIPD; i++ ) std::cout << std::setprecision(17) << "tIPD[i] = " << tIPD[i] << std::endl;
    //for ( int i=0; i<Parameters_MSSM_SLHA2::nBsmIndepParam; i++ ) std::cout << std::setprecision(17) << "m_pars->mdl_bsmIndepParam[i] = " << m_pars->mdl_bsmIndepParam[i] << std::endl;
  }
#else
  // Initialize process (with hardcoded parameters)
  void
  CPPProcess::initProc( const std::string& /*param_card_name*/ )
  {
    // Use hardcoded physics parameters
    if( m_verbose )
    {
      Parameters_MSSM_SLHA2::printIndependentParameters();
      Parameters_MSSM_SLHA2::printIndependentCouplings();
      //Parameters_MSSM_SLHA2::printDependentParameters(); // now computed event-by-event (running alphas #373)
      //Parameters_MSSM_SLHA2::printDependentCouplings(); // now computed event-by-event (running alphas #373)
    }
    // Set external particle masses for this matrix element
    m_masses.push_back( Parameters_MSSM_SLHA2::ZERO );
    m_masses.push_back( Parameters_MSSM_SLHA2::ZERO );
    m_masses.push_back( Parameters_MSSM_SLHA2::mdl_Msu3 );
    m_masses.push_back( Parameters_MSSM_SLHA2::mdl_Msu3 );
  }
#endif

  //--------------------------------------------------------------------------

  // Retrieve the compiler that was used to build this module
  const std::string
  CPPProcess::getCompiler()
  {
    std::stringstream out;
    // HIP version (HIPCC)
    // [Use __HIPCC__ instead of MGONGPUCPP_GPUIMPL here!]
    // [This tests if 'hipcc' was used even to build a .cc file, even if not necessarily 'nvcc -x cu' for a .cu file]
    // [Check 'hipcc -dM -E -x hip -I ../../src CPPProcess.cc | grep HIP']
#ifdef __HIPCC__
#if defined HIP_VERSION_MAJOR && defined HIP_VERSION_MINOR && defined HIP_VERSION_PATCH
    out << "hipcc " << HIP_VERSION_MAJOR << "." << HIP_VERSION_MINOR << "." << HIP_VERSION_PATCH;
#else
    out << "hipcc UNKNOWN";
#endif
    out << " (";
#endif
    // CUDA version (NVCC)
    // [Use __NVCC__ instead of MGONGPUCPP_GPUIMPL here!]
    // [This tests if 'nvcc' was used even to build a .cc file, even if not necessarily 'nvcc -x cu' for a .cu file]
    // [Check 'nvcc --compiler-options -dM -E dummy.c | grep CUDA': see https://stackoverflow.com/a/53713712]
#ifdef __NVCC__
#if defined __CUDACC_VER_MAJOR__ && defined __CUDACC_VER_MINOR__ && defined __CUDACC_VER_BUILD__
    out << "nvcc " << __CUDACC_VER_MAJOR__ << "." << __CUDACC_VER_MINOR__ << "." << __CUDACC_VER_BUILD__;
#else
    out << "nvcc UNKNOWN";
#endif
    out << " (";
#endif
    // ICX version (either as CXX or as host compiler inside NVCC)
#if defined __INTEL_COMPILER
#error "icc is no longer supported: please use icx"
#elif defined __INTEL_LLVM_COMPILER // alternative: __INTEL_CLANG_COMPILER
    out << "icx " << __INTEL_LLVM_COMPILER;
#ifdef __NVCC__
    out << ", ";
#else
    out << " (";
#endif
#endif
    // CLANG version (either as CXX or as host compiler inside NVCC or inside ICX)
#if defined __clang__
#if defined __clang_major__ && defined __clang_minor__ && defined __clang_patchlevel__
#ifdef __APPLE__
    out << "Apple clang " << __clang_major__ << "." << __clang_minor__ << "." << __clang_patchlevel__;
#else
    out << "clang " << __clang_major__ << "." << __clang_minor__ << "." << __clang_patchlevel__;
    /*
    // === AV 26-Jan-2024 DISABLE THIS CODE (START)
    // === AV 26-Jan-2024 First, it is totally wrong to assume that the CXX environment variable is used in the build!
    // === AV 26-Jan-2024 Second and worse, here we need build time values, while CXX in this code is evaluated at runtime!
    // GCC toolchain version inside CLANG
    std::string tchainout;
    std::string tchaincmd = "readelf -p .comment $(${CXX} -print-libgcc-file-name) |& grep 'GCC: (GNU)' | grep -v Warning | sort -u | awk '{print $5}'";
    std::unique_ptr<FILE, decltype( &pclose )> tchainpipe( popen( tchaincmd.c_str(), "r" ), pclose );
    if( !tchainpipe ) throw std::runtime_error( "`readelf ...` failed?" );
    std::array<char, 128> tchainbuf;
    while( fgets( tchainbuf.data(), tchainbuf.size(), tchainpipe.get() ) != nullptr ) tchainout += tchainbuf.data();
    tchainout.pop_back(); // remove trailing newline
#if defined __NVCC__ or defined __INTEL_LLVM_COMPILER
    out << ", gcc " << tchainout;
#else
    out << " (gcc " << tchainout << ")";
#endif
    // === AV 26-Jan-2024 DISABLE THIS CODE (END)
    */
#endif
#else
    out << "clang UNKNOWKN";
#endif
#else
    // GCC version (either as CXX or as host compiler inside NVCC)
#if defined __GNUC__ && defined __GNUC_MINOR__ && defined __GNUC_PATCHLEVEL__
    out << "gcc " << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__;
#else
    out << "gcc UNKNOWKN";
#endif
#endif
#if defined __HIPCC__ or defined __NVCC__ or defined __INTEL_LLVM_COMPILER
    out << ")";
#endif
    return out.str();
  }

  //--------------------------------------------------------------------------

  __global__ void /* clang-format off */
  computeDependentCouplings( const fptype* allgs, // input: Gs[nevt]
                             fptype* allcouplings // output: couplings[nevt*ndcoup*2]
#ifndef MGONGPUCPP_GPUIMPL
                             , const int nevt     // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
#endif
  ) /* clang-format on */
  {
#ifdef MGONGPUCPP_GPUIMPL
    using namespace mg5amcGpu;
    using G_ACCESS = DeviceAccessGs;
    using C_ACCESS = DeviceAccessCouplings;
    G2COUP<G_ACCESS, C_ACCESS>( allgs, allcouplings, bsmIndepParam );
#else
    using namespace mg5amcCpu;
    using G_ACCESS = HostAccessGs;
    using C_ACCESS = HostAccessCouplings;
    for( int ipagV = 0; ipagV < nevt / neppV; ++ipagV )
    {
      const int ievt0 = ipagV * neppV;
      const fptype* gs = MemoryAccessGs::ieventAccessRecordConst( allgs, ievt0 );
      fptype* couplings = MemoryAccessCouplings::ieventAccessRecord( allcouplings, ievt0 );
      G2COUP<G_ACCESS, C_ACCESS>( gs, couplings, bsmIndepParam );
    }
#endif
  }

  //--------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL
  void /* clang-format off */
  sigmaKin_getGoodHel( const fptype* allmomenta,   // input: momenta[nevt*npar*4]
                       const fptype* allcouplings, // input: couplings[nevt*ndcoup*2]
                       fptype* allMEs,             // output: allMEs[nevt], |M|^2 final_avg_over_helicities
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
                       fptype* allNumerators,      // output: multichannel numerators[nevt], running_sum_over_helicities
                       fptype* allDenominators,    // output: multichannel denominators[nevt], running_sum_over_helicities
#endif
                       bool* isGoodHel,            // output: isGoodHel[ncomb] - host array
                       const int nevt )            // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
  { /* clang-format on */
    const int maxtry0 = 16;
    fptype hstMEs[maxtry0];
    const int maxtry = std::min( maxtry0, nevt ); // 16, but at most nevt (avoid invalid memory access if nevt<maxtry0)
    //std::cout << "sigmaKin_getGoodHel nevt=" << nevt << " maxtry=" << maxtry << std::endl;
    for( int ihel = 0; ihel < ncomb; ihel++ )
    {
      const int gpublocks = 1;
      const int gputhreads = maxtry;
      // NEW IMPLEMENTATION OF GETGOODHEL (#630): RESET THE RUNNING SUM OVER HELICITIES TO 0 BEFORE ADDING A NEW HELICITY
      gpuMemset( allMEs, 0, maxtry * sizeof( fptype ) );
      // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum of |M|^2 over helicities for the given event(s)
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      constexpr fptype_sv* allJamp2s = nullptr;        // no need for color selection during helicity filtering
      constexpr unsigned int* allChannelIds = nullptr; // disable multichannel single-diagram enhancement
      gpuLaunchKernel( calculate_wavefunctions, gpublocks, gputhreads, ihel, allmomenta, allcouplings, allMEs, allChannelIds, allNumerators, allDenominators, allJamp2s, gpublocks * gputhreads );
#else
      gpuLaunchKernel( calculate_wavefunctions, gpublocks, gputhreads, ihel, allmomenta, allcouplings, allMEs, gpublocks * gputhreads );
#endif
      gpuMemcpy( hstMEs, allMEs, maxtry * sizeof( fptype ), gpuMemcpyDeviceToHost );
      //std::cout << "sigmaKin_getGoodHel ihel=" << ihel << std::endl;
      for( int ievt = 0; ievt < maxtry; ++ievt )
      {
        //std::cout << "sigmaKin_getGoodHel hstMEs[ievt]=" << hstMEs[ievt] << std::endl;
        if( hstMEs[ievt] != 0 ) // NEW IMPLEMENTATION OF GETGOODHEL (#630): COMPARE EACH HELICITY CONTRIBUTION TO 0
        {
          //if ( !isGoodHel[ihel] ) std::cout << "sigmaKin_getGoodHel ihel=" << ihel << " TRUE" << std::endl;
          isGoodHel[ihel] = true;
        }
      }
    }
  }
#else
  void
  sigmaKin_getGoodHel( const fptype* allmomenta,   // input: momenta[nevt*npar*4]
                       const fptype* allcouplings, // input: couplings[nevt*ndcoup*2]
                       fptype* allMEs,             // output: allMEs[nevt], |M|^2 final_avg_over_helicities
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
                       fptype* allNumerators,      // output: multichannel numerators[nevt], running_sum_over_helicities
                       fptype* allDenominators,    // output: multichannel denominators[nevt], running_sum_over_helicities
#endif
                       bool* isGoodHel,            // output: isGoodHel[ncomb] - host array
                       const int nevt )            // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
  {
    //assert( (size_t)(allmomenta) % mgOnGpu::cppAlign == 0 ); // SANITY CHECK: require SIMD-friendly alignment [COMMENT OUT TO TEST MISALIGNED ACCESS]
    //assert( (size_t)(allMEs) % mgOnGpu::cppAlign == 0 ); // SANITY CHECK: require SIMD-friendly alignment [COMMENT OUT TO TEST MISALIGNED ACCESS]
    // Allocate arrays at build time to contain at least 16 events (or at least neppV events if neppV>16, e.g. in future VPUs)
    constexpr int maxtry0 = std::max( 16, neppV ); // 16, but at least neppV (otherwise the npagV loop does not even start)
    // Loop over only nevt events if nevt is < 16 (note that nevt is always >= neppV)
    assert( nevt >= neppV );
    const int maxtry = std::min( maxtry0, nevt ); // 16, but at most nevt (avoid invalid memory access if nevt<maxtry0)
    //std::cout << "sigmaKin_getGoodHel nevt=" << nevt << " maxtry=" << maxtry << std::endl;
    // HELICITY LOOP: CALCULATE WAVEFUNCTIONS
    const int npagV = maxtry / neppV;
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT /* clang-format off */
    // Mixed fptypes #537: float for color algebra and double elsewhere
    // Delay color algebra and ME updates (only on even pages)
    assert( npagV % 2 == 0 ); // SANITY CHECK for mixed fptypes: two neppV-pages are merged to one 2*neppV-page
    const int npagV2 = npagV / 2; // loop on two SIMD pages (neppV events) at a time
#else
    const int npagV2 = npagV; // loop on one SIMD page (neppV events) at a time
#endif /* clang-format on */
    for( int ipagV2 = 0; ipagV2 < npagV2; ++ipagV2 )
    {
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT /* clang-format off */
      const int ievt00 = ipagV2 * neppV * 2; // loop on two SIMD pages (neppV events) at a time
#else
      const int ievt00 = ipagV2 * neppV; // loop on one SIMD page (neppV events) at a time
#endif /* clang-format on */
      for( int ihel = 0; ihel < ncomb; ihel++ )
      {
        //std::cout << "sigmaKin_getGoodHel ihel=" << ihel << std::endl;
        // NEW IMPLEMENTATION OF GETGOODHEL (#630): RESET THE RUNNING SUM OVER HELICITIES TO 0 BEFORE ADDING A NEW HELICITY
        for( int ieppV = 0; ieppV < neppV; ++ieppV )
        {
          const int ievt = ievt00 + ieppV;
          allMEs[ievt] = 0;
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
          const int ievt2 = ievt00 + ieppV + neppV;
          allMEs[ievt2] = 0;
#endif
        }
        constexpr fptype_sv* jamp2_sv = nullptr; // no need for color selection during helicity filtering
        //std::cout << "sigmaKin_getGoodHel ihel=" << ihel << ( isGoodHel[ihel] ? " true" : " false" ) << std::endl;
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL /* clang-format off */
        constexpr unsigned int channelId = 0; // disable multichannel single-diagram enhancement
        calculate_wavefunctions( ihel, allmomenta, allcouplings, allMEs, channelId, allNumerators, allDenominators, jamp2_sv, ievt00 ); //maxtry?
#else
        calculate_wavefunctions( ihel, allmomenta, allcouplings, allMEs, jamp2_sv, ievt00 ); //maxtry?
#endif /* clang-format on */
        for( int ieppV = 0; ieppV < neppV; ++ieppV )
        {
          const int ievt = ievt00 + ieppV;
          //std::cout << "sigmaKin_getGoodHel allMEs[ievt]=" << allMEs[ievt] << std::endl;
          if( allMEs[ievt] != 0 ) // NEW IMPLEMENTATION OF GETGOODHEL (#630): COMPARE EACH HELICITY CONTRIBUTION TO 0
          {
            //if ( !isGoodHel[ihel] ) std::cout << "sigmaKin_getGoodHel ihel=" << ihel << " TRUE" << std::endl;
            isGoodHel[ihel] = true;
          }
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
          const int ievt2 = ievt00 + ieppV + neppV;
          if( allMEs[ievt2] != 0 ) // NEW IMPLEMENTATION OF GETGOODHEL (#630): COMPARE EACH HELICITY CONTRIBUTION TO 0
          {
            //if ( !isGoodHel[ihel] ) std::cout << "sigmaKin_getGoodHel ihel=" << ihel << " TRUE" << std::endl;
            isGoodHel[ihel] = true;
          }
#endif
        }
      }
    }
  }
#endif

  //--------------------------------------------------------------------------

  int                                          // output: nGoodHel (the number of good helicity combinations out of ncomb)
  sigmaKin_setGoodHel( const bool* isGoodHel ) // input: isGoodHel[ncomb] - host array (CUDA and C++)
  {
    int nGoodHel = 0;
    int goodHel[ncomb] = { 0 }; // all zeros https://en.cppreference.com/w/c/language/array_initialization#Notes
    for( int ihel = 0; ihel < ncomb; ihel++ )
    {
      //std::cout << "sigmaKin_setGoodHel ihel=" << ihel << ( isGoodHel[ihel] ? " true" : " false" ) << std::endl;
      if( isGoodHel[ihel] )
      {
        goodHel[nGoodHel] = ihel;
        nGoodHel++;
      }
    }
#ifdef MGONGPUCPP_GPUIMPL
    gpuMemcpyToSymbol( dcNGoodHel, &nGoodHel, sizeof( int ) );
    gpuMemcpyToSymbol( dcGoodHel, goodHel, ncomb * sizeof( int ) );
#endif
    cNGoodHel = nGoodHel;
    for( int ihel = 0; ihel < ncomb; ihel++ ) cGoodHel[ihel] = goodHel[ihel];
    return nGoodHel;
  }

  //--------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL /* clang-format off */
  __global__ void
  normalise_output( fptype* allMEs,                    // output: allMEs[nevt], |M|^2 running_sum_over_helicities
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
                    const fptype* allNumerators,       // output: multichannel numerators[nevt], running_sum_over_helicities
                    const fptype* allDenominators,     // output: multichannel denominators[nevt], running_sum_over_helicities
                    const unsigned int* allChannelIds, // input: multichannel channelIds[nevt] (1 to #diagrams); nullptr to disable SDE enhancement (fix #899/#911)
#endif
                    const fptype globaldenom ) /* clang-format on */
  {
    const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread)
    allMEs[ievt] /= globaldenom;
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( allChannelIds != nullptr ) allMEs[ievt] *= allNumerators[ievt] / allDenominators[ievt]; // fix segfault #892 (not 'channelIds[0] != 0')
#endif
    return;
  }
#endif

  //--------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL
  __global__ void
  select_hel( int* allselhel,             // output: helicity selection[nevt]
              const fptype* allrndhel,    // input: random numbers[nevt] for helicity selection
              const fptype* allMEs_ighel, // input: allMEs_ighel[nGoodHel][nevt], |M|^2 running_sum_over_helicities
              const int nevt )            // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
  {
    const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread)
    // Event-by-event random choice of helicity #403
    //printf( "select_hel: ievt=%4d rndhel=%f\n", ievt, allrndhel[ievt] );
    for( int ighel = 0; ighel < dcNGoodHel; ighel++ )
    {
      if( allrndhel[ievt] < ( allMEs_ighel[ighel * nevt + ievt] / allMEs_ighel[( dcNGoodHel - 1 ) * nevt + ievt] ) )
      {
        const int ihelF = dcGoodHel[ighel] + 1; // NB Fortran [1,ncomb], cudacpp [0,ncomb-1]
        allselhel[ievt] = ihelF;
        //printf( "select_hel: ievt=%4d ihel=%4d\n", ievt, ihelF );
        break;
      }
    }
    return;
  }
#endif

  //--------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
  __global__ void
  select_col( int* allselcol,                    // output: color selection[nevt]
              const fptype* allrndcol,           // input: random numbers[nevt] for color selection
              const unsigned int* allChannelIds, // input: multichannel channelIds[nevt] (1 to #diagrams); nullptr to disable SDE enhancement (fix #899/#911)
              const fptype_sv* allJamp2s,        // input: jamp2[ncolor][nevt] for color choice (nullptr if disabled)
              const int nevt )                   // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
  {
    const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread)
    // SCALAR channelId for the current event (CUDA)
    unsigned int channelId = gpu_channelId( allChannelIds );
    // Event-by-event random choice of color #402
    if( channelId != 0 ) // no event-by-event choice of color if channelId == 0 (fix FPE #783)
    {
      if( channelId > mgOnGpu::nchannels )
      {
        printf( "INTERNAL ERROR! Cannot choose an event-by-event random color for channelId=%d which is greater than nchannels=%d\n", channelId, mgOnGpu::nchannels );
        assert( channelId <= mgOnGpu::nchannels ); // SANITY CHECK #919 #910
      }
      // Determine the jamp2 for this event (TEMPORARY? could do this with a dedicated memory accessor instead...)
      fptype_sv jamp2_sv[ncolor] = { 0 };
      assert( allJamp2s != nullptr ); // sanity check
      using J_ACCESS = DeviceAccessJamp2;
      for( int icolC = 0; icolC < ncolor; icolC++ )
        jamp2_sv[icolC] = J_ACCESS::kernelAccessIcolConst( allJamp2s, icolC );
      // NB (see #877): in the array channel2iconfig, the input index uses C indexing (channelId -1), the output index uses F indexing (iconfig)
      // NB (see #917): mgOnGpu::channel2iconfig returns an int (which may be -1), not an unsigned int!
      const int iconfig = mgOnGpu::channel2iconfig[channelId - 1]; // map N_diagrams to N_config <= N_diagrams configs (fix LHE color mismatch #856: see also #826, #852, #853)
      if( iconfig <= 0 )
      {
        printf( "INTERNAL ERROR! Cannot choose an event-by-event random color for channelId=%d which has no associated SDE iconfig\n", channelId );
        assert( iconfig > 0 ); // SANITY CHECK #917
      }
      else if( iconfig > (int)mgOnGpu::nconfigSDE )
      {
        printf( "INTERNAL ERROR! Cannot choose an event-by-event random color for channelId=%d (invalid SDE iconfig=%d\n > nconfig=%d)", channelId, iconfig, mgOnGpu::nconfigSDE );
        assert( iconfig <= (int)mgOnGpu::nconfigSDE ); // SANITY CHECK #917
      }
      fptype targetamp[ncolor] = { 0 };
      // NB (see #877): explicitly use 'icolC' rather than 'icol' to indicate that icolC uses C indexing in [0, N_colors-1]
      for( int icolC = 0; icolC < ncolor; icolC++ )
      {
        if( icolC == 0 )
          targetamp[icolC] = 0;
        else
          targetamp[icolC] = targetamp[icolC - 1];
        // NB (see #877): in the array icolamp, the input index uses C indexing (iconfig -1)
        if( mgOnGpu::icolamp[iconfig - 1][icolC] ) targetamp[icolC] += jamp2_sv[icolC];
      }
      //printf( "sigmaKin: ievt=%4d rndcol=%f\n", ievt, allrndcol[ievt] );
      for( int icolC = 0; icolC < ncolor; icolC++ )
      {
        if( allrndcol[ievt] < ( targetamp[icolC] / targetamp[ncolor - 1] ) )
        {
          allselcol[ievt] = icolC + 1; // NB Fortran [1,ncolor], cudacpp [0,ncolor-1]
          //printf( "sigmaKin: ievt=%d icol=%d\n", ievt, icolC+1 );
          break;
        }
      }
    }
    else
    {
      allselcol[ievt] = 0; // no color selected in Fortran range [1,ncolor] if channelId == 0 (see #931)
    }
    return;
  }
#endif
#endif

  //--------------------------------------------------------------------------
  // Evaluate |M|^2, part independent of incoming flavour

  void /* clang-format off */
  sigmaKin( const fptype* allmomenta,           // input: momenta[nevt*npar*4]
            const fptype* allcouplings,         // input: couplings[nevt*ndcoup*2]
            const fptype* allrndhel,            // input: random numbers[nevt] for helicity selection
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
            const fptype* allrndcol,            // input: random numbers[nevt] for color selection
#endif
            fptype* allMEs,                     // output: allMEs[nevt], |M|^2 final_avg_over_helicities
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
            const unsigned int* allChannelIds,  // input: multichannel channelIds[nevt] (1 to #diagrams); nullptr to disable single-diagram enhancement (fix #899/#911)
            fptype* allNumerators,              // output: multichannel numerators[nevt], running_sum_over_helicities
            fptype* allDenominators,            // output: multichannel denominators[nevt], running_sum_over_helicities
#endif
            int* allselhel,                     // output: helicity selection[nevt]
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
            int* allselcol,                     // output: color selection[nevt]
#endif
#ifdef MGONGPUCPP_GPUIMPL
            const int gpublocks,                // input: cuda gpublocks
            const int gputhreads,               // input: cuda gputhreads
            fptype* allMEs_ighel                // tmp: allMEs_ighel[nGoodHel][nevt], |M|^2 running_sum_over_helicities
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
            , fptype* allJamp2s                 // tmp: allJamp2s_icol[ncolor][nevt], |M|^2 running_sum_over_colors
#endif
#else
            , const int nevt                    // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
#endif
            ) /* clang-format on */
  {
    mgDebugInitialise();

    // SANITY CHECKS for cudacpp code generation (see issues #272 and #343 and PRs #619, #626, #360, #396 and #754)
    // These variable are not used anywhere else in the code and their scope is limited to this sanity check
    {
      // nprocesses == 2 may happen for "mirror processes" such as P0_uux_ttx within pp_tt012j (see PR #754)
      constexpr int nprocesses = 1;
      static_assert( nprocesses == 1 || nprocesses == 2, "Assume nprocesses == 1 or 2" );
      constexpr int process_id = 1; // code generation source: standalone_cudacpp
      static_assert( process_id == 1, "Assume process_id == 1" );
    }

    // Denominators: spins, colors and identical particles
    constexpr int helcolDenominators[1] = { 256 }; // assume nprocesses == 1 (#272 and #343)

#ifndef MGONGPUCPP_GPUIMPL
    //assert( (size_t)(allmomenta) % mgOnGpu::cppAlign == 0 ); // SANITY CHECK: require SIMD-friendly alignment [COMMENT OUT TO TEST MISALIGNED ACCESS]
    //assert( (size_t)(allMEs) % mgOnGpu::cppAlign == 0 ); // SANITY CHECK: require SIMD-friendly alignment [COMMENT OUT TO TEST MISALIGNED ACCESS]
    using E_ACCESS = HostAccessMatrixElements; // non-trivial access: buffer includes all events
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    using NUM_ACCESS = HostAccessNumerators;   // non-trivial access: buffer includes all events
    using DEN_ACCESS = HostAccessDenominators; // non-trivial access: buffer includes all events
    using CID_ACCESS = HostAccessChannelIds;   // non-trivial access: buffer includes all events
#endif
#endif

    // Start sigmaKin_lines
#include "GpuAbstraction.h"

    // === PART 0 - INITIALISATION (before calculate_wavefunctions) ===
    // Reset the "matrix elements" - running sums of |M|^2 over helicities for the given event
#ifdef MGONGPUCPP_GPUIMPL
    const int nevt = gpublocks * gputhreads;
    gpuMemset( allMEs, 0, nevt * sizeof( fptype ) );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    gpuMemset( allJamp2s, 0, nevt * ncolor * sizeof( fptype ) );
    gpuMemset( allNumerators, 0, nevt * sizeof( fptype ) );
    gpuMemset( allDenominators, 0, nevt * sizeof( fptype ) );
#endif
#else
    const int npagV = nevt / neppV;
    for( int ipagV = 0; ipagV < npagV; ++ipagV )
    {
      const int ievt0 = ipagV * neppV;
      fptype* MEs = E_ACCESS::ieventAccessRecord( allMEs, ievt0 );
      fptype_sv& MEs_sv = E_ACCESS::kernelAccess( MEs );
      MEs_sv = fptype_sv{ 0 };
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      fptype* numerators = NUM_ACCESS::ieventAccessRecord( allNumerators, ievt0 );
      fptype* denominators = DEN_ACCESS::ieventAccessRecord( allDenominators, ievt0 );
      fptype_sv& numerators_sv = NUM_ACCESS::kernelAccess( numerators );
      fptype_sv& denominators_sv = DEN_ACCESS::kernelAccess( denominators );
      numerators_sv = fptype_sv{ 0 };
      denominators_sv = fptype_sv{ 0 };
#endif
    }
#endif

    // === PART 1 - HELICITY LOOP: CALCULATE WAVEFUNCTIONS ===
    // (in both CUDA and C++, using precomputed good helicities)

#ifdef MGONGPUCPP_GPUIMPL // CUDA OR C++

    // *** START OF PART 1a - CUDA (one event per GPU thread) ***
    // Running sum of partial amplitudes squared for event by event color selection (#402)
    // (for the single event processed in calculate_wavefunctions)
    for( int ighel = 0; ighel < cNGoodHel; ighel++ )
    {
      const int ihel = cGoodHel[ighel];
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      gpuLaunchKernel( calculate_wavefunctions, gpublocks, gputhreads, ihel, allmomenta, allcouplings, allMEs, allChannelIds, allNumerators, allDenominators, allJamp2s, gpublocks * gputhreads );
#else
      gpuLaunchKernel( calculate_wavefunctions, gpublocks, gputhreads, ihel, allmomenta, allcouplings, allMEs, gpublocks * gputhreads );
#endif
      gpuMemcpy( &( allMEs_ighel[ighel * nevt] ), allMEs, nevt * sizeof( fptype ), gpuMemcpyDeviceToDevice );
    }
    checkGpu( gpuDeviceSynchronize() ); // do not start helicity/color selection until the loop over helicities has completed
    // Event-by-event random choice of helicity #403
    gpuLaunchKernel( select_hel, gpublocks, gputhreads, allselhel, allrndhel, allMEs_ighel, gpublocks * gputhreads );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Event-by-event random choice of color #402
    gpuLaunchKernel( select_col, gpublocks, gputhreads, allselcol, allrndcol, allChannelIds, allJamp2s, gpublocks * gputhreads );
#endif
    // *** END OF PART 1a - CUDA (one event per GPU thread) ***

#else // CUDA OR C++

    // *** START OF PART 1b - C++ (loop on event pages)
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
    // Mixed fptypes #537: float for color algebra and double elsewhere
    // Delay color algebra and ME updates (only on even pages)
    assert( npagV % 2 == 0 );     // SANITY CHECK for mixed fptypes: two neppV-pages are merged to one 2*neppV-page
    const int npagV2 = npagV / 2; // loop on two SIMD pages (neppV events) at a time
#else
    const int npagV2 = npagV;            // loop on one SIMD page (neppV events) at a time
#endif
#ifdef _OPENMP
    // OMP multithreading #575 (NB: tested only with gcc11 so far)
    // See https://www.openmp.org/specifications/
    // - default(none): no variables are shared by default
    // - shared: as the name says
    // - private: give each thread its own copy, without initialising
    // - firstprivate: give each thread its own copy, and initialise with value from outside
#define _OMPLIST0 allcouplings, allMEs, allmomenta, allrndcol, allrndhel, allselcol, allselhel, cGoodHel, cNGoodHel, npagV2
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
#define _OMPLIST1 , allDenominators, allNumerators, allChannelIds, mgOnGpu::icolamp, mgOnGpu::channel2iconfig
#else
#define _OMPLIST1
#endif
#pragma omp parallel for default( none ) shared( _OMPLIST0 _OMPLIST1 )
#undef _OMPLIST0
#undef _OMPLIST1
#endif // _OPENMP
    for( int ipagV2 = 0; ipagV2 < npagV2; ++ipagV2 )
    {
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
      const int ievt00 = ipagV2 * neppV * 2; // loop on two SIMD pages (neppV events) at a time
#else
      const int ievt00 = ipagV2 * neppV; // loop on one SIMD page (neppV events) at a time
#endif
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // SCALAR channelId for the whole SIMD neppV2 event page (C++), i.e. one or two neppV event page(s)
      // The cudacpp implementation ASSUMES (and checks! #898) that all channelIds are the same in a neppV2 SIMD event page
      // **NB! in "mixed" precision, using SIMD, calculate_wavefunctions computes MEs for TWO neppV pages with a single channelId! #924
      unsigned int channelId = 0; // disable multichannel single-diagram enhancement unless allChannelIds != nullptr
      if( allChannelIds != nullptr )
      {
        // First - and/or only - neppV page of channels (iParity=0 => ievt0 = ievt00 + 0 * neppV)
        const unsigned int* channelIds = CID_ACCESS::ieventAccessRecordConst( allChannelIds, ievt00 ); // fix bug #899/#911
        uint_sv channelIds_sv = CID_ACCESS::kernelAccessConst( channelIds );                           // fix #895 (compute this only once for all diagrams)
#ifndef MGONGPU_CPPSIMD
        // NB: channelIds_sv is a scalar in no-SIMD C++
        channelId = channelIds_sv;
#else
        // NB: channelIds_sv is a vector in SIMD C++
        channelId = channelIds_sv[0];    // element[0]
        for( int i = 1; i < neppV; ++i ) // elements[1...neppV-1]
        {
          assert( channelId == channelIds_sv[i] ); // SANITY CHECK #898: check that all events in a SIMD vector have the same channelId
        }
#endif
        assert( channelId > 0 ); // SANITY CHECK: scalar channelId must be > 0 if multichannel is enabled (allChannelIds != nullptr)
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
        // Second neppV page of channels (iParity=1 => ievt0 = ievt00 + 1 * neppV)
        const unsigned int* channelIds2 = CID_ACCESS::ieventAccessRecordConst( allChannelIds, ievt00 + neppV ); // fix bug #899/#911
        uint_v channelIds2_v = CID_ACCESS::kernelAccessConst( channelIds2 );                                    // fix #895 (compute this only once for all diagrams)
        // **NB! in "mixed" precision, using SIMD, calculate_wavefunctions computes MEs for TWO neppV pages with a single channelId! #924
        for( int i = 0; i < neppV; ++i )
        {
          assert( channelId == channelIds2_v[i] ); // SANITY CHECKS #898 #924: all events in the 2nd SIMD vector have the same channelId as that of the 1st SIMD vector
        }
#endif
      }
#endif
      // Running sum of partial amplitudes squared for event by event color selection (#402)
      // (jamp2[nParity][ncolor][neppV] for the SIMD vector - or the two SIMD vectors - of events processed in calculate_wavefunctions)
      fptype_sv jamp2_sv[nParity * ncolor] = { 0 };
      fptype_sv MEs_ighel[ncomb] = { 0 };  // sum of MEs for all good helicities up to ighel (for the first - and/or only - neppV page)
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
      fptype_sv MEs_ighel2[ncomb] = { 0 }; // sum of MEs for all good helicities up to ighel (for the second neppV page)
#endif
      for( int ighel = 0; ighel < cNGoodHel; ighel++ )
      {
        const int ihel = cGoodHel[ighel];
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
        // **NB! in "mixed" precision, using SIMD, calculate_wavefunctions computes MEs for TWO neppV pages with a single channelId! #924
        calculate_wavefunctions( ihel, allmomenta, allcouplings, allMEs, channelId, allNumerators, allDenominators, jamp2_sv, ievt00 );
#else
        calculate_wavefunctions( ihel, allmomenta, allcouplings, allMEs, jamp2_sv, ievt00 );
#endif
        MEs_ighel[ighel] = E_ACCESS::kernelAccess( E_ACCESS::ieventAccessRecord( allMEs, ievt00 ) );
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
        MEs_ighel2[ighel] = E_ACCESS::kernelAccess( E_ACCESS::ieventAccessRecord( allMEs, ievt00 + neppV ) );
#endif
      }
      // Event-by-event random choice of helicity #403
      for( int ieppV = 0; ieppV < neppV; ++ieppV )
      {
        const int ievt = ievt00 + ieppV;
        //printf( "sigmaKin: ievt=%4d rndhel=%f\n", ievt, allrndhel[ievt] );
        for( int ighel = 0; ighel < cNGoodHel; ighel++ )
        {
#if defined MGONGPU_CPPSIMD
          const bool okhel = allrndhel[ievt] < ( MEs_ighel[ighel][ieppV] / MEs_ighel[cNGoodHel - 1][ieppV] );
#else
          const bool okhel = allrndhel[ievt] < ( MEs_ighel[ighel] / MEs_ighel[cNGoodHel - 1] );
#endif
          if( okhel )
          {
            const int ihelF = cGoodHel[ighel] + 1; // NB Fortran [1,ncomb], cudacpp [0,ncomb-1]
            allselhel[ievt] = ihelF;
            //printf( "sigmaKin: ievt=%4d ihel=%4d\n", ievt, ihelF );
            break;
          }
        }
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
        const int ievt2 = ievt00 + ieppV + neppV;
        //printf( "sigmaKin: ievt=%4d rndhel=%f\n", ievt2, allrndhel[ievt2] );
        for( int ighel = 0; ighel < cNGoodHel; ighel++ )
        {
          if( allrndhel[ievt2] < ( MEs_ighel2[ighel][ieppV] / MEs_ighel2[cNGoodHel - 1][ieppV] ) )
          {
            const int ihelF = cGoodHel[ighel] + 1; // NB Fortran [1,ncomb], cudacpp [0,ncomb-1]
            allselhel[ievt2] = ihelF;
            //printf( "sigmaKin: ievt=%4d ihel=%4d\n", ievt, ihelF );
            break;
          }
        }
#endif
      }
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL // multichannel enabled (random color choice)
      // Event-by-event random choice of color #402
      if( channelId != 0 ) // no event-by-event choice of color if channelId == 0 (fix FPE #783)
      {
        if( channelId > mgOnGpu::nchannels )
        {
          printf( "INTERNAL ERROR! Cannot choose an event-by-event random color for channelId=%d which is greater than nchannels=%d\n", channelId, mgOnGpu::nchannels );
          assert( channelId <= mgOnGpu::nchannels ); // SANITY CHECK #919 #910
        }
        // NB (see #877): in the array channel2iconfig, the input index uses C indexing (channelId -1), the output index uses F indexing (iconfig)
        // NB (see #917): mgOnGpu::channel2iconfig returns an int (which may be -1), not an unsigned int!
        const int iconfig = mgOnGpu::channel2iconfig[channelId - 1]; // map N_diagrams to N_config <= N_diagrams configs (fix LHE color mismatch #856: see also #826, #852, #853)
        if( iconfig <= 0 )
        {
          printf( "INTERNAL ERROR! Cannot choose an event-by-event random color for channelId=%d which has no associated SDE iconfig\n", channelId );
          assert( iconfig > 0 ); // SANITY CHECK #917
        }
        else if( iconfig > (int)mgOnGpu::nconfigSDE )
        {
          printf( "INTERNAL ERROR! Cannot choose an event-by-event random color for channelId=%d (invalid SDE iconfig=%d\n > nconfig=%d)", channelId, iconfig, mgOnGpu::nconfigSDE );
          assert( iconfig <= (int)mgOnGpu::nconfigSDE ); // SANITY CHECK #917
        }
        fptype_sv targetamp[ncolor] = { 0 };
        // NB (see #877): explicitly use 'icolC' rather than 'icol' to indicate that icolC uses C indexing in [0, N_colors-1]
        for( int icolC = 0; icolC < ncolor; icolC++ )
        {
          if( icolC == 0 )
            targetamp[icolC] = fptype_sv{ 0 };
          else
            targetamp[icolC] = targetamp[icolC - 1];
          if( mgOnGpu::icolamp[iconfig - 1][icolC] ) targetamp[icolC] += jamp2_sv[icolC];
        }
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
        fptype_sv targetamp2[ncolor] = { 0 };
        for( int icolC = 0; icolC < ncolor; icolC++ )
        {
          if( icolC == 0 )
            targetamp2[icolC] = fptype_sv{ 0 };
          else
            targetamp2[icolC] = targetamp2[icolC - 1];
          // NB (see #877): in the array icolamp, the input index uses C indexing (iconfig -1)
          if( mgOnGpu::icolamp[iconfig - 1][icolC] ) targetamp2[icolC] += jamp2_sv[ncolor + icolC];
        }
#endif
        for( int ieppV = 0; ieppV < neppV; ++ieppV )
        {
          const int ievt = ievt00 + ieppV;
          //printf( "sigmaKin: ievt=%4d rndcol=%f\n", ievt, allrndcol[ievt] );
          for( int icolC = 0; icolC < ncolor; icolC++ )
          {
#if defined MGONGPU_CPPSIMD
            // Add volatile here to avoid SIGFPE crashes in FPTYPE=f cpp512z builds (#845)
            volatile const bool okcol = allrndcol[ievt] < ( targetamp[icolC][ieppV] / targetamp[ncolor - 1][ieppV] );
#else
            const bool okcol = allrndcol[ievt] < ( targetamp[icolC] / targetamp[ncolor - 1] );
#endif
            if( okcol )
            {
              allselcol[ievt] = icolC + 1; // NB Fortran [1,ncolor], cudacpp [0,ncolor-1]
              //printf( "sigmaKin: ievt=%d icol=%d\n", ievt, icolC+1 );
              break;
            }
          }
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
          const int ievt2 = ievt00 + ieppV + neppV;
          //printf( "sigmaKin: ievt=%4d rndcol=%f\n", ievt2, allrndcol[ievt2] );
          for( int icolC = 0; icolC < ncolor; icolC++ )
          {
            if( allrndcol[ievt2] < ( targetamp2[icolC][ieppV] / targetamp2[ncolor - 1][ieppV] ) )
            {
              allselcol[ievt2] = icolC + 1; // NB Fortran [1,ncolor], cudacpp [0,ncolor-1]
              //printf( "sigmaKin: ievt2=%d icol=%d\n", ievt2, icolC+1 );
              break;
            }
          }
#endif
        }
      }
      else
      {
        for( int ieppV = 0; ieppV < neppV; ++ieppV )
        {
          const int ievt = ievt00 + ieppV;
          allselcol[ievt] = 0; // no color selected in Fortran range [1,ncolor] if channelId == 0 (see #931)
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
          const int ievt2 = ievt00 + ieppV + neppV;
          allselcol[ievt2] = 0; // no color selected in Fortran range [1,ncolor] if channelId == 0 (see #931)
#endif
        }
      }
#endif // multichannel enabled (random color choice)
    }
    // *** END OF PART 1b - C++ (loop on event pages)

#endif // CUDA or C++

    // PART 2 - FINALISATION (after calculate_wavefunctions)
    // Get the final |M|^2 as an average over helicities/colors of the running sum of |M|^2 over helicities for the given event
    // [NB 'sum over final spins, average over initial spins', eg see
    // https://www.uzh.ch/cmsssl/physik/dam/jcr:2e24b7b1-f4d7-4160-817e-47b13dbf1d7c/Handout_4_2016-UZH.pdf]
    // (TODO OM: see how to handle the renormalization here... dedicated kernel or move it within each calculate_wavefunctions?)
#ifdef MGONGPUCPP_GPUIMPL
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    gpuLaunchKernel( normalise_output, gpublocks, gputhreads, allMEs, allNumerators, allDenominators, allChannelIds, helcolDenominators[0] );
#else
    gpuLaunchKernel( normalise_output, gpublocks, gputhreads, allMEs, helcolDenominators[0] );
#endif
#else
    for( int ipagV = 0; ipagV < npagV; ++ipagV )
    {
      const int ievt0 = ipagV * neppV;
      fptype* MEs = E_ACCESS::ieventAccessRecord( allMEs, ievt0 );
      fptype_sv& MEs_sv = E_ACCESS::kernelAccess( MEs );
      MEs_sv /= helcolDenominators[0];
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( allChannelIds != nullptr ) // fix segfault #892 (not 'channelIds[0] != 0')
      {
        fptype* numerators = NUM_ACCESS::ieventAccessRecord( allNumerators, ievt0 );
        fptype* denominators = DEN_ACCESS::ieventAccessRecord( allDenominators, ievt0 );
        fptype_sv& numerators_sv = NUM_ACCESS::kernelAccess( numerators );
        fptype_sv& denominators_sv = DEN_ACCESS::kernelAccess( denominators );
        MEs_sv *= numerators_sv / denominators_sv;
      }
#endif
      //for( int ieppV = 0; ieppV < neppV; ieppV++ )
      //{
      //  const unsigned int ievt = ipagV * neppV + ieppV;
      //  printf( "sigmaKin: ievt=%2d me=%f\n", ievt, allMEs[ievt] );
      //}
    }
#endif
    mgDebugFinalise();
  }

  //--------------------------------------------------------------------------

} // end namespace

//==========================================================================
