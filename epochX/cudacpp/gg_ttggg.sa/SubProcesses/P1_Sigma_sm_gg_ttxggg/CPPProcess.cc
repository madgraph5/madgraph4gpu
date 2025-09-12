// Copyright (C) 2010 The MadGraph5_aMC@NLO development team and contributors.
// Created by: J. Alwall (Oct 2010) for the MG5aMC CPP backend.
//==========================================================================
// Copyright (C) 2020-2025 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Modified by: S. Roiser (Feb 2020) for the MG5aMC CUDACPP plugin.
// Further modified by: S. Hageboeck, O. Mattelaer, S. Roiser, J. Teig, A. Valassi, Z. Wettersten (2020-2025) for the MG5aMC CUDACPP plugin.
//==========================================================================
// This file has been automatically generated for CUDA/C++ standalone by
// MadGraph5_aMC@NLO v. 3.6.0, 2024-09-30
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "CPPProcess.h"

#include "mgOnGpuConfig.h"

#include "HelAmps_sm.h"
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
// Process: g g > t t~ g g g WEIGHTED<=5 @1

#ifdef MGONGPUCPP_GPUIMPL
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
  constexpr int nw6 = CPPProcess::nw6;       // dimensions of each wavefunction (HELAS KEK 91-11): e.g. 6 for e+ e- -> mu+ mu- (fermions and vectors)
  constexpr int npar = CPPProcess::npar;     // #particles in total (external = initial + final): e.g. 4 for e+ e- -> mu+ mu-
  constexpr int ncomb = CPPProcess::ncomb;   // #helicity combinations: e.g. 16 for e+ e- -> mu+ mu- (2**4 = fermion spin up/down ** npar)
  constexpr int nwf = CPPProcess::nwf;       // #wavefunctions = #external (npar) + #internal: e.g. 5 for e+ e- -> mu+ mu- (1 internal is gamma or Z)
  constexpr int ncolor = CPPProcess::ncolor; // the number of leading colors

  using Parameters_sm_dependentCouplings::ndcoup;   // #couplings that vary event by event (depend on running alphas QCD)
  using Parameters_sm_independentCouplings::nicoup; // #couplings that are fixed for all events (do not depend on running alphas QCD)

  // The number of SIMD vectors of events processed by calculate_jamps
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
  constexpr int nParity = 2;
#else
  constexpr int nParity = 1;
#endif

  // Physics parameters (masses, coupling, etc...)
  // For CUDA performance, hardcoded constexpr's would be better: fewer registers and a tiny throughput increase
  // However, physics parameters are user-defined through card files: use CUDA constant memory instead (issue #39)
  // [NB if hardcoded parameters are used, it's better to define them here to avoid silent shadowing (issue #263)]
  constexpr int nIPD = 2; // SM independent parameters used in this CPPProcess.cc (FIXME? rename as sm_IndepParam?)
  // Note: in the Python code generator, nIPD == nparam, while nIPC <= nicoup, because (see #823)
  // nIPC may vary from one P*/CPPProcess.cc to another, while nicoup is defined in src/Param.h and is common to all P*
  constexpr int nIPC = 0; // SM independent couplings used in this CPPProcess.cc (FIXME? rename as sm_IndepCoupl?)
  static_assert( nIPC <= nicoup );
  static_assert( nIPD >= 0 ); // Hack to avoid build warnings when nIPD==0 is unused
  static_assert( nIPC >= 0 ); // Hack to avoid build warnings when nIPC==0 is unused
#ifdef MGONGPU_HARDCODE_PARAM
  __device__ const fptype cIPD[nIPD] = { (fptype)Parameters_sm::mdl_MT, (fptype)Parameters_sm::mdl_WT };
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
  __device__ const double* bsmIndepParam = Parameters_sm::mdl_bsmIndepParam;
#else
#ifdef MGONGPUCPP_GPUIMPL
  __device__ __constant__ double bsmIndepParam[Parameters_sm::nBsmIndepParam];
#else
  static double bsmIndepParam[Parameters_sm::nBsmIndepParam];
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
  class DeviceAccessJamp
  {
  public:
    static __device__ inline cxtype_ref
    kernelAccessIcol( fptype* buffer, const int icol )
    {
      const int nevt = gridDim.x * blockDim.x;
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x;
      return cxtype_ref( buffer[icol * 2 * nevt + ievt], buffer[icol * 2 * nevt + nevt + ievt] );
    }
    static __device__ inline const cxtype
    kernelAccessIcolConst( const fptype* buffer, const int icol )
    {
      const int nevt = gridDim.x * blockDim.x;
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x;
      return cxtype( buffer[icol * 2 * nevt + ievt], buffer[icol * 2 * nevt + nevt + ievt] );
    }
  };
#else
  class HostAccessJamp
  {
  public:
    static inline cxtype_sv&
    kernelAccessIcol( cxtype_sv* buffer, const int icol )
    {
      return buffer[icol];
    }
    static inline cxtype_sv&
    kernelAccessIcol( fptype* buffer, const int icol )
    {
      return reinterpret_cast<cxtype_sv*>( buffer )[icol];
    }
  };
#endif

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

#include "diagrams.h"

  //--------------------------------------------------------------------------

  // Evaluate QCD partial amplitudes jamps for this given helicity from Feynman diagrams
  // Also compute running sums over helicities adding jamp2, numerator, denominator
  // (NB: this function no longer handles matrix elements as the color sum has now been moved to a separate function/kernel)
  // In CUDA, this function processes a single event
  // ** NB1: NEW Nov2024! In CUDA this is now a kernel function (it used to be a device function)
  // ** NB2: NEW Nov2024! in CUDA this now takes a channelId array as input (it used to take a scalar channelId as input)
  // In C++, this function processes a single event "page" or SIMD vector (or for two in "mixed" precision mode, nParity=2)
  // *** NB: in C++, calculate_jamps accepts a SCALAR channelId because it is GUARANTEED that all events in a SIMD vector have the same channelId #898
#ifdef MGONGPUCPP_GPUIMPL /* clang-format off */
  INLINE void
  calculate_jamps( int ihel,
                   const fptype* allmomenta,          // input: momenta[nevt*npar*4]
                   const fptype* allcouplings,        // input: couplings[nevt*ndcoup*2]
                   fptype* allJamps,                  // output: jamp[ncolor*2*nevt] for this helicity
                   fptype* allWfs,                    // output: wf[nwf*nw6*2*nevt]
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
                   const unsigned int* allChannelIds, // input: multichannel channelIds[nevt] (1 to #diagrams); nullptr to disable SDE (#899/#911)
                   fptype* allNumerators,             // input/output: multichannel numerators[nevt], add helicity ihel
                   fptype* allDenominators,           // input/output: multichannel denominators[nevt], add helicity ihel
#endif
                   gpuStream_t gpustream,             // input: cuda stream for this helicity
                   const int gpublocks,               // input: cuda gpublocks
                   const int gputhreads )             // input: cuda gputhreads
#else
  INLINE void
  calculate_jamps( int ihel,
                   const fptype* allmomenta,          // input: momenta[nevt*npar*4]
                   const fptype* allcouplings,        // input: couplings[nevt*ndcoup*2]
                   cxtype_sv* jamp_sv,                // output: jamp_sv[ncolor] (f/d) or [2*ncolor] (m) for SIMD event page(s) ievt00 and helicity ihel
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
                   const unsigned int channelId,      // input: SCALAR channelId (1 to #diagrams, 0 to disable SDE) for SIMD event page(s) ievt00
                   fptype* allNumerators,             // input/output: multichannel numerators[nevt], add helicity ihel
                   fptype* allDenominators,           // input/output: multichannel denominators[nevt], add helicity ihel
#endif
                   const int ievt00 )                 // input: first event number in current C++ event page (for CUDA, ievt depends on threadid)
#endif
  //ALWAYS_INLINE // attributes are not permitted in a function definition
  {
#ifdef MGONGPUCPP_GPUIMPL
    using CD_ACCESS = DeviceAccessCouplings;      // non-trivial access (dependent couplings): buffer includes all events
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    using NUM_ACCESS = DeviceAccessNumerators;    // non-trivial access: buffer includes all events
    using DEN_ACCESS = DeviceAccessDenominators;  // non-trivial access: buffer includes all events
#endif
#else
    using M_ACCESS = HostAccessMomenta;         // non-trivial access: buffer includes all events
    using CD_ACCESS = HostAccessCouplings;      // non-trivial access (dependent couplings): buffer includes all events
    using CI_ACCESS = HostAccessCouplingsFixed; // TRIVIAL access (independent couplings): buffer for one event
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    using NUM_ACCESS = HostAccessNumerators;    // non-trivial access: buffer includes all events
    using DEN_ACCESS = HostAccessDenominators;  // non-trivial access: buffer includes all events
#endif
#endif /* clang-format on */

    // ----------------------------
    // --- WAVEFUNCTION BUFFERS ---
    // ----------------------------
#ifndef MGONGPUCPP_GPUIMPL
    // Local TEMPORARY variables for a subset of Feynman diagrams in the given CUDA event (ievt) or C++ event page (ipagV)
    // [NB these variables are reused several times (and re-initialised each time) within the same event or event page]
    // ** NB: wavefunctions only need TRIVIAL ACCESS in C++ code
    cxtype_sv w_sv[nwf][nw6]; // particle wavefunctions within Feynman diagrams (nw6 is often 6, the dimension of spin 1/2 or spin 1 particles)
    fptype* wfs = reinterpret_cast<fptype*>( w_sv );
#else
    // Global-memory variables for a subset of Feynman diagrams in the given CUDA event (ievt)
    // ** NB: wavefunctions need non-trivial access in CUDA code because of kernel splitting
    fptype* wfs = allWfs;
#endif

    // === Calculate wavefunctions and amplitudes for all diagrams in all processes         ===
    // === (for one event in CUDA, for one - or two in mixed mode - SIMD event pages in C++ ===

    // *****************************
    // *** START LOOP ON IPARITY ***
    // *****************************
    for( int iParity = 0; iParity < nParity; ++iParity )
    {
#ifndef MGONGPUCPP_GPUIMPL
      const int ievt0 = ievt00 + iParity * neppV;
#endif

      // -----------------
      // --- COUPLINGS ---
      // -----------------
#ifdef MGONGPUCPP_GPUIMPL
      // CUDA diagram kernels take input/output buffers with couplings "fptype* couplings" for all events
      const fptype* couplings = allcouplings;
#else
      // C++ diagram kernels take input/output buffers with couplings "fptype** COUPs" for a single event or SIMD vector
      constexpr size_t nxcoup = ndcoup + nIPC; // both dependent and independent couplings (FIX #823: nIPC instead of nicoup)
      const fptype* allCOUPs[nxcoup];
      const fptype* COUPs[nxcoup];
      // Dependent couplings, vary event-by-event
      for( size_t idcoup = 0; idcoup < ndcoup; idcoup++ )
        allCOUPs[idcoup] = CD_ACCESS::idcoupAccessBufferConst( allcouplings, idcoup );
      // Independent couplings, fixed for all events
      for( size_t iicoup = 0; iicoup < nIPC; iicoup++ ) // (FIX #823: nIPC instead of nicoup)
        allCOUPs[ndcoup + iicoup] = CI_ACCESS::iicoupAccessBufferConst( cIPC, iicoup );
      // Dependent couplings, vary event-by-event
      for( size_t idcoup = 0; idcoup < ndcoup; idcoup++ )
        COUPs[idcoup] = CD_ACCESS::ieventAccessRecordConst( allCOUPs[idcoup], ievt0 );
      // Independent couplings, fixed for all events
      for( size_t iicoup = 0; iicoup < nIPC; iicoup++ ) // (FIX #823: nIPC instead of nicoup)
        COUPs[ndcoup + iicoup] = allCOUPs[ndcoup + iicoup];
#endif

      // ---------------
      // --- MOMENTA ---
      // ---------------
#ifdef MGONGPUCPP_GPUIMPL
      // CUDA diagram kernels take input/output buffers with momenta for all events
      const fptype* momenta = allmomenta;
#else
      // C++ diagram kernels take input/output buffers with momenta for a single event or SIMD vector
      const fptype* momenta = M_ACCESS::ieventAccessRecordConst( allmomenta, ievt0 );
#endif
      // -------------
      // --- JAMPS ---
      // -------------
      // (Note: no need to 'reset color flows' i.e. zero allJamps, this is done in sigmaKin and sigmaKin_getGoodHel)
#ifdef MGONGPUCPP_GPUIMPL
      // In CUDA, write jamps to the output global-memory allJamps [for all events] passed as argument
      // (write directly to J_ACCESS::kernelAccessIcol( allJamps, icol ) instead of writing to jamp_sv[icol])
      fptype* jamps = allJamps;
#else
      // In C++, write jamps to the output array [for one specific event or SIMD vector] passed as argument
      // (write directly to J_ACCESS::kernelAccessIcol( allJamps, icol ) instead of writing to jamp_sv[icol])
      fptype* jamps = reinterpret_cast<fptype*>( iParity == 0 ? jamp_sv : &( jamp_sv[ncolor] ) );
#endif

      // ------------------
      // --- CHANNELIDS ---
      // ------------------
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
#ifdef MGONGPUCPP_GPUIMPL
      // CUDA diagram kernels take input/output buffers with channelIDs for all events
      const unsigned int* channelIds = allChannelIds;
#else
      // C++ diagram kernels take input/output buffers with a single SCALAR channelID for all events in a given SIMD vector
      const unsigned int* channelIds = &channelId;
#endif
#else
      // A uniform interface for diagramXXX including channelIDs, numerators and denominators is used also #ifndef MGONGPU_SUPPORTS_MULTICHANNEL
      // In that case, however, the boilerplate code asserts that all three pointers all nullptr as a sanity check
      const unsigned int* channelIds = nullptr;
#endif

      // -------------------------------
      // --- NUMERATORS/DENOMINATORS ---
      // -------------------------------
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
#ifdef MGONGPUCPP_GPUIMPL
      // CUDA diagram kernels take input/output buffers with numerators/denominators for all events
      fptype* numerators = allNumerators;
      fptype* denominators = allDenominators;
#else
      // C++ diagram kernels take input/output buffers with numerators/denominators for a single event or SIMD vector
      fptype* numerators = NUM_ACCESS::ieventAccessRecord( allNumerators, ievt0 );
      fptype* denominators = DEN_ACCESS::ieventAccessRecord( allDenominators, ievt0 );
#endif
#else
      // A uniform interface for diagramXXX including channelIDs, numerators and denominators is used also #ifndef MGONGPU_SUPPORTS_MULTICHANNEL
      // In that case, however, the boilerplate code asserts that all three pointers all nullptr as a sanity check
      fptype* numerators = nullptr;
      fptype* denominators = nullptr;
#endif

      // ------------------------
      // --- FEYNMAN DIAGRAMS ---
      // ------------------------

      // *** DIAGRAMS 1 TO 1240 ***
#ifdef MGONGPUCPP_GPUIMPL
      gpuLaunchKernelStream( diagram1, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators, momenta, ihel );
      gpuLaunchKernelStream( diagram2, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram3, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram4, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram5, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram6, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram7, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram8, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram9, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram10, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram11, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram12, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram13, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram14, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram15, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram16, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram17, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram18, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram19, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram20, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram21, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram22, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram23, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram24, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram25, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram26, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram27, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram28, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram29, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram30, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram31, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram32, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram33, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram34, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram35, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram36, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram37, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram38, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram39, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram40, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram41, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram42, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram43, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram44, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram45, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram46, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram47, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram48, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram49, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram50, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram51, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram52, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram53, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram54, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram55, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram56, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram57, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram58, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram59, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram60, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram61, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram62, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram63, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram64, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram65, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram66, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram67, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram68, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram69, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram70, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram71, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram72, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram73, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram74, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram75, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram76, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram77, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram78, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram79, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram80, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram81, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram82, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram83, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram84, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram85, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram86, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram87, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram88, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram89, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram90, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram91, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram92, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram93, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram94, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram95, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram96, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram97, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram98, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram99, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram100, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram101, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram102, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram103, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram104, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram105, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram106, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram107, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram108, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram109, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram110, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram111, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram112, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram113, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram114, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram115, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram116, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram117, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram118, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram119, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram120, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram121, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram122, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram123, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram124, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram125, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram126, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram127, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram128, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram129, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram130, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram131, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram132, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram133, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram134, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram135, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram136, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram137, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram138, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram139, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram140, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram141, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram142, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram143, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram144, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram145, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram146, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram147, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram148, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram149, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram150, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram151, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram152, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram153, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram154, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram155, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram156, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram157, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram158, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram159, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram160, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram161, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram162, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram163, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram164, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram165, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram166, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram167, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram168, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram169, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram170, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram171, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram172, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram173, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram174, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram175, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram176, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram177, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram178, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram179, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram180, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram181, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram182, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram183, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram184, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram185, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram186, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram187, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram188, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram189, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram190, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram191, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram192, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram193, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram194, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram195, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram196, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram197, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram198, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram199, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram200, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram201, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram202, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram203, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram204, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram205, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram206, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram207, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram208, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram209, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram210, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram211, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram212, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram213, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram214, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram215, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram216, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram217, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram218, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram219, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram220, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram221, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram222, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram223, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram224, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram225, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram226, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram227, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram228, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram229, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram230, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram231, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram232, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram233, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram234, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram235, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram236, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram237, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram238, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram239, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram240, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram241, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram242, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram243, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram244, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram245, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram246, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram247, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram248, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram249, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram250, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram251, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram252, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram253, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram254, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram255, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram256, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram257, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram258, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram259, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram260, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram261, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram262, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram263, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram264, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram265, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram266, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram267, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram268, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram269, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram270, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram271, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram272, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram273, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram274, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram275, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram276, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram277, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram278, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram279, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram280, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram281, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram282, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram283, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram284, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram285, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram286, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram287, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram288, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram289, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram290, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram291, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram292, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram293, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram294, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram295, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram296, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram297, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram298, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram299, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram300, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram301, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram302, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram303, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram304, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram305, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram306, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram307, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram308, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram309, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram310, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram311, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram312, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram313, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram314, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram315, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram316, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram317, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram318, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram319, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram320, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram321, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram322, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram323, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram324, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram325, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram326, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram327, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram328, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram329, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram330, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram331, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram332, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram333, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram334, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram335, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram336, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram337, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram338, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram339, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram340, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram341, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram342, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram343, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram344, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram345, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram346, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram347, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram348, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram349, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram350, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram351, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram352, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram353, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram354, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram355, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram356, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram357, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram358, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram359, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram360, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram361, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram362, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram363, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram364, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram365, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram366, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram367, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram368, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram369, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram370, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram371, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram372, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram373, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram374, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram375, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram376, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram377, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram378, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram379, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram380, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram381, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram382, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram383, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram384, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram385, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram386, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram387, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram388, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram389, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram390, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram391, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram392, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram393, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram394, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram395, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram396, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram397, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram398, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram399, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram400, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram401, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram402, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram403, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram404, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram405, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram406, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram407, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram408, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram409, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram410, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram411, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram412, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram413, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram414, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram415, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram416, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram417, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram418, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram419, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram420, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram421, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram422, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram423, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram424, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram425, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram426, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram427, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram428, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram429, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram430, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram431, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram432, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram433, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram434, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram435, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram436, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram437, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram438, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram439, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram440, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram441, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram442, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram443, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram444, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram445, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram446, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram447, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram448, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram449, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram450, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram451, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram452, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram453, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram454, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram455, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram456, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram457, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram458, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram459, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram460, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram461, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram462, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram463, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram464, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram465, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram466, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram467, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram468, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram469, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram470, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram471, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram472, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram473, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram474, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram475, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram476, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram477, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram478, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram479, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram480, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram481, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram482, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram483, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram484, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram485, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram486, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram487, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram488, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram489, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram490, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram491, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram492, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram493, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram494, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram495, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram496, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram497, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram498, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram499, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram500, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram501, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram502, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram503, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram504, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram505, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram506, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram507, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram508, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram509, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram510, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram511, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram512, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram513, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram514, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram515, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram516, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram517, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram518, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram519, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram520, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram521, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram522, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram523, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram524, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram525, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram526, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram527, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram528, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram529, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram530, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram531, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram532, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram533, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram534, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram535, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram536, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram537, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram538, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram539, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram540, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram541, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram542, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram543, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram544, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram545, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram546, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram547, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram548, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram549, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram550, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram551, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram552, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram553, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram554, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram555, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram556, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram557, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram558, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram559, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram560, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram561, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram562, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram563, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram564, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram565, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram566, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram567, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram568, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram569, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram570, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram571, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram572, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram573, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram574, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram575, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram576, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram577, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram578, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram579, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram580, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram581, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram582, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram583, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram584, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram585, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram586, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram587, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram588, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram589, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram590, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram591, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram592, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram593, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram594, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram595, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram596, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram597, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram598, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram599, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram600, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram601, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram602, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram603, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram604, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram605, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram606, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram607, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram608, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram609, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram610, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram611, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram612, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram613, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram614, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram615, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram616, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram617, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram618, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram619, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram620, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram621, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram622, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram623, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram624, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram625, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram626, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram627, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram628, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram629, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram630, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram631, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram632, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram633, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram634, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram635, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram636, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram637, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram638, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram639, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram640, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram641, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram642, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram643, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram644, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram645, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram646, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram647, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram648, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram649, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram650, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram651, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram652, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram653, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram654, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram655, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram656, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram657, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram658, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram659, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram660, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram661, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram662, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram663, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram664, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram665, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram666, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram667, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram668, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram669, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram670, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram671, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram672, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram673, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram674, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram675, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram676, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram677, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram678, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram679, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram680, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram681, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram682, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram683, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram684, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram685, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram686, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram687, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram688, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram689, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram690, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram691, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram692, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram693, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram694, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram695, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram696, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram697, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram698, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram699, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram700, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram701, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram702, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram703, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram704, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram705, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram706, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram707, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram708, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram709, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram710, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram711, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram712, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram713, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram714, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram715, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram716, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram717, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram718, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram719, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram720, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram721, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram722, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram723, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram724, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram725, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram726, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram727, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram728, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram729, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram730, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram731, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram732, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram733, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram734, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram735, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram736, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram737, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram738, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram739, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram740, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram741, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram742, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram743, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram744, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram745, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram746, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram747, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram748, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram749, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram750, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram751, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram752, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram753, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram754, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram755, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram756, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram757, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram758, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram759, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram760, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram761, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram762, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram763, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram764, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram765, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram766, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram767, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram768, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram769, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram770, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram771, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram772, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram773, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram774, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram775, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram776, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram777, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram778, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram779, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram780, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram781, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram782, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram783, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram784, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram785, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram786, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram787, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram788, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram789, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram790, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram791, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram792, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram793, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram794, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram795, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram796, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram797, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram798, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram799, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram800, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram801, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram802, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram803, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram804, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram805, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram806, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram807, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram808, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram809, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram810, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram811, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram812, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram813, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram814, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram815, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram816, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram817, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram818, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram819, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram820, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram821, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram822, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram823, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram824, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram825, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram826, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram827, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram828, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram829, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram830, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram831, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram832, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram833, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram834, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram835, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram836, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram837, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram838, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram839, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram840, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram841, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram842, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram843, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram844, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram845, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram846, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram847, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram848, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram849, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram850, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram851, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram852, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram853, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram854, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram855, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram856, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram857, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram858, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram859, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram860, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram861, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram862, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram863, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram864, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram865, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram866, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram867, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram868, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram869, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram870, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram871, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram872, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram873, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram874, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram875, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram876, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram877, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram878, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram879, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram880, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram881, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram882, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram883, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram884, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram885, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram886, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram887, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram888, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram889, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram890, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram891, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram892, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram893, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram894, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram895, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram896, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram897, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram898, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram899, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram900, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram901, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram902, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram903, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram904, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram905, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram906, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram907, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram908, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram909, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram910, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram911, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram912, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram913, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram914, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram915, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram916, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram917, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram918, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram919, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram920, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram921, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram922, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram923, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram924, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram925, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram926, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram927, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram928, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram929, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram930, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram931, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram932, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram933, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram934, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram935, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram936, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram937, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram938, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram939, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram940, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram941, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram942, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram943, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram944, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram945, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram946, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram947, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram948, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram949, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram950, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram951, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram952, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram953, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram954, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram955, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram956, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram957, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram958, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram959, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram960, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram961, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram962, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram963, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram964, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram965, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram966, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram967, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram968, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram969, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram970, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram971, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram972, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram973, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram974, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram975, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram976, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram977, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram978, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram979, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram980, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram981, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram982, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram983, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram984, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram985, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram986, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram987, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram988, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram989, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram990, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram991, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram992, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram993, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram994, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram995, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram996, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram997, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram998, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram999, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1000, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1001, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1002, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1003, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1004, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1005, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1006, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1007, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1008, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1009, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1010, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1011, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1012, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1013, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1014, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1015, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1016, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1017, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1018, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1019, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1020, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1021, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1022, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1023, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1024, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1025, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1026, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1027, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1028, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1029, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1030, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1031, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1032, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1033, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1034, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1035, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1036, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1037, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1038, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1039, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1040, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1041, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1042, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1043, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1044, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1045, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1046, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1047, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1048, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1049, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1050, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1051, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1052, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1053, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1054, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1055, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1056, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1057, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1058, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1059, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1060, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1061, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1062, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1063, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1064, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1065, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1066, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1067, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1068, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1069, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1070, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1071, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1072, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1073, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1074, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1075, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1076, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1077, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1078, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1079, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1080, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1081, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1082, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1083, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1084, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1085, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1086, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1087, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1088, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1089, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1090, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1091, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1092, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1093, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1094, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1095, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1096, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1097, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1098, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1099, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1100, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1101, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1102, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1103, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1104, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1105, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1106, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1107, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1108, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1109, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1110, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1111, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1112, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1113, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1114, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1115, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1116, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1117, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1118, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1119, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1120, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1121, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1122, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1123, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1124, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1125, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1126, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1127, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1128, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1129, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1130, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1131, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1132, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1133, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1134, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1135, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1136, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1137, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1138, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1139, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1140, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1141, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1142, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1143, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1144, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1145, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1146, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1147, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1148, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1149, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1150, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1151, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1152, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1153, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1154, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1155, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1156, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1157, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1158, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1159, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1160, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1161, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1162, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1163, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1164, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1165, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1166, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1167, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1168, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1169, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1170, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1171, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1172, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1173, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1174, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1175, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1176, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1177, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1178, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1179, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1180, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1181, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1182, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1183, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1184, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1185, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1186, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1187, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1188, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1189, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1190, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1191, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1192, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1193, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1194, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1195, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1196, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1197, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1198, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1199, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1200, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1201, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1202, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1203, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1204, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1205, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1206, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1207, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1208, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1209, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1210, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1211, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1212, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1213, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1214, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1215, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1216, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1217, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1218, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1219, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1220, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1221, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1222, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1223, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1224, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1225, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1226, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1227, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1228, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1229, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1230, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1231, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1232, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1233, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1234, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1235, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1236, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1237, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1238, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1239, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuLaunchKernelStream( diagram1240, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
#else
      diagram1( wfs, jamps, channelIds, COUPs, numerators, denominators, momenta, ihel );
      diagram2( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram3( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram4( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram5( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram6( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram7( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram8( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram9( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram10( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram11( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram12( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram13( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram14( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram15( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram16( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram17( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram18( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram19( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram20( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram21( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram22( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram23( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram24( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram25( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram26( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram27( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram28( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram29( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram30( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram31( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram32( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram33( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram34( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram35( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram36( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram37( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram38( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram39( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram40( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram41( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram42( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram43( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram44( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram45( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram46( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram47( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram48( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram49( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram50( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram51( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram52( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram53( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram54( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram55( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram56( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram57( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram58( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram59( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram60( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram61( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram62( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram63( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram64( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram65( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram66( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram67( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram68( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram69( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram70( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram71( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram72( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram73( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram74( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram75( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram76( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram77( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram78( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram79( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram80( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram81( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram82( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram83( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram84( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram85( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram86( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram87( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram88( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram89( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram90( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram91( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram92( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram93( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram94( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram95( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram96( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram97( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram98( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram99( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram100( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram101( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram102( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram103( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram104( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram105( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram106( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram107( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram108( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram109( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram110( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram111( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram112( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram113( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram114( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram115( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram116( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram117( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram118( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram119( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram120( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram121( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram122( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram123( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram124( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram125( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram126( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram127( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram128( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram129( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram130( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram131( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram132( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram133( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram134( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram135( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram136( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram137( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram138( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram139( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram140( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram141( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram142( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram143( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram144( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram145( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram146( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram147( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram148( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram149( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram150( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram151( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram152( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram153( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram154( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram155( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram156( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram157( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram158( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram159( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram160( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram161( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram162( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram163( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram164( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram165( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram166( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram167( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram168( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram169( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram170( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram171( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram172( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram173( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram174( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram175( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram176( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram177( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram178( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram179( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram180( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram181( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram182( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram183( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram184( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram185( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram186( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram187( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram188( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram189( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram190( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram191( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram192( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram193( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram194( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram195( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram196( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram197( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram198( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram199( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram200( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram201( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram202( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram203( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram204( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram205( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram206( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram207( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram208( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram209( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram210( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram211( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram212( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram213( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram214( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram215( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram216( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram217( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram218( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram219( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram220( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram221( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram222( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram223( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram224( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram225( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram226( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram227( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram228( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram229( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram230( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram231( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram232( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram233( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram234( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram235( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram236( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram237( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram238( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram239( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram240( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram241( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram242( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram243( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram244( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram245( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram246( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram247( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram248( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram249( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram250( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram251( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram252( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram253( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram254( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram255( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram256( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram257( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram258( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram259( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram260( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram261( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram262( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram263( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram264( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram265( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram266( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram267( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram268( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram269( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram270( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram271( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram272( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram273( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram274( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram275( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram276( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram277( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram278( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram279( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram280( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram281( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram282( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram283( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram284( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram285( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram286( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram287( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram288( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram289( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram290( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram291( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram292( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram293( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram294( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram295( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram296( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram297( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram298( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram299( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram300( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram301( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram302( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram303( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram304( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram305( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram306( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram307( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram308( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram309( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram310( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram311( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram312( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram313( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram314( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram315( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram316( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram317( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram318( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram319( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram320( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram321( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram322( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram323( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram324( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram325( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram326( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram327( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram328( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram329( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram330( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram331( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram332( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram333( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram334( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram335( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram336( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram337( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram338( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram339( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram340( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram341( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram342( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram343( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram344( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram345( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram346( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram347( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram348( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram349( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram350( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram351( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram352( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram353( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram354( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram355( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram356( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram357( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram358( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram359( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram360( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram361( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram362( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram363( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram364( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram365( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram366( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram367( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram368( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram369( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram370( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram371( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram372( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram373( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram374( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram375( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram376( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram377( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram378( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram379( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram380( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram381( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram382( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram383( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram384( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram385( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram386( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram387( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram388( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram389( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram390( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram391( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram392( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram393( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram394( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram395( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram396( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram397( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram398( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram399( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram400( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram401( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram402( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram403( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram404( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram405( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram406( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram407( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram408( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram409( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram410( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram411( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram412( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram413( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram414( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram415( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram416( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram417( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram418( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram419( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram420( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram421( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram422( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram423( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram424( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram425( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram426( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram427( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram428( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram429( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram430( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram431( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram432( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram433( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram434( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram435( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram436( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram437( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram438( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram439( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram440( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram441( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram442( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram443( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram444( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram445( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram446( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram447( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram448( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram449( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram450( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram451( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram452( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram453( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram454( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram455( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram456( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram457( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram458( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram459( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram460( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram461( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram462( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram463( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram464( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram465( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram466( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram467( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram468( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram469( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram470( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram471( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram472( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram473( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram474( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram475( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram476( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram477( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram478( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram479( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram480( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram481( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram482( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram483( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram484( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram485( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram486( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram487( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram488( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram489( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram490( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram491( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram492( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram493( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram494( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram495( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram496( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram497( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram498( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram499( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram500( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram501( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram502( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram503( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram504( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram505( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram506( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram507( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram508( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram509( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram510( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram511( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram512( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram513( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram514( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram515( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram516( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram517( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram518( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram519( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram520( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram521( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram522( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram523( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram524( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram525( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram526( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram527( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram528( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram529( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram530( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram531( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram532( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram533( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram534( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram535( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram536( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram537( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram538( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram539( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram540( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram541( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram542( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram543( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram544( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram545( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram546( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram547( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram548( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram549( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram550( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram551( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram552( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram553( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram554( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram555( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram556( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram557( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram558( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram559( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram560( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram561( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram562( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram563( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram564( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram565( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram566( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram567( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram568( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram569( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram570( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram571( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram572( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram573( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram574( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram575( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram576( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram577( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram578( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram579( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram580( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram581( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram582( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram583( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram584( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram585( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram586( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram587( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram588( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram589( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram590( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram591( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram592( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram593( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram594( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram595( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram596( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram597( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram598( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram599( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram600( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram601( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram602( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram603( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram604( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram605( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram606( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram607( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram608( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram609( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram610( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram611( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram612( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram613( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram614( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram615( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram616( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram617( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram618( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram619( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram620( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram621( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram622( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram623( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram624( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram625( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram626( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram627( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram628( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram629( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram630( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram631( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram632( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram633( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram634( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram635( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram636( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram637( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram638( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram639( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram640( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram641( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram642( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram643( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram644( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram645( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram646( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram647( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram648( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram649( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram650( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram651( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram652( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram653( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram654( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram655( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram656( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram657( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram658( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram659( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram660( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram661( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram662( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram663( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram664( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram665( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram666( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram667( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram668( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram669( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram670( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram671( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram672( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram673( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram674( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram675( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram676( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram677( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram678( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram679( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram680( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram681( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram682( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram683( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram684( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram685( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram686( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram687( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram688( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram689( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram690( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram691( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram692( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram693( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram694( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram695( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram696( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram697( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram698( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram699( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram700( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram701( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram702( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram703( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram704( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram705( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram706( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram707( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram708( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram709( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram710( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram711( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram712( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram713( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram714( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram715( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram716( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram717( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram718( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram719( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram720( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram721( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram722( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram723( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram724( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram725( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram726( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram727( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram728( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram729( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram730( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram731( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram732( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram733( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram734( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram735( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram736( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram737( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram738( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram739( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram740( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram741( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram742( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram743( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram744( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram745( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram746( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram747( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram748( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram749( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram750( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram751( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram752( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram753( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram754( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram755( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram756( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram757( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram758( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram759( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram760( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram761( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram762( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram763( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram764( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram765( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram766( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram767( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram768( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram769( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram770( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram771( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram772( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram773( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram774( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram775( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram776( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram777( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram778( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram779( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram780( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram781( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram782( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram783( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram784( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram785( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram786( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram787( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram788( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram789( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram790( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram791( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram792( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram793( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram794( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram795( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram796( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram797( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram798( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram799( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram800( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram801( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram802( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram803( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram804( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram805( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram806( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram807( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram808( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram809( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram810( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram811( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram812( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram813( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram814( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram815( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram816( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram817( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram818( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram819( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram820( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram821( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram822( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram823( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram824( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram825( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram826( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram827( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram828( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram829( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram830( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram831( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram832( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram833( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram834( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram835( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram836( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram837( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram838( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram839( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram840( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram841( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram842( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram843( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram844( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram845( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram846( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram847( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram848( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram849( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram850( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram851( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram852( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram853( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram854( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram855( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram856( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram857( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram858( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram859( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram860( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram861( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram862( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram863( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram864( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram865( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram866( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram867( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram868( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram869( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram870( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram871( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram872( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram873( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram874( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram875( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram876( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram877( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram878( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram879( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram880( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram881( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram882( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram883( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram884( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram885( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram886( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram887( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram888( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram889( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram890( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram891( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram892( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram893( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram894( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram895( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram896( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram897( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram898( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram899( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram900( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram901( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram902( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram903( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram904( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram905( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram906( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram907( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram908( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram909( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram910( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram911( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram912( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram913( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram914( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram915( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram916( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram917( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram918( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram919( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram920( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram921( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram922( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram923( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram924( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram925( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram926( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram927( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram928( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram929( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram930( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram931( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram932( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram933( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram934( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram935( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram936( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram937( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram938( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram939( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram940( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram941( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram942( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram943( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram944( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram945( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram946( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram947( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram948( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram949( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram950( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram951( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram952( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram953( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram954( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram955( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram956( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram957( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram958( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram959( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram960( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram961( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram962( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram963( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram964( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram965( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram966( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram967( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram968( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram969( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram970( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram971( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram972( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram973( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram974( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram975( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram976( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram977( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram978( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram979( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram980( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram981( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram982( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram983( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram984( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram985( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram986( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram987( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram988( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram989( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram990( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram991( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram992( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram993( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram994( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram995( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram996( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram997( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram998( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram999( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1000( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1001( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1002( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1003( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1004( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1005( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1006( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1007( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1008( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1009( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1010( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1011( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1012( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1013( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1014( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1015( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1016( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1017( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1018( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1019( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1020( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1021( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1022( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1023( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1024( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1025( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1026( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1027( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1028( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1029( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1030( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1031( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1032( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1033( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1034( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1035( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1036( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1037( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1038( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1039( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1040( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1041( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1042( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1043( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1044( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1045( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1046( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1047( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1048( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1049( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1050( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1051( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1052( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1053( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1054( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1055( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1056( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1057( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1058( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1059( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1060( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1061( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1062( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1063( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1064( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1065( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1066( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1067( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1068( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1069( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1070( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1071( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1072( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1073( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1074( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1075( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1076( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1077( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1078( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1079( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1080( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1081( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1082( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1083( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1084( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1085( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1086( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1087( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1088( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1089( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1090( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1091( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1092( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1093( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1094( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1095( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1096( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1097( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1098( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1099( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1100( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1101( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1102( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1103( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1104( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1105( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1106( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1107( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1108( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1109( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1110( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1111( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1112( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1113( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1114( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1115( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1116( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1117( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1118( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1119( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1120( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1121( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1122( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1123( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1124( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1125( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1126( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1127( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1128( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1129( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1130( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1131( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1132( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1133( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1134( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1135( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1136( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1137( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1138( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1139( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1140( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1141( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1142( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1143( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1144( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1145( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1146( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1147( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1148( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1149( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1150( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1151( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1152( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1153( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1154( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1155( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1156( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1157( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1158( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1159( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1160( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1161( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1162( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1163( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1164( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1165( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1166( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1167( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1168( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1169( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1170( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1171( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1172( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1173( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1174( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1175( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1176( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1177( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1178( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1179( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1180( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1181( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1182( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1183( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1184( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1185( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1186( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1187( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1188( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1189( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1190( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1191( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1192( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1193( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1194( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1195( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1196( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1197( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1198( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1199( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1200( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1201( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1202( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1203( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1204( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1205( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1206( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1207( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1208( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1209( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1210( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1211( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1212( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1213( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1214( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1215( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1216( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1217( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1218( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1219( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1220( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1221( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1222( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1223( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1224( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1225( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1226( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1227( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1228( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1229( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1230( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1231( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1232( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1233( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1234( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1235( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1236( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1237( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1238( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1239( wfs, jamps, channelIds, COUPs, numerators, denominators );
      diagram1240( wfs, jamps, channelIds, COUPs, numerators, denominators );
#endif
    }
    // *****************************
    // *** END LOOP ON IPARITY ***
    // *****************************

    return;
  }

  //--------------------------------------------------------------------------

  __global__ INLINE void /* clang-format off */
  color_sum( fptype* allMEs,              // output: allMEs[nevt], add |M|^2 for this specific helicity
#ifdef MGONGPUCPP_GPUIMPL
             const fptype* allJamps       // input: jamp[ncolor*2*nevt] for one specific helicity
#else
             const cxtype_sv* jamp_sv,    // input: jamp_sv[ncolor] (f/d) or [2*ncolor] (m) for SIMD event page(s) ievt00 and helicity ihel
             const int ievt00             // input: first event number in current C++ event page (for CUDA, ievt depends on threadid)
#endif
             ) /* clang-format on */
  {
#ifdef MGONGPUCPP_GPUIMPL
    //using namespace mg5amcGpu;
    using E_ACCESS = DeviceAccessMatrixElements; // non-trivial access: buffer includes all events
#else
    //using namespace mg5amcCpu;
    using E_ACCESS = HostAccessMatrixElements; // non-trivial access: buffer includes all events
#endif

    // *** COLOR MATRIX BELOW ***

    // The color denominators (initialize all array elements, with ncolor=120)
    // [NB do keep 'static' for these constexpr arrays, see issue #283]
    static constexpr fptype2 denom[ncolor] = { 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324, 324 }; // 1-D array[120]

    // The color matrix (initialize all array elements, with ncolor=120)
    // [NB do keep 'static' for these constexpr arrays, see issue #283]
    static constexpr fptype2 cf[ncolor][ncolor] = {
      { 4096, -512, -512, 64, 64, 640, -512, 64, 64, -8, -8, -80, 64, -8, 640, -80, 568, 496, -8, -80, -80, 496, 496, -224, -512, 64, 64, -8, -8, -80, 64, -8, -8, 1, 1, 10, -8, 1, -80, 10, -71, -62, 1, 10, 10, -62, -62, 28, 64, -8, -8, 1, 1, 10, 640, -80, -80, 10, 10, 100, 568, -71, 496, -62, 505, 514, -71, 19, -62, -53, -134, -44, -8, 1, -80, 10, -71, -62, -80, 10, 496, -62, 19, -53, 496, -62, -224, 28, -134, -44, 505, -134, -134, 442, 442, -116, 1, 10, 10, -62, -62, 28, 10, 100, -62, 514, -53, -44, -62, -53, 28, -44, 442, -116, 514, -44, -44, -116, -116, 136 },
      { -512, 4096, 64, 640, -512, 64, 64, -512, -8, -80, 64, -8, -8, -80, -80, 496, 496, -224, 64, -8, 640, -80, 568, 496, 64, -512, -8, -80, 64, -8, -8, 64, 1, 10, -8, 1, 1, 10, 10, -62, -62, 28, -8, 1, -80, 10, -71, -62, -8, 64, 1, 10, -8, 1, -80, 640, 10, 100, -80, 10, -71, 19, -62, -53, -134, -44, 568, -71, 496, -62, 505, 514, 1, 10, 10, -62, -62, 28, 10, 100, -62, 514, -53, -44, -62, -53, 28, -44, 442, -116, 514, -44, -44, -116, -116, 136, -8, 1, -80, 10, -71, -62, -80, 10, 496, -62, 19, -53, 496, -62, -224, 28, -134, -44, 505, -134, -134, 442, 442, -116 },
      { -512, 64, 4096, -512, 640, 64, 64, -8, 640, -80, 568, 496, -512, 64, 64, -8, -8, -80, -80, -8, 496, -224, -80, 496, 64, -8, -512, 64, -80, -8, -8, 1, -80, 10, -71, -62, 64, -8, -8, 1, 1, 10, 10, 1, -62, 28, 10, -62, -8, 1, -80, 10, -71, -62, -80, 10, 496, -62, 19, -53, 496, -62, -224, 28, -134, -44, 505, -134, -134, 442, 442, -116, 64, -8, -8, 1, 1, 10, 640, -80, -80, 10, 10, 100, 568, -71, 496, -62, 505, 514, -71, 19, -62, -53, -134, -44, 10, 1, -62, 28, 10, -62, 100, 10, -53, -44, -62, 514, 514, -44, -44, -116, -116, 136, -62, -53, 28, -44, 442, -116 },
      { 64, 640, -512, 4096, 64, -512, -8, -80, -80, 496, 496, -224, 64, -512, -8, -80, 64, -8, -8, 64, 568, 496, 640, -80, -8, -80, 64, -512, -8, 64, 1, 10, 10, -62, -62, 28, -8, 64, 1, 10, -8, 1, 1, -8, -71, -62, -80, 10, 1, 10, 10, -62, -62, 28, 10, 100, -62, 514, -53, -44, -62, -53, 28, -44, 442, -116, 514, -44, -44, -116, -116, 136, -8, 64, 1, 10, -8, 1, -80, 640, 10, 100, -80, 10, -71, 19, -62, -53, -134, -44, 568, -71, 496, -62, 505, 514, 1, -8, -71, -62, -80, 10, 10, -80, 19, -53, 496, -62, 505, -134, -134, 442, 442, -116, 496, -62, -224, 28, -134, -44 },
      { 64, -512, 640, 64, 4096, -512, -8, 64, 568, 496, 640, -80, -80, -8, 496, -224, -80, 496, -512, 64, 64, -8, -8, -80, -8, 64, -80, -8, -512, 64, 1, -8, -71, -62, -80, 10, 10, 1, -62, 28, 10, -62, 64, -8, -8, 1, 1, 10, 1, -8, -71, -62, -80, 10, 10, -80, 19, -53, 496, -62, 505, -134, -134, 442, 442, -116, 496, -62, -224, 28, -134, -44, 10, 1, -62, 28, 10, -62, 100, 10, -53, -44, -62, 514, 514, -44, -44, -116, -116, 136, -62, -53, 28, -44, 442, -116, 64, -8, -8, 1, 1, 10, 640, -80, -80, 10, 10, 100, 568, -71, 496, -62, 505, 514, -71, 19, -62, -53, -134, -44 },
      { 640, 64, 64, -512, -512, 4096, -80, -8, 496, -224, -80, 496, -8, 64, 568, 496, 640, -80, 64, -512, -8, -80, 64, -8, -80, -8, -8, 64, 64, -512, 10, 1, -62, 28, 10, -62, 1, -8, -71, -62, -80, 10, -8, 64, 1, 10, -8, 1, 10, 1, -62, 28, 10, -62, 100, 10, -53, -44, -62, 514, 514, -44, -44, -116, -116, 136, -62, -53, 28, -44, 442, -116, 1, -8, -71, -62, -80, 10, 10, -80, 19, -53, 496, -62, 505, -134, -134, 442, 442, -116, 496, -62, -224, 28, -134, -44, -8, 64, 1, 10, -8, 1, -80, 640, 10, 100, -80, 10, -71, 19, -62, -53, -134, -44, 568, -71, 496, -62, 505, 514 },
      { -512, 64, 64, -8, -8, -80, 4096, -512, -512, 64, 64, 640, 640, -80, 64, -8, 496, 568, -80, 496, -8, -80, -224, 496, 64, -8, -8, 1, 1, 10, 640, -80, -80, 10, 10, 100, 568, -71, 496, -62, 505, 514, -71, 19, -62, -53, -134, -44, -512, 64, 64, -8, -8, -80, 64, -8, -8, 1, 1, 10, -8, 1, -80, 10, -71, -62, 1, 10, 10, -62, -62, 28, -80, 10, -8, 1, -62, -71, 496, -62, -224, 28, -134, -44, -80, 10, 496, -62, 19, -53, -134, 505, 442, -116, -134, 442, 10, -62, 1, 10, 28, -62, -62, -53, 28, -44, 442, -116, 10, 100, -62, 514, -53, -44, -44, 514, -116, 136, -44, -116 },
      { 64, -512, -8, -80, 64, -8, -512, 4096, 64, 640, -512, 64, -80, 496, -8, -80, -224, 496, 640, -80, 64, -8, 496, 568, -8, 64, 1, 10, -8, 1, -80, 640, 10, 100, -80, 10, -71, 19, -62, -53, -134, -44, 568, -71, 496, -62, 505, 514, 64, -512, -8, -80, 64, -8, -8, 64, 1, 10, -8, 1, 1, 10, 10, -62, -62, 28, -8, 1, -80, 10, -71, -62, 10, -62, 1, 10, 28, -62, -62, -53, 28, -44, 442, -116, 10, 100, -62, 514, -53, -44, -44, 514, -116, 136, -44, -116, -80, 10, -8, 1, -62, -71, 496, -62, -224, 28, -134, -44, -80, 10, 496, -62, 19, -53, -134, 505, 442, -116, -134, 442 },
      { 64, -8, 640, -80, 568, 496, -512, 64, 4096, -512, 640, 64, 64, -8, -512, 64, -80, -8, 496, -224, -80, -8, 496, -80, -8, 1, -80, 10, -71, -62, -80, 10, 496, -62, 19, -53, 496, -62, -224, 28, -134, -44, 505, -134, -134, 442, 442, -116, 64, -8, -512, 64, -80, -8, -8, 1, -80, 10, -71, -62, 64, -8, -8, 1, 1, 10, 10, 1, -62, 28, 10, -62, -8, 1, 64, -8, 10, 1, 568, -71, 496, -62, 505, 514, 640, -80, -80, 10, 10, 100, 19, -71, -134, -44, -62, -53, -62, 28, 10, 1, -62, 10, 514, -44, -44, -116, -116, 136, 100, 10, -53, -44, -62, 514, -53, -62, 442, -116, 28, -44 },
      { -8, -80, -80, 496, 496, -224, 64, 640, -512, 4096, 64, -512, -8, -80, 64, -512, -8, 64, 568, 496, -8, 64, -80, 640, 1, 10, 10, -62, -62, 28, 10, 100, -62, 514, -53, -44, -62, -53, 28, -44, 442, -116, 514, -44, -44, -116, -116, 136, -8, -80, 64, -512, -8, 64, 1, 10, 10, -62, -62, 28, -8, 64, 1, 10, -8, 1, 1, -8, -71, -62, -80, 10, 1, 10, -8, 64, 1, -8, -71, 19, -62, -53, -134, -44, -80, 640, 10, 100, -80, 10, -71, 568, 505, 514, 496, -62, -71, -62, 1, -8, 10, -80, 505, -134, -134, 442, 442, -116, 10, -80, 19, -53, 496, -62, -62, 496, -134, -44, -224, 28 },
      { -8, 64, 568, 496, 640, -80, 64, -512, 640, 64, 4096, -512, 496, -224, -80, -8, 496, -80, 64, -8, -512, 64, -80, -8, 1, -8, -71, -62, -80, 10, 10, -80, 19, -53, 496, -62, 505, -134, -134, 442, 442, -116, 496, -62, -224, 28, -134, -44, -8, 64, -80, -8, -512, 64, 1, -8, -71, -62, -80, 10, 10, 1, -62, 28, 10, -62, 64, -8, -8, 1, 1, 10, -62, 28, 10, 1, -62, 10, 514, -44, -44, -116, -116, 136, 100, 10, -53, -44, -62, 514, -53, -62, 442, -116, 28, -44, -8, 1, 64, -8, 10, 1, 568, -71, 496, -62, 505, 514, 640, -80, -80, 10, 10, 100, 19, -71, -134, -44, -62, -53 },
      { -80, -8, 496, -224, -80, 496, 640, 64, 64, -512, -512, 4096, 568, 496, -8, 64, -80, 640, -8, -80, 64, -512, -8, 64, 10, 1, -62, 28, 10, -62, 100, 10, -53, -44, -62, 514, 514, -44, -44, -116, -116, 136, -62, -53, 28, -44, 442, -116, -80, -8, -8, 64, 64, -512, 10, 1, -62, 28, 10, -62, 1, -8, -71, -62, -80, 10, -8, 64, 1, 10, -8, 1, -71, -62, 1, -8, 10, -80, 505, -134, -134, 442, 442, -116, 10, -80, 19, -53, 496, -62, -62, 496, -134, -44, -224, 28, 1, 10, -8, 64, 1, -8, -71, 19, -62, -53, -134, -44, -80, 640, 10, 100, -80, 10, -71, 568, 505, 514, 496, -62 },
      { 64, -8, -512, 64, -80, -8, 640, -80, 64, -8, 496, 568, 4096, -512, -512, 64, 64, 640, 496, -80, -224, 496, -8, -80, -8, 1, 64, -8, 10, 1, 568, -71, 496, -62, 505, 514, 640, -80, -80, 10, 10, 100, 19, -71, -134, -44, -62, -53, -80, 10, -8, 1, -62, -71, 496, -62, -224, 28, -134, -44, -80, 10, 496, -62, 19, -53, -134, 505, 442, -116, -134, 442, -512, 64, 64, -8, -8, -80, 64, -8, -8, 1, 1, 10, -8, 1, -80, 10, -71, -62, 1, 10, 10, -62, -62, 28, -62, 10, 28, -62, 1, 10, -53, -62, 442, -116, 28, -44, -44, 514, -116, 136, -44, -116, 10, 100, -62, 514, -53, -44 },
      { -8, -80, 64, -512, -8, 64, -80, 496, -8, -80, -224, 496, -512, 4096, 64, 640, -512, 64, -80, 640, 496, 568, 64, -8, 1, 10, -8, 64, 1, -8, -71, 19, -62, -53, -134, -44, -80, 640, 10, 100, -80, 10, -71, 568, 505, 514, 496, -62, 10, -62, 1, 10, 28, -62, -62, -53, 28, -44, 442, -116, 10, 100, -62, 514, -53, -44, -44, 514, -116, 136, -44, -116, 64, -512, -8, -80, 64, -8, -8, 64, 1, 10, -8, 1, 1, 10, 10, -62, -62, 28, -8, 1, -80, 10, -71, -62, 10, -80, -62, -71, -8, 1, -62, 496, -134, -44, -224, 28, -134, 505, 442, -116, -134, 442, -80, 10, 496, -62, 19, -53 },
      { 640, -80, 64, -8, 496, 568, 64, -8, -512, 64, -80, -8, -512, 64, 4096, -512, 640, 64, -224, 496, 496, -80, -80, -8, -80, 10, -8, 1, -62, -71, 496, -62, -224, 28, -134, -44, -80, 10, 496, -62, 19, -53, -134, 505, 442, -116, -134, 442, -8, 1, 64, -8, 10, 1, 568, -71, 496, -62, 505, 514, 640, -80, -80, 10, 10, 100, 19, -71, -134, -44, -62, -53, 64, -8, -512, 64, -80, -8, -8, 1, -80, 10, -71, -62, 64, -8, -8, 1, 1, 10, 10, 1, -62, 28, 10, -62, 28, -62, -62, 10, 10, 1, -44, 514, -116, 136, -44, -116, -53, -62, 442, -116, 28, -44, 100, 10, -53, -44, -62, 514 },
      { -80, 496, -8, -80, -224, 496, -8, -80, 64, -512, -8, 64, 64, 640, -512, 4096, 64, -512, 496, 568, -80, 640, -8, 64, 10, -62, 1, 10, 28, -62, -62, -53, 28, -44, 442, -116, 10, 100, -62, 514, -53, -44, -44, 514, -116, 136, -44, -116, 1, 10, -8, 64, 1, -8, -71, 19, -62, -53, -134, -44, -80, 640, 10, 100, -80, 10, -71, 568, 505, 514, 496, -62, -8, -80, 64, -512, -8, 64, 1, 10, 10, -62, -62, 28, -8, 64, 1, 10, -8, 1, 1, -8, -71, -62, -80, 10, -62, -71, 10, -80, 1, -8, -134, 505, 442, -116, -134, 442, -62, 496, -134, -44, -224, 28, 10, -80, 19, -53, 496, -62 },
      { 568, 496, -8, 64, -80, 640, 496, -224, -80, -8, 496, -80, 64, -512, 640, 64, 4096, -512, -8, 64, -80, -8, -512, 64, -71, -62, 1, -8, 10, -80, 505, -134, -134, 442, 442, -116, 10, -80, 19, -53, 496, -62, -62, 496, -134, -44, -224, 28, -62, 28, 10, 1, -62, 10, 514, -44, -44, -116, -116, 136, 100, 10, -53, -44, -62, 514, -53, -62, 442, -116, 28, -44, -8, 64, -80, -8, -512, 64, 1, -8, -71, -62, -80, 10, 10, 1, -62, 28, 10, -62, 64, -8, -8, 1, 1, 10, 1, -8, 10, 1, 64, -8, -71, 568, 505, 514, 496, -62, 19, -71, -134, -44, -62, -53, 640, -80, -80, 10, 10, 100 },
      { 496, -224, -80, -8, 496, -80, 568, 496, -8, 64, -80, 640, 640, 64, 64, -512, -512, 4096, -80, -8, -8, 64, 64, -512, -62, 28, 10, 1, -62, 10, 514, -44, -44, -116, -116, 136, 100, 10, -53, -44, -62, 514, -53, -62, 442, -116, 28, -44, -71, -62, 1, -8, 10, -80, 505, -134, -134, 442, 442, -116, 10, -80, 19, -53, 496, -62, -62, 496, -134, -44, -224, 28, -80, -8, -8, 64, 64, -512, 10, 1, -62, 28, 10, -62, 1, -8, -71, -62, -80, 10, -8, 64, 1, 10, -8, 1, 10, 1, 1, -8, -8, 64, 19, -71, -134, -44, -62, -53, -71, 568, 505, 514, 496, -62, -80, 640, 10, 100, -80, 10 },
      { -8, 64, -80, -8, -512, 64, -80, 640, 496, 568, 64, -8, 496, -80, -224, 496, -8, -80, 4096, -512, -512, 64, 64, 640, 1, -8, 10, 1, 64, -8, -71, 568, 505, 514, 496, -62, 19, -71, -134, -44, -62, -53, 640, -80, -80, 10, 10, 100, 10, -80, -62, -71, -8, 1, -62, 496, -134, -44, -224, 28, -134, 505, 442, -116, -134, 442, -80, 10, 496, -62, 19, -53, -62, 10, 28, -62, 1, 10, -53, -62, 442, -116, 28, -44, -44, 514, -116, 136, -44, -116, 10, 100, -62, 514, -53, -44, -512, 64, 64, -8, -8, -80, 64, -8, -8, 1, 1, 10, -8, 1, -80, 10, -71, -62, 1, 10, 10, -62, -62, 28 },
      { -80, -8, -8, 64, 64, -512, 496, -80, -224, 496, -8, -80, -80, 640, 496, 568, 64, -8, -512, 4096, 64, 640, -512, 64, 10, 1, 1, -8, -8, 64, 19, -71, -134, -44, -62, -53, -71, 568, 505, 514, 496, -62, -80, 640, 10, 100, -80, 10, -62, 10, 28, -62, 1, 10, -53, -62, 442, -116, 28, -44, -44, 514, -116, 136, -44, -116, 10, 100, -62, 514, -53, -44, 10, -80, -62, -71, -8, 1, -62, 496, -134, -44, -224, 28, -134, 505, 442, -116, -134, 442, -80, 10, 496, -62, 19, -53, 64, -512, -8, -80, 64, -8, -8, 64, 1, 10, -8, 1, 1, 10, 10, -62, -62, 28, -8, 1, -80, 10, -71, -62 },
      { -80, 640, 496, 568, 64, -8, -8, 64, -80, -8, -512, 64, -224, 496, 496, -80, -80, -8, -512, 64, 4096, -512, 640, 64, 10, -80, -62, -71, -8, 1, -62, 496, -134, -44, -224, 28, -134, 505, 442, -116, -134, 442, -80, 10, 496, -62, 19, -53, 1, -8, 10, 1, 64, -8, -71, 568, 505, 514, 496, -62, 19, -71, -134, -44, -62, -53, 640, -80, -80, 10, 10, 100, 28, -62, -62, 10, 10, 1, -44, 514, -116, 136, -44, -116, -53, -62, 442, -116, 28, -44, 100, 10, -53, -44, -62, 514, 64, -8, -512, 64, -80, -8, -8, 1, -80, 10, -71, -62, 64, -8, -8, 1, 1, 10, 10, 1, -62, 28, 10, -62 },
      { 496, -80, -224, 496, -8, -80, -80, -8, -8, 64, 64, -512, 496, 568, -80, 640, -8, 64, 64, 640, -512, 4096, 64, -512, -62, 10, 28, -62, 1, 10, -53, -62, 442, -116, 28, -44, -44, 514, -116, 136, -44, -116, 10, 100, -62, 514, -53, -44, 10, 1, 1, -8, -8, 64, 19, -71, -134, -44, -62, -53, -71, 568, 505, 514, 496, -62, -80, 640, 10, 100, -80, 10, -62, -71, 10, -80, 1, -8, -134, 505, 442, -116, -134, 442, -62, 496, -134, -44, -224, 28, 10, -80, 19, -53, 496, -62, -8, -80, 64, -512, -8, 64, 1, 10, 10, -62, -62, 28, -8, 64, 1, 10, -8, 1, 1, -8, -71, -62, -80, 10 },
      { 496, 568, -80, 640, -8, 64, -224, 496, 496, -80, -80, -8, -8, 64, -80, -8, -512, 64, 64, -512, 640, 64, 4096, -512, -62, -71, 10, -80, 1, -8, -134, 505, 442, -116, -134, 442, -62, 496, -134, -44, -224, 28, 10, -80, 19, -53, 496, -62, 28, -62, -62, 10, 10, 1, -44, 514, -116, 136, -44, -116, -53, -62, 442, -116, 28, -44, 100, 10, -53, -44, -62, 514, 1, -8, 10, 1, 64, -8, -71, 568, 505, 514, 496, -62, 19, -71, -134, -44, -62, -53, 640, -80, -80, 10, 10, 100, -8, 64, -80, -8, -512, 64, 1, -8, -71, -62, -80, 10, 10, 1, -62, 28, 10, -62, 64, -8, -8, 1, 1, 10 },
      { -224, 496, 496, -80, -80, -8, 496, 568, -80, 640, -8, 64, -80, -8, -8, 64, 64, -512, 640, 64, 64, -512, -512, 4096, 28, -62, -62, 10, 10, 1, -44, 514, -116, 136, -44, -116, -53, -62, 442, -116, 28, -44, 100, 10, -53, -44, -62, 514, -62, -71, 10, -80, 1, -8, -134, 505, 442, -116, -134, 442, -62, 496, -134, -44, -224, 28, 10, -80, 19, -53, 496, -62, 10, 1, 1, -8, -8, 64, 19, -71, -134, -44, -62, -53, -71, 568, 505, 514, 496, -62, -80, 640, 10, 100, -80, 10, -80, -8, -8, 64, 64, -512, 10, 1, -62, 28, 10, -62, 1, -8, -71, -62, -80, 10, -8, 64, 1, 10, -8, 1 },
      { -512, 64, 64, -8, -8, -80, 64, -8, -8, 1, 1, 10, -8, 1, -80, 10, -71, -62, 1, 10, 10, -62, -62, 28, 4096, -512, -512, 64, 64, 640, -512, 64, 64, -8, -8, -80, 64, -8, 640, -80, 568, 496, -8, -80, -80, 496, 496, -224, 640, -80, -80, 10, 10, 100, 64, -8, -8, 1, 1, 10, 496, -62, 568, -71, 514, 505, -62, -53, -71, 19, -44, -134, -80, 10, 496, -62, 19, -53, -8, 1, -80, 10, -71, -62, -224, 28, 496, -62, -44, -134, -134, 442, 505, -134, -116, 442, 10, 100, -62, 514, -53, -44, 1, 10, 10, -62, -62, 28, 28, -44, -62, -53, -116, 442, -44, -116, 514, -44, 136, -116 },
      { 64, -512, -8, -80, 64, -8, -8, 64, 1, 10, -8, 1, 1, 10, 10, -62, -62, 28, -8, 1, -80, 10, -71, -62, -512, 4096, 64, 640, -512, 64, 64, -512, -8, -80, 64, -8, -8, -80, -80, 496, 496, -224, 64, -8, 640, -80, 568, 496, -80, 640, 10, 100, -80, 10, -8, 64, 1, 10, -8, 1, -62, -53, -71, 19, -44, -134, 496, -62, 568, -71, 514, 505, 10, 100, -62, 514, -53, -44, 1, 10, 10, -62, -62, 28, 28, -44, -62, -53, -116, 442, -44, -116, 514, -44, 136, -116, -80, 10, 496, -62, 19, -53, -8, 1, -80, 10, -71, -62, -224, 28, 496, -62, -44, -134, -134, 442, 505, -134, -116, 442 },
      { 64, -8, -512, 64, -80, -8, -8, 1, -80, 10, -71, -62, 64, -8, -8, 1, 1, 10, 10, 1, -62, 28, 10, -62, -512, 64, 4096, -512, 640, 64, 64, -8, 640, -80, 568, 496, -512, 64, 64, -8, -8, -80, -80, -8, 496, -224, -80, 496, -80, 10, 496, -62, 19, -53, -8, 1, -80, 10, -71, -62, -224, 28, 496, -62, -44, -134, -134, 442, 505, -134, -116, 442, 640, -80, -80, 10, 10, 100, 64, -8, -8, 1, 1, 10, 496, -62, 568, -71, 514, 505, -62, -53, -71, 19, -44, -134, 100, 10, -53, -44, -62, 514, 10, 1, -62, 28, 10, -62, -44, -116, 514, -44, 136, -116, 28, -44, -62, -53, -116, 442 },
      { -8, -80, 64, -512, -8, 64, 1, 10, 10, -62, -62, 28, -8, 64, 1, 10, -8, 1, 1, -8, -71, -62, -80, 10, 64, 640, -512, 4096, 64, -512, -8, -80, -80, 496, 496, -224, 64, -512, -8, -80, 64, -8, -8, 64, 568, 496, 640, -80, 10, 100, -62, 514, -53, -44, 1, 10, 10, -62, -62, 28, 28, -44, -62, -53, -116, 442, -44, -116, 514, -44, 136, -116, -80, 640, 10, 100, -80, 10, -8, 64, 1, 10, -8, 1, -62, -53, -71, 19, -44, -134, 496, -62, 568, -71, 514, 505, 10, -80, 19, -53, 496, -62, 1, -8, -71, -62, -80, 10, -134, 442, 505, -134, -116, 442, -224, 28, 496, -62, -44, -134 },
      { -8, 64, -80, -8, -512, 64, 1, -8, -71, -62, -80, 10, 10, 1, -62, 28, 10, -62, 64, -8, -8, 1, 1, 10, 64, -512, 640, 64, 4096, -512, -8, 64, 568, 496, 640, -80, -80, -8, 496, -224, -80, 496, -512, 64, 64, -8, -8, -80, 10, -80, 19, -53, 496, -62, 1, -8, -71, -62, -80, 10, -134, 442, 505, -134, -116, 442, -224, 28, 496, -62, -44, -134, 100, 10, -53, -44, -62, 514, 10, 1, -62, 28, 10, -62, -44, -116, 514, -44, 136, -116, 28, -44, -62, -53, -116, 442, 640, -80, -80, 10, 10, 100, 64, -8, -8, 1, 1, 10, 496, -62, 568, -71, 514, 505, -62, -53, -71, 19, -44, -134 },
      { -80, -8, -8, 64, 64, -512, 10, 1, -62, 28, 10, -62, 1, -8, -71, -62, -80, 10, -8, 64, 1, 10, -8, 1, 640, 64, 64, -512, -512, 4096, -80, -8, 496, -224, -80, 496, -8, 64, 568, 496, 640, -80, 64, -512, -8, -80, 64, -8, 100, 10, -53, -44, -62, 514, 10, 1, -62, 28, 10, -62, -44, -116, 514, -44, 136, -116, 28, -44, -62, -53, -116, 442, 10, -80, 19, -53, 496, -62, 1, -8, -71, -62, -80, 10, -134, 442, 505, -134, -116, 442, -224, 28, 496, -62, -44, -134, -80, 640, 10, 100, -80, 10, -8, 64, 1, 10, -8, 1, -62, -53, -71, 19, -44, -134, 496, -62, 568, -71, 514, 505 },
      { 64, -8, -8, 1, 1, 10, 640, -80, -80, 10, 10, 100, 568, -71, 496, -62, 505, 514, -71, 19, -62, -53, -134, -44, -512, 64, 64, -8, -8, -80, 4096, -512, -512, 64, 64, 640, 640, -80, 64, -8, 496, 568, -80, 496, -8, -80, -224, 496, 64, -8, -8, 1, 1, 10, -512, 64, 64, -8, -8, -80, -80, 10, -8, 1, -62, -71, 10, -62, 1, 10, 28, -62, 496, -62, -224, 28, -134, -44, -80, 10, -8, 1, -62, -71, 496, -62, -80, 10, -53, 19, 442, -116, -134, 505, 442, -134, -62, -53, 28, -44, 442, -116, 10, -62, 1, 10, 28, -62, -62, 514, 10, 100, -44, -53, -116, 136, -44, 514, -116, -44 },
      { -8, 64, 1, 10, -8, 1, -80, 640, 10, 100, -80, 10, -71, 19, -62, -53, -134, -44, 568, -71, 496, -62, 505, 514, 64, -512, -8, -80, 64, -8, -512, 4096, 64, 640, -512, 64, -80, 496, -8, -80, -224, 496, 640, -80, 64, -8, 496, 568, -8, 64, 1, 10, -8, 1, 64, -512, -8, -80, 64, -8, 10, -62, 1, 10, 28, -62, -80, 10, -8, 1, -62, -71, -62, -53, 28, -44, 442, -116, 10, -62, 1, 10, 28, -62, -62, 514, 10, 100, -44, -53, -116, 136, -44, 514, -116, -44, 496, -62, -224, 28, -134, -44, -80, 10, -8, 1, -62, -71, 496, -62, -80, 10, -53, 19, 442, -116, -134, 505, 442, -134 },
      { -8, 1, -80, 10, -71, -62, -80, 10, 496, -62, 19, -53, 496, -62, -224, 28, -134, -44, 505, -134, -134, 442, 442, -116, 64, -8, 640, -80, 568, 496, -512, 64, 4096, -512, 640, 64, 64, -8, -512, 64, -80, -8, 496, -224, -80, -8, 496, -80, -8, 1, -80, 10, -71, -62, 64, -8, -512, 64, -80, -8, -8, 1, 64, -8, 10, 1, -62, 28, 10, 1, -62, 10, 568, -71, 496, -62, 505, 514, -8, 1, 64, -8, 10, 1, -80, 10, 640, -80, 100, 10, -134, -44, 19, -71, -53, -62, 514, -44, -44, -116, -116, 136, -62, 28, 10, 1, -62, 10, -53, -44, 100, 10, 514, -62, 442, -116, -53, -62, -44, 28 },
      { 1, 10, 10, -62, -62, 28, 10, 100, -62, 514, -53, -44, -62, -53, 28, -44, 442, -116, 514, -44, -44, -116, -116, 136, -8, -80, -80, 496, 496, -224, 64, 640, -512, 4096, 64, -512, -8, -80, 64, -512, -8, 64, 568, 496, -8, 64, -80, 640, 1, 10, 10, -62, -62, 28, -8, -80, 64, -512, -8, 64, 1, 10, -8, 64, 1, -8, -71, -62, 1, -8, 10, -80, -71, 19, -62, -53, -134, -44, 1, 10, -8, 64, 1, -8, 10, 100, -80, 640, 10, -80, 505, 514, -71, 568, -62, 496, 505, -134, -134, 442, 442, -116, -71, -62, 1, -8, 10, -80, 19, -53, 10, -80, -62, 496, -134, -44, -62, 496, 28, -224 },
      { 1, -8, -71, -62, -80, 10, 10, -80, 19, -53, 496, -62, 505, -134, -134, 442, 442, -116, 496, -62, -224, 28, -134, -44, -8, 64, 568, 496, 640, -80, 64, -512, 640, 64, 4096, -512, 496, -224, -80, -8, 496, -80, 64, -8, -512, 64, -80, -8, 1, -8, -71, -62, -80, 10, -8, 64, -80, -8, -512, 64, -62, 28, 10, 1, -62, 10, -8, 1, 64, -8, 10, 1, 514, -44, -44, -116, -116, 136, -62, 28, 10, 1, -62, 10, -53, -44, 100, 10, 514, -62, 442, -116, -53, -62, -44, 28, 568, -71, 496, -62, 505, 514, -8, 1, 64, -8, 10, 1, -80, 10, 640, -80, 100, 10, -134, -44, 19, -71, -53, -62 },
      { 10, 1, -62, 28, 10, -62, 100, 10, -53, -44, -62, 514, 514, -44, -44, -116, -116, 136, -62, -53, 28, -44, 442, -116, -80, -8, 496, -224, -80, 496, 640, 64, 64, -512, -512, 4096, 568, 496, -8, 64, -80, 640, -8, -80, 64, -512, -8, 64, 10, 1, -62, 28, 10, -62, -80, -8, -8, 64, 64, -512, -71, -62, 1, -8, 10, -80, 1, 10, -8, 64, 1, -8, 505, -134, -134, 442, 442, -116, -71, -62, 1, -8, 10, -80, 19, -53, 10, -80, -62, 496, -134, -44, -62, 496, 28, -224, -71, 19, -62, -53, -134, -44, 1, 10, -8, 64, 1, -8, 10, 100, -80, 640, 10, -80, 505, 514, -71, 568, -62, 496 },
      { -8, 1, 64, -8, 10, 1, 568, -71, 496, -62, 505, 514, 640, -80, -80, 10, 10, 100, 19, -71, -134, -44, -62, -53, 64, -8, -512, 64, -80, -8, 640, -80, 64, -8, 496, 568, 4096, -512, -512, 64, 64, 640, 496, -80, -224, 496, -8, -80, 496, -62, -224, 28, -134, -44, -80, 10, -8, 1, -62, -71, 496, -62, -80, 10, -53, 19, 442, -116, -134, 505, 442, -134, 64, -8, -8, 1, 1, 10, -512, 64, 64, -8, -8, -80, -80, 10, -8, 1, -62, -71, 10, -62, 1, 10, 28, -62, -53, -62, 442, -116, 28, -44, -62, 10, 28, -62, 1, 10, -116, 136, -44, 514, -116, -44, -62, 514, 10, 100, -44, -53 },
      { 1, 10, -8, 64, 1, -8, -71, 19, -62, -53, -134, -44, -80, 640, 10, 100, -80, 10, -71, 568, 505, 514, 496, -62, -8, -80, 64, -512, -8, 64, -80, 496, -8, -80, -224, 496, -512, 4096, 64, 640, -512, 64, -80, 640, 496, 568, 64, -8, -62, -53, 28, -44, 442, -116, 10, -62, 1, 10, 28, -62, -62, 514, 10, 100, -44, -53, -116, 136, -44, 514, -116, -44, -8, 64, 1, 10, -8, 1, 64, -512, -8, -80, 64, -8, 10, -62, 1, 10, 28, -62, -80, 10, -8, 1, -62, -71, -62, 496, -134, -44, -224, 28, 10, -80, -62, -71, -8, 1, 442, -116, -134, 505, 442, -134, 496, -62, -80, 10, -53, 19 },
      { -80, 10, -8, 1, -62, -71, 496, -62, -224, 28, -134, -44, -80, 10, 496, -62, 19, -53, -134, 505, 442, -116, -134, 442, 640, -80, 64, -8, 496, 568, 64, -8, -512, 64, -80, -8, -512, 64, 4096, -512, 640, 64, -224, 496, 496, -80, -80, -8, 568, -71, 496, -62, 505, 514, -8, 1, 64, -8, 10, 1, -80, 10, 640, -80, 100, 10, -134, -44, 19, -71, -53, -62, -8, 1, -80, 10, -71, -62, 64, -8, -512, 64, -80, -8, -8, 1, 64, -8, 10, 1, -62, 28, 10, 1, -62, 10, -44, 514, -116, 136, -44, -116, 28, -62, -62, 10, 10, 1, 442, -116, -53, -62, -44, 28, -53, -44, 100, 10, 514, -62 },
      { 10, -62, 1, 10, 28, -62, -62, -53, 28, -44, 442, -116, 10, 100, -62, 514, -53, -44, -44, 514, -116, 136, -44, -116, -80, 496, -8, -80, -224, 496, -8, -80, 64, -512, -8, 64, 64, 640, -512, 4096, 64, -512, 496, 568, -80, 640, -8, 64, -71, 19, -62, -53, -134, -44, 1, 10, -8, 64, 1, -8, 10, 100, -80, 640, 10, -80, 505, 514, -71, 568, -62, 496, 1, 10, 10, -62, -62, 28, -8, -80, 64, -512, -8, 64, 1, 10, -8, 64, 1, -8, -71, -62, 1, -8, 10, -80, -134, 505, 442, -116, -134, 442, -62, -71, 10, -80, 1, -8, -134, -44, -62, 496, 28, -224, 19, -53, 10, -80, -62, 496 },
      { -71, -62, 1, -8, 10, -80, 505, -134, -134, 442, 442, -116, 10, -80, 19, -53, 496, -62, -62, 496, -134, -44, -224, 28, 568, 496, -8, 64, -80, 640, 496, -224, -80, -8, 496, -80, 64, -512, 640, 64, 4096, -512, -8, 64, -80, -8, -512, 64, 514, -44, -44, -116, -116, 136, -62, 28, 10, 1, -62, 10, -53, -44, 100, 10, 514, -62, 442, -116, -53, -62, -44, 28, 1, -8, -71, -62, -80, 10, -8, 64, -80, -8, -512, 64, -62, 28, 10, 1, -62, 10, -8, 1, 64, -8, 10, 1, -71, 568, 505, 514, 496, -62, 1, -8, 10, 1, 64, -8, -134, -44, 19, -71, -53, -62, -80, 10, 640, -80, 100, 10 },
      { -62, 28, 10, 1, -62, 10, 514, -44, -44, -116, -116, 136, 100, 10, -53, -44, -62, 514, -53, -62, 442, -116, 28, -44, 496, -224, -80, -8, 496, -80, 568, 496, -8, 64, -80, 640, 640, 64, 64, -512, -512, 4096, -80, -8, -8, 64, 64, -512, 505, -134, -134, 442, 442, -116, -71, -62, 1, -8, 10, -80, 19, -53, 10, -80, -62, 496, -134, -44, -62, 496, 28, -224, 10, 1, -62, 28, 10, -62, -80, -8, -8, 64, 64, -512, -71, -62, 1, -8, 10, -80, 1, 10, -8, 64, 1, -8, 19, -71, -134, -44, -62, -53, 10, 1, 1, -8, -8, 64, 505, 514, -71, 568, -62, 496, 10, 100, -80, 640, 10, -80 },
      { 1, -8, 10, 1, 64, -8, -71, 568, 505, 514, 496, -62, 19, -71, -134, -44, -62, -53, 640, -80, -80, 10, 10, 100, -8, 64, -80, -8, -512, 64, -80, 640, 496, 568, 64, -8, 496, -80, -224, 496, -8, -80, 4096, -512, -512, 64, 64, 640, -62, 496, -134, -44, -224, 28, 10, -80, -62, -71, -8, 1, 442, -116, -134, 505, 442, -134, 496, -62, -80, 10, -53, 19, -53, -62, 442, -116, 28, -44, -62, 10, 28, -62, 1, 10, -116, 136, -44, 514, -116, -44, -62, 514, 10, 100, -44, -53, 64, -8, -8, 1, 1, 10, -512, 64, 64, -8, -8, -80, -80, 10, -8, 1, -62, -71, 10, -62, 1, 10, 28, -62 },
      { 10, 1, 1, -8, -8, 64, 19, -71, -134, -44, -62, -53, -71, 568, 505, 514, 496, -62, -80, 640, 10, 100, -80, 10, -80, -8, -8, 64, 64, -512, 496, -80, -224, 496, -8, -80, -80, 640, 496, 568, 64, -8, -512, 4096, 64, 640, -512, 64, -53, -62, 442, -116, 28, -44, -62, 10, 28, -62, 1, 10, -116, 136, -44, 514, -116, -44, -62, 514, 10, 100, -44, -53, -62, 496, -134, -44, -224, 28, 10, -80, -62, -71, -8, 1, 442, -116, -134, 505, 442, -134, 496, -62, -80, 10, -53, 19, -8, 64, 1, 10, -8, 1, 64, -512, -8, -80, 64, -8, 10, -62, 1, 10, 28, -62, -80, 10, -8, 1, -62, -71 },
      { 10, -80, -62, -71, -8, 1, -62, 496, -134, -44, -224, 28, -134, 505, 442, -116, -134, 442, -80, 10, 496, -62, 19, -53, -80, 640, 496, 568, 64, -8, -8, 64, -80, -8, -512, 64, -224, 496, 496, -80, -80, -8, -512, 64, 4096, -512, 640, 64, -71, 568, 505, 514, 496, -62, 1, -8, 10, 1, 64, -8, -134, -44, 19, -71, -53, -62, -80, 10, 640, -80, 100, 10, -44, 514, -116, 136, -44, -116, 28, -62, -62, 10, 10, 1, 442, -116, -53, -62, -44, 28, -53, -44, 100, 10, 514, -62, -8, 1, -80, 10, -71, -62, 64, -8, -512, 64, -80, -8, -8, 1, 64, -8, 10, 1, -62, 28, 10, 1, -62, 10 },
      { -62, 10, 28, -62, 1, 10, -53, -62, 442, -116, 28, -44, -44, 514, -116, 136, -44, -116, 10, 100, -62, 514, -53, -44, 496, -80, -224, 496, -8, -80, -80, -8, -8, 64, 64, -512, 496, 568, -80, 640, -8, 64, 64, 640, -512, 4096, 64, -512, 19, -71, -134, -44, -62, -53, 10, 1, 1, -8, -8, 64, 505, 514, -71, 568, -62, 496, 10, 100, -80, 640, 10, -80, -134, 505, 442, -116, -134, 442, -62, -71, 10, -80, 1, -8, -134, -44, -62, 496, 28, -224, 19, -53, 10, -80, -62, 496, 1, 10, 10, -62, -62, 28, -8, -80, 64, -512, -8, 64, 1, 10, -8, 64, 1, -8, -71, -62, 1, -8, 10, -80 },
      { -62, -71, 10, -80, 1, -8, -134, 505, 442, -116, -134, 442, -62, 496, -134, -44, -224, 28, 10, -80, 19, -53, 496, -62, 496, 568, -80, 640, -8, 64, -224, 496, 496, -80, -80, -8, -8, 64, -80, -8, -512, 64, 64, -512, 640, 64, 4096, -512, -44, 514, -116, 136, -44, -116, 28, -62, -62, 10, 10, 1, 442, -116, -53, -62, -44, 28, -53, -44, 100, 10, 514, -62, -71, 568, 505, 514, 496, -62, 1, -8, 10, 1, 64, -8, -134, -44, 19, -71, -53, -62, -80, 10, 640, -80, 100, 10, 1, -8, -71, -62, -80, 10, -8, 64, -80, -8, -512, 64, -62, 28, 10, 1, -62, 10, -8, 1, 64, -8, 10, 1 },
      { 28, -62, -62, 10, 10, 1, -44, 514, -116, 136, -44, -116, -53, -62, 442, -116, 28, -44, 100, 10, -53, -44, -62, 514, -224, 496, 496, -80, -80, -8, 496, 568, -80, 640, -8, 64, -80, -8, -8, 64, 64, -512, 640, 64, 64, -512, -512, 4096, -134, 505, 442, -116, -134, 442, -62, -71, 10, -80, 1, -8, -134, -44, -62, 496, 28, -224, 19, -53, 10, -80, -62, 496, 19, -71, -134, -44, -62, -53, 10, 1, 1, -8, -8, 64, 505, 514, -71, 568, -62, 496, 10, 100, -80, 640, 10, -80, 10, 1, -62, 28, 10, -62, -80, -8, -8, 64, 64, -512, -71, -62, 1, -8, 10, -80, 1, 10, -8, 64, 1, -8 },
      { 64, -8, -8, 1, 1, 10, -512, 64, 64, -8, -8, -80, -80, 10, -8, 1, -62, -71, 10, -62, 1, 10, 28, -62, 640, -80, -80, 10, 10, 100, 64, -8, -8, 1, 1, 10, 496, -62, 568, -71, 514, 505, -62, -53, -71, 19, -44, -134, 4096, -512, -512, 64, 64, 640, -512, 64, 64, -8, -8, -80, 64, -8, 640, -80, 568, 496, -8, -80, -80, 496, 496, -224, 496, -62, -80, 10, -53, 19, -224, 28, 496, -62, -44, -134, -8, 1, -80, 10, -71, -62, 442, -134, -116, 442, 505, -134, -62, 514, 10, 100, -44, -53, 28, -44, -62, -53, -116, 442, 1, 10, 10, -62, -62, 28, -116, -44, 136, -116, 514, -44 },
      { -8, 64, 1, 10, -8, 1, 64, -512, -8, -80, 64, -8, 10, -62, 1, 10, 28, -62, -80, 10, -8, 1, -62, -71, -80, 640, 10, 100, -80, 10, -8, 64, 1, 10, -8, 1, -62, -53, -71, 19, -44, -134, 496, -62, 568, -71, 514, 505, -512, 4096, 64, 640, -512, 64, 64, -512, -8, -80, 64, -8, -8, -80, -80, 496, 496, -224, 64, -8, 640, -80, 568, 496, -62, 514, 10, 100, -44, -53, 28, -44, -62, -53, -116, 442, 1, 10, 10, -62, -62, 28, -116, -44, 136, -116, 514, -44, 496, -62, -80, 10, -53, 19, -224, 28, 496, -62, -44, -134, -8, 1, -80, 10, -71, -62, 442, -134, -116, 442, 505, -134 },
      { -8, 1, -80, 10, -71, -62, 64, -8, -512, 64, -80, -8, -8, 1, 64, -8, 10, 1, -62, 28, 10, 1, -62, 10, -80, 10, 496, -62, 19, -53, -8, 1, -80, 10, -71, -62, -224, 28, 496, -62, -44, -134, -134, 442, 505, -134, -116, 442, -512, 64, 4096, -512, 640, 64, 64, -8, 640, -80, 568, 496, -512, 64, 64, -8, -8, -80, -80, -8, 496, -224, -80, 496, -80, 10, 640, -80, 100, 10, 496, -62, 568, -71, 514, 505, 64, -8, -8, 1, 1, 10, -53, -62, -44, -134, -71, 19, -53, -44, 100, 10, 514, -62, -44, -116, 514, -44, 136, -116, 10, 1, -62, 28, 10, -62, -44, 28, -116, 442, -62, -53 },
      { 1, 10, 10, -62, -62, 28, -8, -80, 64, -512, -8, 64, 1, 10, -8, 64, 1, -8, -71, -62, 1, -8, 10, -80, 10, 100, -62, 514, -53, -44, 1, 10, 10, -62, -62, 28, 28, -44, -62, -53, -116, 442, -44, -116, 514, -44, 136, -116, 64, 640, -512, 4096, 64, -512, -8, -80, -80, 496, 496, -224, 64, -512, -8, -80, 64, -8, -8, 64, 568, 496, 640, -80, 10, 100, -80, 640, 10, -80, -62, -53, -71, 19, -44, -134, -8, 64, 1, 10, -8, 1, -62, 496, 514, 505, 568, -71, 19, -53, 10, -80, -62, 496, -134, 442, 505, -134, -116, 442, 1, -8, -71, -62, -80, 10, 28, -224, -44, -134, 496, -62 },
      { 1, -8, -71, -62, -80, 10, -8, 64, -80, -8, -512, 64, -62, 28, 10, 1, -62, 10, -8, 1, 64, -8, 10, 1, 10, -80, 19, -53, 496, -62, 1, -8, -71, -62, -80, 10, -134, 442, 505, -134, -116, 442, -224, 28, 496, -62, -44, -134, 64, -512, 640, 64, 4096, -512, -8, 64, 568, 496, 640, -80, -80, -8, 496, -224, -80, 496, -512, 64, 64, -8, -8, -80, -53, -44, 100, 10, 514, -62, -44, -116, 514, -44, 136, -116, 10, 1, -62, 28, 10, -62, -44, 28, -116, 442, -62, -53, -80, 10, 640, -80, 100, 10, 496, -62, 568, -71, 514, 505, 64, -8, -8, 1, 1, 10, -53, -62, -44, -134, -71, 19 },
      { 10, 1, -62, 28, 10, -62, -80, -8, -8, 64, 64, -512, -71, -62, 1, -8, 10, -80, 1, 10, -8, 64, 1, -8, 100, 10, -53, -44, -62, 514, 10, 1, -62, 28, 10, -62, -44, -116, 514, -44, 136, -116, 28, -44, -62, -53, -116, 442, 640, 64, 64, -512, -512, 4096, -80, -8, 496, -224, -80, 496, -8, 64, 568, 496, 640, -80, 64, -512, -8, -80, 64, -8, 19, -53, 10, -80, -62, 496, -134, 442, 505, -134, -116, 442, 1, -8, -71, -62, -80, 10, 28, -224, -44, -134, 496, -62, 10, 100, -80, 640, 10, -80, -62, -53, -71, 19, -44, -134, -8, 64, 1, 10, -8, 1, -62, 496, 514, 505, 568, -71 },
      { 640, -80, -80, 10, 10, 100, 64, -8, -8, 1, 1, 10, 496, -62, 568, -71, 514, 505, -62, -53, -71, 19, -44, -134, 64, -8, -8, 1, 1, 10, -512, 64, 64, -8, -8, -80, -80, 10, -8, 1, -62, -71, 10, -62, 1, 10, 28, -62, -512, 64, 64, -8, -8, -80, 4096, -512, -512, 64, 64, 640, 640, -80, 64, -8, 496, 568, -80, 496, -8, -80, -224, 496, -224, 28, 496, -62, -44, -134, 496, -62, -80, 10, -53, 19, -80, 10, -8, 1, -62, -71, -116, 442, 442, -134, -134, 505, 28, -44, -62, -53, -116, 442, -62, 514, 10, 100, -44, -53, 10, -62, 1, 10, 28, -62, 136, -116, -116, -44, -44, 514 },
      { -80, 640, 10, 100, -80, 10, -8, 64, 1, 10, -8, 1, -62, -53, -71, 19, -44, -134, 496, -62, 568, -71, 514, 505, -8, 64, 1, 10, -8, 1, 64, -512, -8, -80, 64, -8, 10, -62, 1, 10, 28, -62, -80, 10, -8, 1, -62, -71, 64, -512, -8, -80, 64, -8, -512, 4096, 64, 640, -512, 64, -80, 496, -8, -80, -224, 496, 640, -80, 64, -8, 496, 568, 28, -44, -62, -53, -116, 442, -62, 514, 10, 100, -44, -53, 10, -62, 1, 10, 28, -62, 136, -116, -116, -44, -44, 514, -224, 28, 496, -62, -44, -134, 496, -62, -80, 10, -53, 19, -80, 10, -8, 1, -62, -71, -116, 442, 442, -134, -134, 505 },
      { -80, 10, 496, -62, 19, -53, -8, 1, -80, 10, -71, -62, -224, 28, 496, -62, -44, -134, -134, 442, 505, -134, -116, 442, -8, 1, -80, 10, -71, -62, 64, -8, -512, 64, -80, -8, -8, 1, 64, -8, 10, 1, -62, 28, 10, 1, -62, 10, 64, -8, 640, -80, 568, 496, -512, 64, 4096, -512, 640, 64, 64, -8, -512, 64, -80, -8, 496, -224, -80, -8, 496, -80, 496, -62, 568, -71, 514, 505, -80, 10, 640, -80, 100, 10, -8, 1, 64, -8, 10, 1, -44, -134, -53, -62, 19, -71, -44, -116, 514, -44, 136, -116, -53, -44, 100, 10, 514, -62, -62, 28, 10, 1, -62, 10, -116, 442, -44, 28, -53, -62 },
      { 10, 100, -62, 514, -53, -44, 1, 10, 10, -62, -62, 28, 28, -44, -62, -53, -116, 442, -44, -116, 514, -44, 136, -116, 1, 10, 10, -62, -62, 28, -8, -80, 64, -512, -8, 64, 1, 10, -8, 64, 1, -8, -71, -62, 1, -8, 10, -80, -8, -80, -80, 496, 496, -224, 64, 640, -512, 4096, 64, -512, -8, -80, 64, -512, -8, 64, 568, 496, -8, 64, -80, 640, -62, -53, -71, 19, -44, -134, 10, 100, -80, 640, 10, -80, 1, 10, -8, 64, 1, -8, 514, 505, -62, 496, -71, 568, -134, 442, 505, -134, -116, 442, 19, -53, 10, -80, -62, 496, -71, -62, 1, -8, 10, -80, -44, -134, 28, -224, -62, 496 },
      { 10, -80, 19, -53, 496, -62, 1, -8, -71, -62, -80, 10, -134, 442, 505, -134, -116, 442, -224, 28, 496, -62, -44, -134, 1, -8, -71, -62, -80, 10, -8, 64, -80, -8, -512, 64, -62, 28, 10, 1, -62, 10, -8, 1, 64, -8, 10, 1, -8, 64, 568, 496, 640, -80, 64, -512, 640, 64, 4096, -512, 496, -224, -80, -8, 496, -80, 64, -8, -512, 64, -80, -8, -44, -116, 514, -44, 136, -116, -53, -44, 100, 10, 514, -62, -62, 28, 10, 1, -62, 10, -116, 442, -44, 28, -53, -62, 496, -62, 568, -71, 514, 505, -80, 10, 640, -80, 100, 10, -8, 1, 64, -8, 10, 1, -44, -134, -53, -62, 19, -71 },
      { 100, 10, -53, -44, -62, 514, 10, 1, -62, 28, 10, -62, -44, -116, 514, -44, 136, -116, 28, -44, -62, -53, -116, 442, 10, 1, -62, 28, 10, -62, -80, -8, -8, 64, 64, -512, -71, -62, 1, -8, 10, -80, 1, 10, -8, 64, 1, -8, -80, -8, 496, -224, -80, 496, 640, 64, 64, -512, -512, 4096, 568, 496, -8, 64, -80, 640, -8, -80, 64, -512, -8, 64, -134, 442, 505, -134, -116, 442, 19, -53, 10, -80, -62, 496, -71, -62, 1, -8, 10, -80, -44, -134, 28, -224, -62, 496, -62, -53, -71, 19, -44, -134, 10, 100, -80, 640, 10, -80, 1, 10, -8, 64, 1, -8, 514, 505, -62, 496, -71, 568 },
      { 568, -71, 496, -62, 505, 514, -8, 1, 64, -8, 10, 1, -80, 10, 640, -80, 100, 10, -134, -44, 19, -71, -53, -62, 496, -62, -224, 28, -134, -44, -80, 10, -8, 1, -62, -71, 496, -62, -80, 10, -53, 19, 442, -116, -134, 505, 442, -134, 64, -8, -512, 64, -80, -8, 640, -80, 64, -8, 496, 568, 4096, -512, -512, 64, 64, 640, 496, -80, -224, 496, -8, -80, -8, 1, 64, -8, 10, 1, -80, 10, -8, 1, -62, -71, -512, 64, 64, -8, -8, -80, -62, 10, 28, -62, 1, 10, 442, -116, -53, -62, -44, 28, -116, 136, -44, 514, -116, -44, -62, 10, 28, -62, 1, 10, 514, -62, -44, -53, 10, 100 },
      { -71, 19, -62, -53, -134, -44, 1, 10, -8, 64, 1, -8, 10, 100, -80, 640, 10, -80, 505, 514, -71, 568, -62, 496, -62, -53, 28, -44, 442, -116, 10, -62, 1, 10, 28, -62, -62, 514, 10, 100, -44, -53, -116, 136, -44, 514, -116, -44, -8, -80, 64, -512, -8, 64, -80, 496, -8, -80, -224, 496, -512, 4096, 64, 640, -512, 64, -80, 640, 496, 568, 64, -8, 1, 10, -8, 64, 1, -8, 10, -62, 1, 10, 28, -62, 64, -512, -8, -80, 64, -8, 10, -80, -62, -71, -8, 1, -134, -44, -62, 496, 28, -224, 442, -116, -134, 505, 442, -134, 10, -80, -62, -71, -8, 1, -62, 496, -53, 19, -80, 10 },
      { 496, -62, -224, 28, -134, -44, -80, 10, -8, 1, -62, -71, 496, -62, -80, 10, -53, 19, 442, -116, -134, 505, 442, -134, 568, -71, 496, -62, 505, 514, -8, 1, 64, -8, 10, 1, -80, 10, 640, -80, 100, 10, -134, -44, 19, -71, -53, -62, 640, -80, 64, -8, 496, 568, 64, -8, -512, 64, -80, -8, -512, 64, 4096, -512, 640, 64, -224, 496, 496, -80, -80, -8, -80, 10, -8, 1, -62, -71, -8, 1, 64, -8, 10, 1, 64, -8, -512, 64, -80, -8, 28, -62, -62, 10, 10, 1, -116, 136, -44, 514, -116, -44, 442, -116, -53, -62, -44, 28, 28, -62, -62, 10, 10, 1, -44, -53, 514, -62, 100, 10 },
      { -62, -53, 28, -44, 442, -116, 10, -62, 1, 10, 28, -62, -62, 514, 10, 100, -44, -53, -116, 136, -44, 514, -116, -44, -71, 19, -62, -53, -134, -44, 1, 10, -8, 64, 1, -8, 10, 100, -80, 640, 10, -80, 505, 514, -71, 568, -62, 496, -80, 496, -8, -80, -224, 496, -8, -80, 64, -512, -8, 64, 64, 640, -512, 4096, 64, -512, 496, 568, -80, 640, -8, 64, 10, -62, 1, 10, 28, -62, 1, 10, -8, 64, 1, -8, -8, -80, 64, -512, -8, 64, -62, -71, 10, -80, 1, -8, 442, -116, -134, 505, 442, -134, -134, -44, -62, 496, 28, -224, -62, -71, 10, -80, 1, -8, -53, 19, -62, 496, 10, -80 },
      { 505, -134, -134, 442, 442, -116, -71, -62, 1, -8, 10, -80, 19, -53, 10, -80, -62, 496, -134, -44, -62, 496, 28, -224, 514, -44, -44, -116, -116, 136, -62, 28, 10, 1, -62, 10, -53, -44, 100, 10, 514, -62, 442, -116, -53, -62, -44, 28, 568, 496, -8, 64, -80, 640, 496, -224, -80, -8, 496, -80, 64, -512, 640, 64, 4096, -512, -8, 64, -80, -8, -512, 64, -71, -62, 1, -8, 10, -80, -62, 28, 10, 1, -62, 10, -8, 64, -80, -8, -512, 64, 1, -8, 10, 1, 64, -8, 505, 514, -71, 568, -62, 496, -134, -44, 19, -71, -53, -62, 1, -8, 10, 1, 64, -8, 10, -80, 100, 10, 640, -80 },
      { 514, -44, -44, -116, -116, 136, -62, 28, 10, 1, -62, 10, -53, -44, 100, 10, 514, -62, 442, -116, -53, -62, -44, 28, 505, -134, -134, 442, 442, -116, -71, -62, 1, -8, 10, -80, 19, -53, 10, -80, -62, 496, -134, -44, -62, 496, 28, -224, 496, -224, -80, -8, 496, -80, 568, 496, -8, 64, -80, 640, 640, 64, 64, -512, -512, 4096, -80, -8, -8, 64, 64, -512, -62, 28, 10, 1, -62, 10, -71, -62, 1, -8, 10, -80, -80, -8, -8, 64, 64, -512, 10, 1, 1, -8, -8, 64, -134, -44, 19, -71, -53, -62, 505, 514, -71, 568, -62, 496, 10, 1, 1, -8, -8, 64, 100, 10, 10, -80, -80, 640 },
      { -71, 568, 505, 514, 496, -62, 1, -8, 10, 1, 64, -8, -134, -44, 19, -71, -53, -62, -80, 10, 640, -80, 100, 10, -62, 496, -134, -44, -224, 28, 10, -80, -62, -71, -8, 1, 442, -116, -134, 505, 442, -134, 496, -62, -80, 10, -53, 19, -8, 64, -80, -8, -512, 64, -80, 640, 496, 568, 64, -8, 496, -80, -224, 496, -8, -80, 4096, -512, -512, 64, 64, 640, 442, -116, -53, -62, -44, 28, -116, 136, -44, 514, -116, -44, -62, 10, 28, -62, 1, 10, 514, -62, -44, -53, 10, 100, -8, 1, 64, -8, 10, 1, -80, 10, -8, 1, -62, -71, -512, 64, 64, -8, -8, -80, -62, 10, 28, -62, 1, 10 },
      { 19, -71, -134, -44, -62, -53, 10, 1, 1, -8, -8, 64, 505, 514, -71, 568, -62, 496, 10, 100, -80, 640, 10, -80, -53, -62, 442, -116, 28, -44, -62, 10, 28, -62, 1, 10, -116, 136, -44, 514, -116, -44, -62, 514, 10, 100, -44, -53, -80, -8, -8, 64, 64, -512, 496, -80, -224, 496, -8, -80, -80, 640, 496, 568, 64, -8, -512, 4096, 64, 640, -512, 64, -134, -44, -62, 496, 28, -224, 442, -116, -134, 505, 442, -134, 10, -80, -62, -71, -8, 1, -62, 496, -53, 19, -80, 10, 1, 10, -8, 64, 1, -8, 10, -62, 1, 10, 28, -62, 64, -512, -8, -80, 64, -8, 10, -80, -62, -71, -8, 1 },
      { -62, 496, -134, -44, -224, 28, 10, -80, -62, -71, -8, 1, 442, -116, -134, 505, 442, -134, 496, -62, -80, 10, -53, 19, -71, 568, 505, 514, 496, -62, 1, -8, 10, 1, 64, -8, -134, -44, 19, -71, -53, -62, -80, 10, 640, -80, 100, 10, -80, 640, 496, 568, 64, -8, -8, 64, -80, -8, -512, 64, -224, 496, 496, -80, -80, -8, -512, 64, 4096, -512, 640, 64, -116, 136, -44, 514, -116, -44, 442, -116, -53, -62, -44, 28, 28, -62, -62, 10, 10, 1, -44, -53, 514, -62, 100, 10, -80, 10, -8, 1, -62, -71, -8, 1, 64, -8, 10, 1, 64, -8, -512, 64, -80, -8, 28, -62, -62, 10, 10, 1 },
      { -53, -62, 442, -116, 28, -44, -62, 10, 28, -62, 1, 10, -116, 136, -44, 514, -116, -44, -62, 514, 10, 100, -44, -53, 19, -71, -134, -44, -62, -53, 10, 1, 1, -8, -8, 64, 505, 514, -71, 568, -62, 496, 10, 100, -80, 640, 10, -80, 496, -80, -224, 496, -8, -80, -80, -8, -8, 64, 64, -512, 496, 568, -80, 640, -8, 64, 64, 640, -512, 4096, 64, -512, 442, -116, -134, 505, 442, -134, -134, -44, -62, 496, 28, -224, -62, -71, 10, -80, 1, -8, -53, 19, -62, 496, 10, -80, 10, -62, 1, 10, 28, -62, 1, 10, -8, 64, 1, -8, -8, -80, 64, -512, -8, 64, -62, -71, 10, -80, 1, -8 },
      { -134, 505, 442, -116, -134, 442, -62, -71, 10, -80, 1, -8, -134, -44, -62, 496, 28, -224, 19, -53, 10, -80, -62, 496, -44, 514, -116, 136, -44, -116, 28, -62, -62, 10, 10, 1, 442, -116, -53, -62, -44, 28, -53, -44, 100, 10, 514, -62, 496, 568, -80, 640, -8, 64, -224, 496, 496, -80, -80, -8, -8, 64, -80, -8, -512, 64, 64, -512, 640, 64, 4096, -512, 505, 514, -71, 568, -62, 496, -134, -44, 19, -71, -53, -62, 1, -8, 10, 1, 64, -8, 10, -80, 100, 10, 640, -80, -71, -62, 1, -8, 10, -80, -62, 28, 10, 1, -62, 10, -8, 64, -80, -8, -512, 64, 1, -8, 10, 1, 64, -8 },
      { -44, 514, -116, 136, -44, -116, 28, -62, -62, 10, 10, 1, 442, -116, -53, -62, -44, 28, -53, -44, 100, 10, 514, -62, -134, 505, 442, -116, -134, 442, -62, -71, 10, -80, 1, -8, -134, -44, -62, 496, 28, -224, 19, -53, 10, -80, -62, 496, -224, 496, 496, -80, -80, -8, 496, 568, -80, 640, -8, 64, -80, -8, -8, 64, 64, -512, 640, 64, 64, -512, -512, 4096, -134, -44, 19, -71, -53, -62, 505, 514, -71, 568, -62, 496, 10, 1, 1, -8, -8, 64, 100, 10, 10, -80, -80, 640, -62, 28, 10, 1, -62, 10, -71, -62, 1, -8, 10, -80, -80, -8, -8, 64, 64, -512, 10, 1, 1, -8, -8, 64 },
      { -8, 1, 64, -8, 10, 1, -80, 10, -8, 1, -62, -71, -512, 64, 64, -8, -8, -80, -62, 10, 28, -62, 1, 10, -80, 10, 640, -80, 100, 10, 496, -62, 568, -71, 514, 505, 64, -8, -8, 1, 1, 10, -53, -62, -44, -134, -71, 19, 496, -62, -80, 10, -53, 19, -224, 28, 496, -62, -44, -134, -8, 1, -80, 10, -71, -62, 442, -134, -116, 442, 505, -134, 4096, -512, -512, 64, 64, 640, -512, 64, 64, -8, -8, -80, 64, -8, 640, -80, 568, 496, -8, -80, -80, 496, 496, -224, 514, -62, -44, -53, 10, 100, -44, 28, -116, 442, -62, -53, -116, -44, 136, -116, 514, -44, 1, 10, 10, -62, -62, 28 },
      { 1, 10, -8, 64, 1, -8, 10, -62, 1, 10, 28, -62, 64, -512, -8, -80, 64, -8, 10, -80, -62, -71, -8, 1, 10, 100, -80, 640, 10, -80, -62, -53, -71, 19, -44, -134, -8, 64, 1, 10, -8, 1, -62, 496, 514, 505, 568, -71, -62, 514, 10, 100, -44, -53, 28, -44, -62, -53, -116, 442, 1, 10, 10, -62, -62, 28, -116, -44, 136, -116, 514, -44, -512, 4096, 64, 640, -512, 64, 64, -512, -8, -80, 64, -8, -8, -80, -80, 496, 496, -224, 64, -8, 640, -80, 568, 496, -62, 496, -53, 19, -80, 10, 28, -224, -44, -134, 496, -62, 442, -134, -116, 442, 505, -134, -8, 1, -80, 10, -71, -62 },
      { -80, 10, -8, 1, -62, -71, -8, 1, 64, -8, 10, 1, 64, -8, -512, 64, -80, -8, 28, -62, -62, 10, 10, 1, 496, -62, -80, 10, -53, 19, -224, 28, 496, -62, -44, -134, -8, 1, -80, 10, -71, -62, 442, -134, -116, 442, 505, -134, -80, 10, 640, -80, 100, 10, 496, -62, 568, -71, 514, 505, 64, -8, -8, 1, 1, 10, -53, -62, -44, -134, -71, 19, -512, 64, 4096, -512, 640, 64, 64, -8, 640, -80, 568, 496, -512, 64, 64, -8, -8, -80, -80, -8, 496, -224, -80, 496, -44, -53, 514, -62, 100, 10, -116, -44, 136, -116, 514, -44, -44, 28, -116, 442, -62, -53, 10, 1, -62, 28, 10, -62 },
      { 10, -62, 1, 10, 28, -62, 1, 10, -8, 64, 1, -8, -8, -80, 64, -512, -8, 64, -62, -71, 10, -80, 1, -8, -62, 514, 10, 100, -44, -53, 28, -44, -62, -53, -116, 442, 1, 10, 10, -62, -62, 28, -116, -44, 136, -116, 514, -44, 10, 100, -80, 640, 10, -80, -62, -53, -71, 19, -44, -134, -8, 64, 1, 10, -8, 1, -62, 496, 514, 505, 568, -71, 64, 640, -512, 4096, 64, -512, -8, -80, -80, 496, 496, -224, 64, -512, -8, -80, 64, -8, -8, 64, 568, 496, 640, -80, -53, 19, -62, 496, 10, -80, 442, -134, -116, 442, 505, -134, 28, -224, -44, -134, 496, -62, 1, -8, -71, -62, -80, 10 },
      { -71, -62, 1, -8, 10, -80, -62, 28, 10, 1, -62, 10, -8, 64, -80, -8, -512, 64, 1, -8, 10, 1, 64, -8, 19, -53, 10, -80, -62, 496, -134, 442, 505, -134, -116, 442, 1, -8, -71, -62, -80, 10, 28, -224, -44, -134, 496, -62, -53, -44, 100, 10, 514, -62, -44, -116, 514, -44, 136, -116, 10, 1, -62, 28, 10, -62, -44, 28, -116, 442, -62, -53, 64, -512, 640, 64, 4096, -512, -8, 64, 568, 496, 640, -80, -80, -8, 496, -224, -80, 496, -512, 64, 64, -8, -8, -80, 10, -80, 100, 10, 640, -80, -62, 496, 514, 505, 568, -71, -53, -62, -44, -134, -71, 19, 64, -8, -8, 1, 1, 10 },
      { -62, 28, 10, 1, -62, 10, -71, -62, 1, -8, 10, -80, -80, -8, -8, 64, 64, -512, 10, 1, 1, -8, -8, 64, -53, -44, 100, 10, 514, -62, -44, -116, 514, -44, 136, -116, 10, 1, -62, 28, 10, -62, -44, 28, -116, 442, -62, -53, 19, -53, 10, -80, -62, 496, -134, 442, 505, -134, -116, 442, 1, -8, -71, -62, -80, 10, 28, -224, -44, -134, 496, -62, 640, 64, 64, -512, -512, 4096, -80, -8, 496, -224, -80, 496, -8, 64, 568, 496, 640, -80, 64, -512, -8, -80, 64, -8, 100, 10, 10, -80, -80, 640, -53, -62, -44, -134, -71, 19, -62, 496, 514, 505, 568, -71, -8, 64, 1, 10, -8, 1 },
      { -80, 10, 640, -80, 100, 10, 496, -62, 568, -71, 514, 505, 64, -8, -8, 1, 1, 10, -53, -62, -44, -134, -71, 19, -8, 1, 64, -8, 10, 1, -80, 10, -8, 1, -62, -71, -512, 64, 64, -8, -8, -80, -62, 10, 28, -62, 1, 10, -224, 28, 496, -62, -44, -134, 496, -62, -80, 10, -53, 19, -80, 10, -8, 1, -62, -71, -116, 442, 442, -134, -134, 505, -512, 64, 64, -8, -8, -80, 4096, -512, -512, 64, 64, 640, 640, -80, 64, -8, 496, 568, -80, 496, -8, -80, -224, 496, -44, 28, -116, 442, -62, -53, 514, -62, -44, -53, 10, 100, 136, -116, -116, -44, -44, 514, 10, -62, 1, 10, 28, -62 },
      { 10, 100, -80, 640, 10, -80, -62, -53, -71, 19, -44, -134, -8, 64, 1, 10, -8, 1, -62, 496, 514, 505, 568, -71, 1, 10, -8, 64, 1, -8, 10, -62, 1, 10, 28, -62, 64, -512, -8, -80, 64, -8, 10, -80, -62, -71, -8, 1, 28, -44, -62, -53, -116, 442, -62, 514, 10, 100, -44, -53, 10, -62, 1, 10, 28, -62, 136, -116, -116, -44, -44, 514, 64, -512, -8, -80, 64, -8, -512, 4096, 64, 640, -512, 64, -80, 496, -8, -80, -224, 496, 640, -80, 64, -8, 496, 568, 28, -224, -44, -134, 496, -62, -62, 496, -53, 19, -80, 10, -116, 442, 442, -134, -134, 505, -80, 10, -8, 1, -62, -71 },
      { 496, -62, -80, 10, -53, 19, -224, 28, 496, -62, -44, -134, -8, 1, -80, 10, -71, -62, 442, -134, -116, 442, 505, -134, -80, 10, -8, 1, -62, -71, -8, 1, 64, -8, 10, 1, 64, -8, -512, 64, -80, -8, 28, -62, -62, 10, 10, 1, 496, -62, 568, -71, 514, 505, -80, 10, 640, -80, 100, 10, -8, 1, 64, -8, 10, 1, -44, -134, -53, -62, 19, -71, 64, -8, 640, -80, 568, 496, -512, 64, 4096, -512, 640, 64, 64, -8, -512, 64, -80, -8, 496, -224, -80, -8, 496, -80, -116, -44, 136, -116, 514, -44, -44, -53, 514, -62, 100, 10, -116, 442, -44, 28, -53, -62, -62, 28, 10, 1, -62, 10 },
      { -62, 514, 10, 100, -44, -53, 28, -44, -62, -53, -116, 442, 1, 10, 10, -62, -62, 28, -116, -44, 136, -116, 514, -44, 10, -62, 1, 10, 28, -62, 1, 10, -8, 64, 1, -8, -8, -80, 64, -512, -8, 64, -62, -71, 10, -80, 1, -8, -62, -53, -71, 19, -44, -134, 10, 100, -80, 640, 10, -80, 1, 10, -8, 64, 1, -8, 514, 505, -62, 496, -71, 568, -8, -80, -80, 496, 496, -224, 64, 640, -512, 4096, 64, -512, -8, -80, 64, -512, -8, 64, 568, 496, -8, 64, -80, 640, 442, -134, -116, 442, 505, -134, -53, 19, -62, 496, 10, -80, -44, -134, 28, -224, -62, 496, -71, -62, 1, -8, 10, -80 },
      { 19, -53, 10, -80, -62, 496, -134, 442, 505, -134, -116, 442, 1, -8, -71, -62, -80, 10, 28, -224, -44, -134, 496, -62, -71, -62, 1, -8, 10, -80, -62, 28, 10, 1, -62, 10, -8, 64, -80, -8, -512, 64, 1, -8, 10, 1, 64, -8, -44, -116, 514, -44, 136, -116, -53, -44, 100, 10, 514, -62, -62, 28, 10, 1, -62, 10, -116, 442, -44, 28, -53, -62, -8, 64, 568, 496, 640, -80, 64, -512, 640, 64, 4096, -512, 496, -224, -80, -8, 496, -80, 64, -8, -512, 64, -80, -8, -62, 496, 514, 505, 568, -71, 10, -80, 100, 10, 640, -80, -44, -134, -53, -62, 19, -71, -8, 1, 64, -8, 10, 1 },
      { -53, -44, 100, 10, 514, -62, -44, -116, 514, -44, 136, -116, 10, 1, -62, 28, 10, -62, -44, 28, -116, 442, -62, -53, -62, 28, 10, 1, -62, 10, -71, -62, 1, -8, 10, -80, -80, -8, -8, 64, 64, -512, 10, 1, 1, -8, -8, 64, -134, 442, 505, -134, -116, 442, 19, -53, 10, -80, -62, 496, -71, -62, 1, -8, 10, -80, -44, -134, 28, -224, -62, 496, -80, -8, 496, -224, -80, 496, 640, 64, 64, -512, -512, 4096, 568, 496, -8, 64, -80, 640, -8, -80, 64, -512, -8, 64, -53, -62, -44, -134, -71, 19, 100, 10, 10, -80, -80, 640, 514, 505, -62, 496, -71, 568, 1, 10, -8, 64, 1, -8 },
      { 496, -62, 568, -71, 514, 505, -80, 10, 640, -80, 100, 10, -8, 1, 64, -8, 10, 1, -44, -134, -53, -62, 19, -71, -224, 28, 496, -62, -44, -134, 496, -62, -80, 10, -53, 19, -80, 10, -8, 1, -62, -71, -116, 442, 442, -134, -134, 505, -8, 1, 64, -8, 10, 1, -80, 10, -8, 1, -62, -71, -512, 64, 64, -8, -8, -80, -62, 10, 28, -62, 1, 10, 64, -8, -512, 64, -80, -8, 640, -80, 64, -8, 496, 568, 4096, -512, -512, 64, 64, 640, 496, -80, -224, 496, -8, -80, -116, 442, -44, 28, -53, -62, 136, -116, -116, -44, -44, 514, 514, -62, -44, -53, 10, 100, -62, 10, 28, -62, 1, 10 },
      { -62, -53, -71, 19, -44, -134, 10, 100, -80, 640, 10, -80, 1, 10, -8, 64, 1, -8, 514, 505, -62, 496, -71, 568, 28, -44, -62, -53, -116, 442, -62, 514, 10, 100, -44, -53, 10, -62, 1, 10, 28, -62, 136, -116, -116, -44, -44, 514, 1, 10, -8, 64, 1, -8, 10, -62, 1, 10, 28, -62, 64, -512, -8, -80, 64, -8, 10, -80, -62, -71, -8, 1, -8, -80, 64, -512, -8, 64, -80, 496, -8, -80, -224, 496, -512, 4096, 64, 640, -512, 64, -80, 640, 496, 568, 64, -8, -44, -134, 28, -224, -62, 496, -116, 442, 442, -134, -134, 505, -62, 496, -53, 19, -80, 10, 10, -80, -62, -71, -8, 1 },
      { -224, 28, 496, -62, -44, -134, 496, -62, -80, 10, -53, 19, -80, 10, -8, 1, -62, -71, -116, 442, 442, -134, -134, 505, 496, -62, 568, -71, 514, 505, -80, 10, 640, -80, 100, 10, -8, 1, 64, -8, 10, 1, -44, -134, -53, -62, 19, -71, -80, 10, -8, 1, -62, -71, -8, 1, 64, -8, 10, 1, 64, -8, -512, 64, -80, -8, 28, -62, -62, 10, 10, 1, 640, -80, 64, -8, 496, 568, 64, -8, -512, 64, -80, -8, -512, 64, 4096, -512, 640, 64, -224, 496, 496, -80, -80, -8, 136, -116, -116, -44, -44, 514, -116, 442, -44, 28, -53, -62, -44, -53, 514, -62, 100, 10, 28, -62, -62, 10, 10, 1 },
      { 28, -44, -62, -53, -116, 442, -62, 514, 10, 100, -44, -53, 10, -62, 1, 10, 28, -62, 136, -116, -116, -44, -44, 514, -62, -53, -71, 19, -44, -134, 10, 100, -80, 640, 10, -80, 1, 10, -8, 64, 1, -8, 514, 505, -62, 496, -71, 568, 10, -62, 1, 10, 28, -62, 1, 10, -8, 64, 1, -8, -8, -80, 64, -512, -8, 64, -62, -71, 10, -80, 1, -8, -80, 496, -8, -80, -224, 496, -8, -80, 64, -512, -8, 64, 64, 640, -512, 4096, 64, -512, 496, 568, -80, 640, -8, 64, -116, 442, 442, -134, -134, 505, -44, -134, 28, -224, -62, 496, -53, 19, -62, 496, 10, -80, -62, -71, 10, -80, 1, -8 },
      { -134, 442, 505, -134, -116, 442, 19, -53, 10, -80, -62, 496, -71, -62, 1, -8, 10, -80, -44, -134, 28, -224, -62, 496, -44, -116, 514, -44, 136, -116, -53, -44, 100, 10, 514, -62, -62, 28, 10, 1, -62, 10, -116, 442, -44, 28, -53, -62, -71, -62, 1, -8, 10, -80, -62, 28, 10, 1, -62, 10, -8, 64, -80, -8, -512, 64, 1, -8, 10, 1, 64, -8, 568, 496, -8, 64, -80, 640, 496, -224, -80, -8, 496, -80, 64, -512, 640, 64, 4096, -512, -8, 64, -80, -8, -512, 64, 514, 505, -62, 496, -71, 568, -44, -134, -53, -62, 19, -71, 10, -80, 100, 10, 640, -80, 1, -8, 10, 1, 64, -8 },
      { -44, -116, 514, -44, 136, -116, -53, -44, 100, 10, 514, -62, -62, 28, 10, 1, -62, 10, -116, 442, -44, 28, -53, -62, -134, 442, 505, -134, -116, 442, 19, -53, 10, -80, -62, 496, -71, -62, 1, -8, 10, -80, -44, -134, 28, -224, -62, 496, -62, 28, 10, 1, -62, 10, -71, -62, 1, -8, 10, -80, -80, -8, -8, 64, 64, -512, 10, 1, 1, -8, -8, 64, 496, -224, -80, -8, 496, -80, 568, 496, -8, 64, -80, 640, 640, 64, 64, -512, -512, 4096, -80, -8, -8, 64, 64, -512, -44, -134, -53, -62, 19, -71, 514, 505, -62, 496, -71, 568, 100, 10, 10, -80, -80, 640, 10, 1, 1, -8, -8, 64 },
      { 505, 514, -71, 568, -62, 496, -134, -44, 19, -71, -53, -62, 1, -8, 10, 1, 64, -8, 10, -80, 100, 10, 640, -80, -134, -44, -62, 496, 28, -224, 442, -116, -134, 505, 442, -134, 10, -80, -62, -71, -8, 1, -62, 496, -53, 19, -80, 10, 442, -116, -53, -62, -44, 28, -116, 136, -44, 514, -116, -44, -62, 10, 28, -62, 1, 10, 514, -62, -44, -53, 10, 100, -8, 64, -80, -8, -512, 64, -80, 640, 496, 568, 64, -8, 496, -80, -224, 496, -8, -80, 4096, -512, -512, 64, 64, 640, 1, -8, 10, 1, 64, -8, 10, -80, -62, -71, -8, 1, -62, 10, 28, -62, 1, 10, -512, 64, 64, -8, -8, -80 },
      { -134, -44, 19, -71, -53, -62, 505, 514, -71, 568, -62, 496, 10, 1, 1, -8, -8, 64, 100, 10, 10, -80, -80, 640, 442, -116, -53, -62, -44, 28, -116, 136, -44, 514, -116, -44, -62, 10, 28, -62, 1, 10, 514, -62, -44, -53, 10, 100, -134, -44, -62, 496, 28, -224, 442, -116, -134, 505, 442, -134, 10, -80, -62, -71, -8, 1, -62, 496, -53, 19, -80, 10, -80, -8, -8, 64, 64, -512, 496, -80, -224, 496, -8, -80, -80, 640, 496, 568, 64, -8, -512, 4096, 64, 640, -512, 64, 10, 1, 1, -8, -8, 64, -62, 10, 28, -62, 1, 10, 10, -80, -62, -71, -8, 1, 64, -512, -8, -80, 64, -8 },
      { -134, -44, -62, 496, 28, -224, 442, -116, -134, 505, 442, -134, 10, -80, -62, -71, -8, 1, -62, 496, -53, 19, -80, 10, 505, 514, -71, 568, -62, 496, -134, -44, 19, -71, -53, -62, 1, -8, 10, 1, 64, -8, 10, -80, 100, 10, 640, -80, -116, 136, -44, 514, -116, -44, 442, -116, -53, -62, -44, 28, 28, -62, -62, 10, 10, 1, -44, -53, 514, -62, 100, 10, -80, 640, 496, 568, 64, -8, -8, 64, -80, -8, -512, 64, -224, 496, 496, -80, -80, -8, -512, 64, 4096, -512, 640, 64, 10, -80, -62, -71, -8, 1, 1, -8, 10, 1, 64, -8, 28, -62, -62, 10, 10, 1, 64, -8, -512, 64, -80, -8 },
      { 442, -116, -53, -62, -44, 28, -116, 136, -44, 514, -116, -44, -62, 10, 28, -62, 1, 10, 514, -62, -44, -53, 10, 100, -134, -44, 19, -71, -53, -62, 505, 514, -71, 568, -62, 496, 10, 1, 1, -8, -8, 64, 100, 10, 10, -80, -80, 640, 442, -116, -134, 505, 442, -134, -134, -44, -62, 496, 28, -224, -62, -71, 10, -80, 1, -8, -53, 19, -62, 496, 10, -80, 496, -80, -224, 496, -8, -80, -80, -8, -8, 64, 64, -512, 496, 568, -80, 640, -8, 64, 64, 640, -512, 4096, 64, -512, -62, 10, 28, -62, 1, 10, 10, 1, 1, -8, -8, 64, -62, -71, 10, -80, 1, -8, -8, -80, 64, -512, -8, 64 },
      { 442, -116, -134, 505, 442, -134, -134, -44, -62, 496, 28, -224, -62, -71, 10, -80, 1, -8, -53, 19, -62, 496, 10, -80, -116, 136, -44, 514, -116, -44, 442, -116, -53, -62, -44, 28, 28, -62, -62, 10, 10, 1, -44, -53, 514, -62, 100, 10, 505, 514, -71, 568, -62, 496, -134, -44, 19, -71, -53, -62, 1, -8, 10, 1, 64, -8, 10, -80, 100, 10, 640, -80, 496, 568, -80, 640, -8, 64, -224, 496, 496, -80, -80, -8, -8, 64, -80, -8, -512, 64, 64, -512, 640, 64, 4096, -512, -62, -71, 10, -80, 1, -8, 28, -62, -62, 10, 10, 1, 1, -8, 10, 1, 64, -8, -8, 64, -80, -8, -512, 64 },
      { -116, 136, -44, 514, -116, -44, 442, -116, -53, -62, -44, 28, 28, -62, -62, 10, 10, 1, -44, -53, 514, -62, 100, 10, 442, -116, -134, 505, 442, -134, -134, -44, -62, 496, 28, -224, -62, -71, 10, -80, 1, -8, -53, 19, -62, 496, 10, -80, -134, -44, 19, -71, -53, -62, 505, 514, -71, 568, -62, 496, 10, 1, 1, -8, -8, 64, 100, 10, 10, -80, -80, 640, -224, 496, 496, -80, -80, -8, 496, 568, -80, 640, -8, 64, -80, -8, -8, 64, 64, -512, 640, 64, 64, -512, -512, 4096, 28, -62, -62, 10, 10, 1, -62, -71, 10, -80, 1, -8, 10, 1, 1, -8, -8, 64, -80, -8, -8, 64, 64, -512 },
      { 1, -8, 10, 1, 64, -8, 10, -80, -62, -71, -8, 1, -62, 10, 28, -62, 1, 10, -512, 64, 64, -8, -8, -80, 10, -80, 100, 10, 640, -80, -62, 496, 514, 505, 568, -71, -53, -62, -44, -134, -71, 19, 64, -8, -8, 1, 1, 10, -62, 496, -53, 19, -80, 10, 28, -224, -44, -134, 496, -62, 442, -134, -116, 442, 505, -134, -8, 1, -80, 10, -71, -62, 514, -62, -44, -53, 10, 100, -44, 28, -116, 442, -62, -53, -116, -44, 136, -116, 514, -44, 1, 10, 10, -62, -62, 28, 4096, -512, -512, 64, 64, 640, -512, 64, 64, -8, -8, -80, 64, -8, 640, -80, 568, 496, -8, -80, -80, 496, 496, -224 },
      { 10, 1, 1, -8, -8, 64, -62, 10, 28, -62, 1, 10, 10, -80, -62, -71, -8, 1, 64, -512, -8, -80, 64, -8, 100, 10, 10, -80, -80, 640, -53, -62, -44, -134, -71, 19, -62, 496, 514, 505, 568, -71, -8, 64, 1, 10, -8, 1, 514, -62, -44, -53, 10, 100, -44, 28, -116, 442, -62, -53, -116, -44, 136, -116, 514, -44, 1, 10, 10, -62, -62, 28, -62, 496, -53, 19, -80, 10, 28, -224, -44, -134, 496, -62, 442, -134, -116, 442, 505, -134, -8, 1, -80, 10, -71, -62, -512, 4096, 64, 640, -512, 64, 64, -512, -8, -80, 64, -8, -8, -80, -80, 496, 496, -224, 64, -8, 640, -80, 568, 496 },
      { 10, -80, -62, -71, -8, 1, 1, -8, 10, 1, 64, -8, 28, -62, -62, 10, 10, 1, 64, -8, -512, 64, -80, -8, -62, 496, -53, 19, -80, 10, 28, -224, -44, -134, 496, -62, 442, -134, -116, 442, 505, -134, -8, 1, -80, 10, -71, -62, 10, -80, 100, 10, 640, -80, -62, 496, 514, 505, 568, -71, -53, -62, -44, -134, -71, 19, 64, -8, -8, 1, 1, 10, -44, -53, 514, -62, 100, 10, -116, -44, 136, -116, 514, -44, -44, 28, -116, 442, -62, -53, 10, 1, -62, 28, 10, -62, -512, 64, 4096, -512, 640, 64, 64, -8, 640, -80, 568, 496, -512, 64, 64, -8, -8, -80, -80, -8, 496, -224, -80, 496 },
      { -62, 10, 28, -62, 1, 10, 10, 1, 1, -8, -8, 64, -62, -71, 10, -80, 1, -8, -8, -80, 64, -512, -8, 64, 514, -62, -44, -53, 10, 100, -44, 28, -116, 442, -62, -53, -116, -44, 136, -116, 514, -44, 1, 10, 10, -62, -62, 28, 100, 10, 10, -80, -80, 640, -53, -62, -44, -134, -71, 19, -62, 496, 514, 505, 568, -71, -8, 64, 1, 10, -8, 1, -53, 19, -62, 496, 10, -80, 442, -134, -116, 442, 505, -134, 28, -224, -44, -134, 496, -62, 1, -8, -71, -62, -80, 10, 64, 640, -512, 4096, 64, -512, -8, -80, -80, 496, 496, -224, 64, -512, -8, -80, 64, -8, -8, 64, 568, 496, 640, -80 },
      { -62, -71, 10, -80, 1, -8, 28, -62, -62, 10, 10, 1, 1, -8, 10, 1, 64, -8, -8, 64, -80, -8, -512, 64, -53, 19, -62, 496, 10, -80, 442, -134, -116, 442, 505, -134, 28, -224, -44, -134, 496, -62, 1, -8, -71, -62, -80, 10, -44, -53, 514, -62, 100, 10, -116, -44, 136, -116, 514, -44, -44, 28, -116, 442, -62, -53, 10, 1, -62, 28, 10, -62, 10, -80, 100, 10, 640, -80, -62, 496, 514, 505, 568, -71, -53, -62, -44, -134, -71, 19, 64, -8, -8, 1, 1, 10, 64, -512, 640, 64, 4096, -512, -8, 64, 568, 496, 640, -80, -80, -8, 496, -224, -80, 496, -512, 64, 64, -8, -8, -80 },
      { 28, -62, -62, 10, 10, 1, -62, -71, 10, -80, 1, -8, 10, 1, 1, -8, -8, 64, -80, -8, -8, 64, 64, -512, -44, -53, 514, -62, 100, 10, -116, -44, 136, -116, 514, -44, -44, 28, -116, 442, -62, -53, 10, 1, -62, 28, 10, -62, -53, 19, -62, 496, 10, -80, 442, -134, -116, 442, 505, -134, 28, -224, -44, -134, 496, -62, 1, -8, -71, -62, -80, 10, 100, 10, 10, -80, -80, 640, -53, -62, -44, -134, -71, 19, -62, 496, 514, 505, 568, -71, -8, 64, 1, 10, -8, 1, 640, 64, 64, -512, -512, 4096, -80, -8, 496, -224, -80, 496, -8, 64, 568, 496, 640, -80, 64, -512, -8, -80, 64, -8 },
      { 10, -80, 100, 10, 640, -80, -62, 496, 514, 505, 568, -71, -53, -62, -44, -134, -71, 19, 64, -8, -8, 1, 1, 10, 1, -8, 10, 1, 64, -8, 10, -80, -62, -71, -8, 1, -62, 10, 28, -62, 1, 10, -512, 64, 64, -8, -8, -80, 28, -224, -44, -134, 496, -62, -62, 496, -53, 19, -80, 10, -116, 442, 442, -134, -134, 505, -80, 10, -8, 1, -62, -71, -44, 28, -116, 442, -62, -53, 514, -62, -44, -53, 10, 100, 136, -116, -116, -44, -44, 514, 10, -62, 1, 10, 28, -62, -512, 64, 64, -8, -8, -80, 4096, -512, -512, 64, 64, 640, 640, -80, 64, -8, 496, 568, -80, 496, -8, -80, -224, 496 },
      { 100, 10, 10, -80, -80, 640, -53, -62, -44, -134, -71, 19, -62, 496, 514, 505, 568, -71, -8, 64, 1, 10, -8, 1, 10, 1, 1, -8, -8, 64, -62, 10, 28, -62, 1, 10, 10, -80, -62, -71, -8, 1, 64, -512, -8, -80, 64, -8, -44, 28, -116, 442, -62, -53, 514, -62, -44, -53, 10, 100, 136, -116, -116, -44, -44, 514, 10, -62, 1, 10, 28, -62, 28, -224, -44, -134, 496, -62, -62, 496, -53, 19, -80, 10, -116, 442, 442, -134, -134, 505, -80, 10, -8, 1, -62, -71, 64, -512, -8, -80, 64, -8, -512, 4096, 64, 640, -512, 64, -80, 496, -8, -80, -224, 496, 640, -80, 64, -8, 496, 568 },
      { -62, 496, -53, 19, -80, 10, 28, -224, -44, -134, 496, -62, 442, -134, -116, 442, 505, -134, -8, 1, -80, 10, -71, -62, 10, -80, -62, -71, -8, 1, 1, -8, 10, 1, 64, -8, 28, -62, -62, 10, 10, 1, 64, -8, -512, 64, -80, -8, -62, 496, 514, 505, 568, -71, 10, -80, 100, 10, 640, -80, -44, -134, -53, -62, 19, -71, -8, 1, 64, -8, 10, 1, -116, -44, 136, -116, 514, -44, -44, -53, 514, -62, 100, 10, -116, 442, -44, 28, -53, -62, -62, 28, 10, 1, -62, 10, 64, -8, 640, -80, 568, 496, -512, 64, 4096, -512, 640, 64, 64, -8, -512, 64, -80, -8, 496, -224, -80, -8, 496, -80 },
      { 514, -62, -44, -53, 10, 100, -44, 28, -116, 442, -62, -53, -116, -44, 136, -116, 514, -44, 1, 10, 10, -62, -62, 28, -62, 10, 28, -62, 1, 10, 10, 1, 1, -8, -8, 64, -62, -71, 10, -80, 1, -8, -8, -80, 64, -512, -8, 64, -53, -62, -44, -134, -71, 19, 100, 10, 10, -80, -80, 640, 514, 505, -62, 496, -71, 568, 1, 10, -8, 64, 1, -8, 442, -134, -116, 442, 505, -134, -53, 19, -62, 496, 10, -80, -44, -134, 28, -224, -62, 496, -71, -62, 1, -8, 10, -80, -8, -80, -80, 496, 496, -224, 64, 640, -512, 4096, 64, -512, -8, -80, 64, -512, -8, 64, 568, 496, -8, 64, -80, 640 },
      { -53, 19, -62, 496, 10, -80, 442, -134, -116, 442, 505, -134, 28, -224, -44, -134, 496, -62, 1, -8, -71, -62, -80, 10, -62, -71, 10, -80, 1, -8, 28, -62, -62, 10, 10, 1, 1, -8, 10, 1, 64, -8, -8, 64, -80, -8, -512, 64, -116, -44, 136, -116, 514, -44, -44, -53, 514, -62, 100, 10, -116, 442, -44, 28, -53, -62, -62, 28, 10, 1, -62, 10, -62, 496, 514, 505, 568, -71, 10, -80, 100, 10, 640, -80, -44, -134, -53, -62, 19, -71, -8, 1, 64, -8, 10, 1, -8, 64, 568, 496, 640, -80, 64, -512, 640, 64, 4096, -512, 496, -224, -80, -8, 496, -80, 64, -8, -512, 64, -80, -8 },
      { -44, -53, 514, -62, 100, 10, -116, -44, 136, -116, 514, -44, -44, 28, -116, 442, -62, -53, 10, 1, -62, 28, 10, -62, 28, -62, -62, 10, 10, 1, -62, -71, 10, -80, 1, -8, 10, 1, 1, -8, -8, 64, -80, -8, -8, 64, 64, -512, 442, -134, -116, 442, 505, -134, -53, 19, -62, 496, 10, -80, -44, -134, 28, -224, -62, 496, -71, -62, 1, -8, 10, -80, -53, -62, -44, -134, -71, 19, 100, 10, 10, -80, -80, 640, 514, 505, -62, 496, -71, 568, 1, 10, -8, 64, 1, -8, -80, -8, 496, -224, -80, 496, 640, 64, 64, -512, -512, 4096, 568, 496, -8, 64, -80, 640, -8, -80, 64, -512, -8, 64 },
      { -62, 496, 514, 505, 568, -71, 10, -80, 100, 10, 640, -80, -44, -134, -53, -62, 19, -71, -8, 1, 64, -8, 10, 1, 28, -224, -44, -134, 496, -62, -62, 496, -53, 19, -80, 10, -116, 442, 442, -134, -134, 505, -80, 10, -8, 1, -62, -71, 1, -8, 10, 1, 64, -8, 10, -80, -62, -71, -8, 1, -62, 10, 28, -62, 1, 10, -512, 64, 64, -8, -8, -80, -116, 442, -44, 28, -53, -62, 136, -116, -116, -44, -44, 514, 514, -62, -44, -53, 10, 100, -62, 10, 28, -62, 1, 10, 64, -8, -512, 64, -80, -8, 640, -80, 64, -8, 496, 568, 4096, -512, -512, 64, 64, 640, 496, -80, -224, 496, -8, -80 },
      { -53, -62, -44, -134, -71, 19, 100, 10, 10, -80, -80, 640, 514, 505, -62, 496, -71, 568, 1, 10, -8, 64, 1, -8, -44, 28, -116, 442, -62, -53, 514, -62, -44, -53, 10, 100, 136, -116, -116, -44, -44, 514, 10, -62, 1, 10, 28, -62, 10, 1, 1, -8, -8, 64, -62, 10, 28, -62, 1, 10, 10, -80, -62, -71, -8, 1, 64, -512, -8, -80, 64, -8, -44, -134, 28, -224, -62, 496, -116, 442, 442, -134, -134, 505, -62, 496, -53, 19, -80, 10, 10, -80, -62, -71, -8, 1, -8, -80, 64, -512, -8, 64, -80, 496, -8, -80, -224, 496, -512, 4096, 64, 640, -512, 64, -80, 640, 496, 568, 64, -8 },
      { 28, -224, -44, -134, 496, -62, -62, 496, -53, 19, -80, 10, -116, 442, 442, -134, -134, 505, -80, 10, -8, 1, -62, -71, -62, 496, 514, 505, 568, -71, 10, -80, 100, 10, 640, -80, -44, -134, -53, -62, 19, -71, -8, 1, 64, -8, 10, 1, 10, -80, -62, -71, -8, 1, 1, -8, 10, 1, 64, -8, 28, -62, -62, 10, 10, 1, 64, -8, -512, 64, -80, -8, 136, -116, -116, -44, -44, 514, -116, 442, -44, 28, -53, -62, -44, -53, 514, -62, 100, 10, 28, -62, -62, 10, 10, 1, 640, -80, 64, -8, 496, 568, 64, -8, -512, 64, -80, -8, -512, 64, 4096, -512, 640, 64, -224, 496, 496, -80, -80, -8 },
      { -44, 28, -116, 442, -62, -53, 514, -62, -44, -53, 10, 100, 136, -116, -116, -44, -44, 514, 10, -62, 1, 10, 28, -62, -53, -62, -44, -134, -71, 19, 100, 10, 10, -80, -80, 640, 514, 505, -62, 496, -71, 568, 1, 10, -8, 64, 1, -8, -62, 10, 28, -62, 1, 10, 10, 1, 1, -8, -8, 64, -62, -71, 10, -80, 1, -8, -8, -80, 64, -512, -8, 64, -116, 442, 442, -134, -134, 505, -44, -134, 28, -224, -62, 496, -53, 19, -62, 496, 10, -80, -62, -71, 10, -80, 1, -8, -80, 496, -8, -80, -224, 496, -8, -80, 64, -512, -8, 64, 64, 640, -512, 4096, 64, -512, 496, 568, -80, 640, -8, 64 },
      { 442, -134, -116, 442, 505, -134, -53, 19, -62, 496, 10, -80, -44, -134, 28, -224, -62, 496, -71, -62, 1, -8, 10, -80, -116, -44, 136, -116, 514, -44, -44, -53, 514, -62, 100, 10, -116, 442, -44, 28, -53, -62, -62, 28, 10, 1, -62, 10, -62, -71, 10, -80, 1, -8, 28, -62, -62, 10, 10, 1, 1, -8, 10, 1, 64, -8, -8, 64, -80, -8, -512, 64, 514, 505, -62, 496, -71, 568, -44, -134, -53, -62, 19, -71, 10, -80, 100, 10, 640, -80, 1, -8, 10, 1, 64, -8, 568, 496, -8, 64, -80, 640, 496, -224, -80, -8, 496, -80, 64, -512, 640, 64, 4096, -512, -8, 64, -80, -8, -512, 64 },
      { -116, -44, 136, -116, 514, -44, -44, -53, 514, -62, 100, 10, -116, 442, -44, 28, -53, -62, -62, 28, 10, 1, -62, 10, 442, -134, -116, 442, 505, -134, -53, 19, -62, 496, 10, -80, -44, -134, 28, -224, -62, 496, -71, -62, 1, -8, 10, -80, 28, -62, -62, 10, 10, 1, -62, -71, 10, -80, 1, -8, 10, 1, 1, -8, -8, 64, -80, -8, -8, 64, 64, -512, -44, -134, -53, -62, 19, -71, 514, 505, -62, 496, -71, 568, 100, 10, 10, -80, -80, 640, 10, 1, 1, -8, -8, 64, 496, -224, -80, -8, 496, -80, 568, 496, -8, 64, -80, 640, 640, 64, 64, -512, -512, 4096, -80, -8, -8, 64, 64, -512 },
      { 514, 505, -62, 496, -71, 568, -44, -134, -53, -62, 19, -71, 10, -80, 100, 10, 640, -80, 1, -8, 10, 1, 64, -8, -44, -134, 28, -224, -62, 496, -116, 442, 442, -134, -134, 505, -62, 496, -53, 19, -80, 10, 10, -80, -62, -71, -8, 1, -116, 442, -44, 28, -53, -62, 136, -116, -116, -44, -44, 514, 514, -62, -44, -53, 10, 100, -62, 10, 28, -62, 1, 10, 1, -8, 10, 1, 64, -8, 10, -80, -62, -71, -8, 1, -62, 10, 28, -62, 1, 10, -512, 64, 64, -8, -8, -80, -8, 64, -80, -8, -512, 64, -80, 640, 496, 568, 64, -8, 496, -80, -224, 496, -8, -80, 4096, -512, -512, 64, 64, 640 },
      { -44, -134, -53, -62, 19, -71, 514, 505, -62, 496, -71, 568, 100, 10, 10, -80, -80, 640, 10, 1, 1, -8, -8, 64, -116, 442, -44, 28, -53, -62, 136, -116, -116, -44, -44, 514, 514, -62, -44, -53, 10, 100, -62, 10, 28, -62, 1, 10, -44, -134, 28, -224, -62, 496, -116, 442, 442, -134, -134, 505, -62, 496, -53, 19, -80, 10, 10, -80, -62, -71, -8, 1, 10, 1, 1, -8, -8, 64, -62, 10, 28, -62, 1, 10, 10, -80, -62, -71, -8, 1, 64, -512, -8, -80, 64, -8, -80, -8, -8, 64, 64, -512, 496, -80, -224, 496, -8, -80, -80, 640, 496, 568, 64, -8, -512, 4096, 64, 640, -512, 64 },
      { -44, -134, 28, -224, -62, 496, -116, 442, 442, -134, -134, 505, -62, 496, -53, 19, -80, 10, 10, -80, -62, -71, -8, 1, 514, 505, -62, 496, -71, 568, -44, -134, -53, -62, 19, -71, 10, -80, 100, 10, 640, -80, 1, -8, 10, 1, 64, -8, 136, -116, -116, -44, -44, 514, -116, 442, -44, 28, -53, -62, -44, -53, 514, -62, 100, 10, 28, -62, -62, 10, 10, 1, 10, -80, -62, -71, -8, 1, 1, -8, 10, 1, 64, -8, 28, -62, -62, 10, 10, 1, 64, -8, -512, 64, -80, -8, -80, 640, 496, 568, 64, -8, -8, 64, -80, -8, -512, 64, -224, 496, 496, -80, -80, -8, -512, 64, 4096, -512, 640, 64 },
      { -116, 442, -44, 28, -53, -62, 136, -116, -116, -44, -44, 514, 514, -62, -44, -53, 10, 100, -62, 10, 28, -62, 1, 10, -44, -134, -53, -62, 19, -71, 514, 505, -62, 496, -71, 568, 100, 10, 10, -80, -80, 640, 10, 1, 1, -8, -8, 64, -116, 442, 442, -134, -134, 505, -44, -134, 28, -224, -62, 496, -53, 19, -62, 496, 10, -80, -62, -71, 10, -80, 1, -8, -62, 10, 28, -62, 1, 10, 10, 1, 1, -8, -8, 64, -62, -71, 10, -80, 1, -8, -8, -80, 64, -512, -8, 64, 496, -80, -224, 496, -8, -80, -80, -8, -8, 64, 64, -512, 496, 568, -80, 640, -8, 64, 64, 640, -512, 4096, 64, -512 },
      { -116, 442, 442, -134, -134, 505, -44, -134, 28, -224, -62, 496, -53, 19, -62, 496, 10, -80, -62, -71, 10, -80, 1, -8, 136, -116, -116, -44, -44, 514, -116, 442, -44, 28, -53, -62, -44, -53, 514, -62, 100, 10, 28, -62, -62, 10, 10, 1, 514, 505, -62, 496, -71, 568, -44, -134, -53, -62, 19, -71, 10, -80, 100, 10, 640, -80, 1, -8, 10, 1, 64, -8, -62, -71, 10, -80, 1, -8, 28, -62, -62, 10, 10, 1, 1, -8, 10, 1, 64, -8, -8, 64, -80, -8, -512, 64, 496, 568, -80, 640, -8, 64, -224, 496, 496, -80, -80, -8, -8, 64, -80, -8, -512, 64, 64, -512, 640, 64, 4096, -512 },
      { 136, -116, -116, -44, -44, 514, -116, 442, -44, 28, -53, -62, -44, -53, 514, -62, 100, 10, 28, -62, -62, 10, 10, 1, -116, 442, 442, -134, -134, 505, -44, -134, 28, -224, -62, 496, -53, 19, -62, 496, 10, -80, -62, -71, 10, -80, 1, -8, -44, -134, -53, -62, 19, -71, 514, 505, -62, 496, -71, 568, 100, 10, 10, -80, -80, 640, 10, 1, 1, -8, -8, 64, 28, -62, -62, 10, 10, 1, -62, -71, 10, -80, 1, -8, 10, 1, 1, -8, -8, 64, -80, -8, -8, 64, 64, -512, -224, 496, 496, -80, -80, -8, 496, 568, -80, 640, -8, 64, -80, -8, -8, 64, 64, -512, 640, 64, 64, -512, -512, 4096 } }; // 2-D array[120][120]

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

    // Use the property that M is a real matrix (see #475):
    // we can rewrite the quadratic form (A-iB)(M)(A+iB) as AMA - iBMA + iBMA + BMB = AMA + BMB
    // In addition, on C++ use the property that M is symmetric (see #475),
    // and also use constexpr to compute "2*" and "/denom[icol]" once and for all at compile time:
    // we gain (not a factor 2...) in speed here as we only loop over the up diagonal part of the matrix.
    // Strangely, CUDA is slower instead, so keep the old implementation for the moment.

#ifndef MGONGPUCPP_GPUIMPL

    // === C++ START ===
    fptype_sv deltaMEs = { 0 };
    fptype2_sv jampR_sv[ncolor];
    fptype2_sv jampI_sv[ncolor];
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
    // Mixed mode: merge two neppV vectors into one neppV2 vector
    fptype_sv deltaMEs_next = { 0 };
    for( int icol = 0; icol < ncolor; icol++ )
    {
      jampR_sv[icol] = fpvmerge( cxreal( jamp_sv[icol] ), cxreal( jamp_sv[ncolor + icol] ) );
      jampI_sv[icol] = fpvmerge( cximag( jamp_sv[icol] ), cximag( jamp_sv[ncolor + icol] ) );
    }
#else
    // Double/Float mode: one neppV vector is one neppV2 vector
    for( int icol = 0; icol < ncolor; icol++ )
    {
      jampR_sv[icol] = (fptype2_sv)( cxreal( jamp_sv[icol] ) );
      jampI_sv[icol] = (fptype2_sv)( cximag( jamp_sv[icol] ) );
    }
#endif
    // Loop over icol
    for( int icol = 0; icol < ncolor; icol++ )
    {
      // Diagonal terms
      fptype2_sv& jampRi_sv = jampR_sv[icol];
      fptype2_sv& jampIi_sv = jampI_sv[icol];
      fptype2_sv ztempR_sv = cf2.value[icol][icol] * jampRi_sv;
      fptype2_sv ztempI_sv = cf2.value[icol][icol] * jampIi_sv;
      // Loop over jcol
      for( int jcol = icol + 1; jcol < ncolor; jcol++ )
      {
        // Off-diagonal terms
        fptype2_sv& jampRj_sv = jampR_sv[jcol];
        fptype2_sv& jampIj_sv = jampI_sv[jcol];
        ztempR_sv += cf2.value[icol][jcol] * jampRj_sv;
        ztempI_sv += cf2.value[icol][jcol] * jampIj_sv;
      }
      fptype2_sv deltaMEs2 = ( jampRi_sv * ztempR_sv + jampIi_sv * ztempI_sv ); // may underflow #831
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
      deltaMEs += fpvsplit0( deltaMEs2 );
      deltaMEs_next += fpvsplit1( deltaMEs2 );
#else
      deltaMEs += deltaMEs2;
#endif
    }
    // === C++ END ===

#else

    // === CUDA START ===
    fptype_sv deltaMEs = { 0 };
    using J_ACCESS = DeviceAccessJamp;
    cxtype jamp_sv[ncolor];
    for( int icol = 0; icol < ncolor; icol++ )
      jamp_sv[icol] = J_ACCESS::kernelAccessIcolConst( allJamps, icol );
    // Loop over icol
    for( int icol = 0; icol < ncolor; icol++ )
    {
      fptype2_sv ztempR_sv = { 0 };
      fptype2_sv ztempI_sv = { 0 };
      // Loop over jcol
      for( int jcol = 0; jcol < ncolor; jcol++ )
      {
        fptype2_sv jampRj_sv = cxreal( jamp_sv[jcol] );
        fptype2_sv jampIj_sv = cximag( jamp_sv[jcol] );
        ztempR_sv += cf[icol][jcol] * jampRj_sv;
        ztempI_sv += cf[icol][jcol] * jampIj_sv;
      }
      deltaMEs += ( ztempR_sv * cxreal( jamp_sv[icol] ) + ztempI_sv * cximag( jamp_sv[icol] ) ) / denom[icol];
    }
    // === CUDA END ===

#endif

    // *** STORE THE RESULTS ***
#ifndef MGONGPUCPP_GPUIMPL
    fptype* MEs = E_ACCESS::ieventAccessRecord( allMEs, ievt00 );
#else
    fptype* MEs = allMEs;
#endif
    // NB: color_sum ADDS |M|^2 for one helicity to the running sum of |M|^2 over helicities for the given event(s)
    fptype_sv& MEs_sv = E_ACCESS::kernelAccess( MEs );
    MEs_sv += deltaMEs; // fix #435
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
    fptype* MEs_next = E_ACCESS::ieventAccessRecord( allMEs, ievt00 + neppV );
    fptype_sv& MEs_sv_next = E_ACCESS::kernelAccess( MEs_next );
    MEs_sv_next += deltaMEs_next;
#endif
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
      { -1, -1, -1, 1, -1, -1, -1 },
      { -1, -1, -1, 1, -1, -1, 1 },
      { -1, -1, -1, 1, -1, 1, -1 },
      { -1, -1, -1, 1, -1, 1, 1 },
      { -1, -1, -1, 1, 1, -1, -1 },
      { -1, -1, -1, 1, 1, -1, 1 },
      { -1, -1, -1, 1, 1, 1, -1 },
      { -1, -1, -1, 1, 1, 1, 1 },
      { -1, -1, -1, -1, -1, -1, -1 },
      { -1, -1, -1, -1, -1, -1, 1 },
      { -1, -1, -1, -1, -1, 1, -1 },
      { -1, -1, -1, -1, -1, 1, 1 },
      { -1, -1, -1, -1, 1, -1, -1 },
      { -1, -1, -1, -1, 1, -1, 1 },
      { -1, -1, -1, -1, 1, 1, -1 },
      { -1, -1, -1, -1, 1, 1, 1 },
      { -1, -1, 1, 1, -1, -1, -1 },
      { -1, -1, 1, 1, -1, -1, 1 },
      { -1, -1, 1, 1, -1, 1, -1 },
      { -1, -1, 1, 1, -1, 1, 1 },
      { -1, -1, 1, 1, 1, -1, -1 },
      { -1, -1, 1, 1, 1, -1, 1 },
      { -1, -1, 1, 1, 1, 1, -1 },
      { -1, -1, 1, 1, 1, 1, 1 },
      { -1, -1, 1, -1, -1, -1, -1 },
      { -1, -1, 1, -1, -1, -1, 1 },
      { -1, -1, 1, -1, -1, 1, -1 },
      { -1, -1, 1, -1, -1, 1, 1 },
      { -1, -1, 1, -1, 1, -1, -1 },
      { -1, -1, 1, -1, 1, -1, 1 },
      { -1, -1, 1, -1, 1, 1, -1 },
      { -1, -1, 1, -1, 1, 1, 1 },
      { -1, 1, -1, 1, -1, -1, -1 },
      { -1, 1, -1, 1, -1, -1, 1 },
      { -1, 1, -1, 1, -1, 1, -1 },
      { -1, 1, -1, 1, -1, 1, 1 },
      { -1, 1, -1, 1, 1, -1, -1 },
      { -1, 1, -1, 1, 1, -1, 1 },
      { -1, 1, -1, 1, 1, 1, -1 },
      { -1, 1, -1, 1, 1, 1, 1 },
      { -1, 1, -1, -1, -1, -1, -1 },
      { -1, 1, -1, -1, -1, -1, 1 },
      { -1, 1, -1, -1, -1, 1, -1 },
      { -1, 1, -1, -1, -1, 1, 1 },
      { -1, 1, -1, -1, 1, -1, -1 },
      { -1, 1, -1, -1, 1, -1, 1 },
      { -1, 1, -1, -1, 1, 1, -1 },
      { -1, 1, -1, -1, 1, 1, 1 },
      { -1, 1, 1, 1, -1, -1, -1 },
      { -1, 1, 1, 1, -1, -1, 1 },
      { -1, 1, 1, 1, -1, 1, -1 },
      { -1, 1, 1, 1, -1, 1, 1 },
      { -1, 1, 1, 1, 1, -1, -1 },
      { -1, 1, 1, 1, 1, -1, 1 },
      { -1, 1, 1, 1, 1, 1, -1 },
      { -1, 1, 1, 1, 1, 1, 1 },
      { -1, 1, 1, -1, -1, -1, -1 },
      { -1, 1, 1, -1, -1, -1, 1 },
      { -1, 1, 1, -1, -1, 1, -1 },
      { -1, 1, 1, -1, -1, 1, 1 },
      { -1, 1, 1, -1, 1, -1, -1 },
      { -1, 1, 1, -1, 1, -1, 1 },
      { -1, 1, 1, -1, 1, 1, -1 },
      { -1, 1, 1, -1, 1, 1, 1 },
      { 1, -1, -1, 1, -1, -1, -1 },
      { 1, -1, -1, 1, -1, -1, 1 },
      { 1, -1, -1, 1, -1, 1, -1 },
      { 1, -1, -1, 1, -1, 1, 1 },
      { 1, -1, -1, 1, 1, -1, -1 },
      { 1, -1, -1, 1, 1, -1, 1 },
      { 1, -1, -1, 1, 1, 1, -1 },
      { 1, -1, -1, 1, 1, 1, 1 },
      { 1, -1, -1, -1, -1, -1, -1 },
      { 1, -1, -1, -1, -1, -1, 1 },
      { 1, -1, -1, -1, -1, 1, -1 },
      { 1, -1, -1, -1, -1, 1, 1 },
      { 1, -1, -1, -1, 1, -1, -1 },
      { 1, -1, -1, -1, 1, -1, 1 },
      { 1, -1, -1, -1, 1, 1, -1 },
      { 1, -1, -1, -1, 1, 1, 1 },
      { 1, -1, 1, 1, -1, -1, -1 },
      { 1, -1, 1, 1, -1, -1, 1 },
      { 1, -1, 1, 1, -1, 1, -1 },
      { 1, -1, 1, 1, -1, 1, 1 },
      { 1, -1, 1, 1, 1, -1, -1 },
      { 1, -1, 1, 1, 1, -1, 1 },
      { 1, -1, 1, 1, 1, 1, -1 },
      { 1, -1, 1, 1, 1, 1, 1 },
      { 1, -1, 1, -1, -1, -1, -1 },
      { 1, -1, 1, -1, -1, -1, 1 },
      { 1, -1, 1, -1, -1, 1, -1 },
      { 1, -1, 1, -1, -1, 1, 1 },
      { 1, -1, 1, -1, 1, -1, -1 },
      { 1, -1, 1, -1, 1, -1, 1 },
      { 1, -1, 1, -1, 1, 1, -1 },
      { 1, -1, 1, -1, 1, 1, 1 },
      { 1, 1, -1, 1, -1, -1, -1 },
      { 1, 1, -1, 1, -1, -1, 1 },
      { 1, 1, -1, 1, -1, 1, -1 },
      { 1, 1, -1, 1, -1, 1, 1 },
      { 1, 1, -1, 1, 1, -1, -1 },
      { 1, 1, -1, 1, 1, -1, 1 },
      { 1, 1, -1, 1, 1, 1, -1 },
      { 1, 1, -1, 1, 1, 1, 1 },
      { 1, 1, -1, -1, -1, -1, -1 },
      { 1, 1, -1, -1, -1, -1, 1 },
      { 1, 1, -1, -1, -1, 1, -1 },
      { 1, 1, -1, -1, -1, 1, 1 },
      { 1, 1, -1, -1, 1, -1, -1 },
      { 1, 1, -1, -1, 1, -1, 1 },
      { 1, 1, -1, -1, 1, 1, -1 },
      { 1, 1, -1, -1, 1, 1, 1 },
      { 1, 1, 1, 1, -1, -1, -1 },
      { 1, 1, 1, 1, -1, -1, 1 },
      { 1, 1, 1, 1, -1, 1, -1 },
      { 1, 1, 1, 1, -1, 1, 1 },
      { 1, 1, 1, 1, 1, -1, -1 },
      { 1, 1, 1, 1, 1, -1, 1 },
      { 1, 1, 1, 1, 1, 1, -1 },
      { 1, 1, 1, 1, 1, 1, 1 },
      { 1, 1, 1, -1, -1, -1, -1 },
      { 1, 1, 1, -1, -1, -1, 1 },
      { 1, 1, 1, -1, -1, 1, -1 },
      { 1, 1, 1, -1, -1, 1, 1 },
      { 1, 1, 1, -1, 1, -1, -1 },
      { 1, 1, 1, -1, 1, -1, 1 },
      { 1, 1, 1, -1, 1, 1, -1 },
      { 1, 1, 1, -1, 1, 1, 1 } };
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
    m_pars = Parameters_sm::getInstance();
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
    m_masses.push_back( m_pars->mdl_MT );
    m_masses.push_back( m_pars->mdl_MT );
    m_masses.push_back( m_pars->ZERO );
    m_masses.push_back( m_pars->ZERO );
    m_masses.push_back( m_pars->ZERO );
    // Read physics parameters like masses and couplings from user configuration files (static: initialize once)
    // Then copy them to CUDA constant memory (issue #39) or its C++ emulation in file-scope static memory
    const fptype tIPD[nIPD] = { (fptype)m_pars->mdl_MT, (fptype)m_pars->mdl_WT };
    //const cxtype tIPC[0] = { ... }; // nIPC=0
#ifdef MGONGPUCPP_GPUIMPL
    gpuMemcpyToSymbol( cIPD, tIPD, nIPD * sizeof( fptype ) );
    //gpuMemcpyToSymbol( cIPC, tIPC, 0 * sizeof( cxtype ) ); // nIPC=0
#ifdef MGONGPUCPP_NBSMINDEPPARAM_GT_0
    if( Parameters_sm::nBsmIndepParam > 0 )
      gpuMemcpyToSymbol( bsmIndepParam, m_pars->mdl_bsmIndepParam, Parameters_sm::nBsmIndepParam * sizeof( double ) );
#endif
#else
    memcpy( cIPD, tIPD, nIPD * sizeof( fptype ) );
    //memcpy( cIPC, tIPC, nIPC * sizeof( cxtype ) ); // nIPC=0
#ifdef MGONGPUCPP_NBSMINDEPPARAM_GT_0
    if( Parameters_sm::nBsmIndepParam > 0 )
      memcpy( bsmIndepParam, m_pars->mdl_bsmIndepParam, Parameters_sm::nBsmIndepParam * sizeof( double ) );
#endif
#endif
    //for ( int i=0; i<nIPD; i++ ) std::cout << std::setprecision(17) << "tIPD[i] = " << tIPD[i] << std::endl;
    //for ( int i=0; i<Parameters_sm::nBsmIndepParam; i++ ) std::cout << std::setprecision(17) << "m_pars->mdl_bsmIndepParam[i] = " << m_pars->mdl_bsmIndepParam[i] << std::endl;
  }
#else
  // Initialize process (with hardcoded parameters)
  void
  CPPProcess::initProc( const std::string& /*param_card_name*/ )
  {
    // Use hardcoded physics parameters
    if( m_verbose )
    {
      Parameters_sm::printIndependentParameters();
      Parameters_sm::printIndependentCouplings();
      //Parameters_sm::printDependentParameters(); // now computed event-by-event (running alphas #373)
      //Parameters_sm::printDependentCouplings(); // now computed event-by-event (running alphas #373)
    }
    // Set external particle masses for this matrix element
    m_masses.push_back( Parameters_sm::ZERO );
    m_masses.push_back( Parameters_sm::ZERO );
    m_masses.push_back( Parameters_sm::mdl_MT );
    m_masses.push_back( Parameters_sm::mdl_MT );
    m_masses.push_back( Parameters_sm::ZERO );
    m_masses.push_back( Parameters_sm::ZERO );
    m_masses.push_back( Parameters_sm::ZERO );
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
    using CD_ACCESS = DeviceAccessCouplings;
    G2COUP<G_ACCESS, CD_ACCESS>( allgs, allcouplings, bsmIndepParam );
#else
    using namespace mg5amcCpu;
    using G_ACCESS = HostAccessGs;
    using CD_ACCESS = HostAccessCouplings;
    for( int ipagV = 0; ipagV < nevt / neppV; ++ipagV )
    {
      const int ievt0 = ipagV * neppV;
      const fptype* gs = MemoryAccessGs::ieventAccessRecordConst( allgs, ievt0 );
      fptype* couplings = MemoryAccessCouplings::ieventAccessRecord( allcouplings, ievt0 );
      G2COUP<G_ACCESS, CD_ACCESS>( gs, couplings, bsmIndepParam );
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
                       fptype* allJamps,           // tmp: jamp[ncolor*2*nevt] _for one helicity_ (reused in the helicity loop)
                       fptype* allWfs,             // tmp: wf[nwf*nw6*2*nevt]
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
      gpuMemset( allJamps, 0, maxtry * ncolor * mgOnGpu::nx2 * sizeof( fptype ) );
      // NB: color_sum ADDS |M|^2 for one helicity to the running sum of |M|^2 over helicities for the given event(s)
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      constexpr unsigned int* allChannelIds = nullptr; // disable multichannel single-diagram enhancement
      calculate_jamps( ihel, allmomenta, allcouplings, allJamps, allWfs, allChannelIds, allNumerators, allDenominators, 0, gpublocks, gputhreads );
#else
      calculate_jamps( ihel, allmomenta, allcouplings, allJamps, allWfs, 0, gpublocks, gputhreads );
#endif
      gpuLaunchKernel( color_sum, gpublocks, gputhreads, allMEs, allJamps );
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
        //std::cout << "sigmaKin_getGoodHel ihel=" << ihel << ( isGoodHel[ihel] ? " true" : " false" ) << std::endl;
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
        cxtype_sv jamp_sv[2 * ncolor] = {}; // all zeros
#else
        cxtype_sv jamp_sv[ncolor] = {};  // all zeros
#endif
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL /* clang-format off */
        constexpr unsigned int channelId = 0; // disable multichannel single-diagram enhancement
        calculate_jamps( ihel, allmomenta, allcouplings, jamp_sv, channelId, allNumerators, allDenominators, ievt00 ); //maxtry?
#else
        calculate_jamps( ihel, allmomenta, allcouplings, jamp_sv, ievt00 ); //maxtry?
#endif /* clang-format on */
        color_sum( allMEs, jamp_sv, ievt00 );
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
                    fptype* ghelAllNumerators,         // input/tmp: allNumerators super-buffer for nGoodHel <= ncomb individual helicities (index is ighel)
                    fptype* ghelAllDenominators,       // input/tmp: allNumerators super-buffer for nGoodHel <= ncomb individual helicities (index is ighel)
                    const unsigned int* allChannelIds, // input: multichannel channelIds[nevt] (1 to #diagrams); nullptr to disable SDE enhancement (fix #899/#911)
#endif
                    const fptype globaldenom ) /* clang-format on */
  {
    const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread)
    allMEs[ievt] /= globaldenom;
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    const int nevt = gridDim.x * blockDim.x;
    if( allChannelIds != nullptr ) // fix segfault #892 (not 'channelIds[0] != 0')
    {
      fptype* totAllNumerators = ghelAllNumerators;     // reuse "helicity #0" buffer to compute the total over all helicities
      fptype* totAllDenominators = ghelAllDenominators; // reuse "helicity #0" buffer to compute the total over all helicities
      for( int ighel = 1; ighel < dcNGoodHel; ighel++ ) // NB: the loop starts at ighel=1
      {
        fptype* hAllNumerators = ghelAllNumerators + ighel * nevt;
        fptype* hAllDenominators = ghelAllDenominators + ighel * nevt;
        totAllNumerators[ievt] += hAllNumerators[ievt];
        totAllDenominators[ievt] += hAllDenominators[ievt];
      }
      allMEs[ievt] *= totAllNumerators[ievt] / totAllDenominators[ievt];
    }
#endif
    return;
  }
#endif

  //--------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL
  __global__ void
  add_and_select_hel( int* allselhel,          // output: helicity selection[nevt]
                      const fptype* allrndhel, // input: random numbers[nevt] for helicity selection
                      fptype* ghelAllMEs,      // input/tmp: allMEs for nGoodHel <= ncomb individual/runningsum helicities (index is ighel)
                      fptype* allMEs,          // output: allMEs[nevt], final sum over helicities
                      const int nevt )         // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
  {
    const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread)
    // Compute the sum of MEs over all good helicities (defer this after the helicity loop to avoid breaking streams parall>
    for( int ighel = 0; ighel < dcNGoodHel; ighel++ )
    {
      allMEs[ievt] += ghelAllMEs[ighel * nevt + ievt];
      ghelAllMEs[ighel * nevt + ievt] = allMEs[ievt]; // reuse the buffer to store the running sum for helicity selection
    }
    // Event-by-event random choice of helicity #403
    //printf( "select_hel: ievt=%4d rndhel=%f\n", ievt, allrndhel[ievt] );
    for( int ighel = 0; ighel < dcNGoodHel; ighel++ )
    {
      if( allrndhel[ievt] < ( ghelAllMEs[ighel * nevt + ievt] / allMEs[ievt] ) )
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
  update_jamp2s( const fptype_sv* allJamps, // input: jamp[ncolor*2*nevt] for this helicity
                 fptype* colAllJamp2s )     // output: allJamp2s[ncolor][nevt] super-buffer, sum over col/hel (nullptr to disable)
  {
    using J_ACCESS = DeviceAccessJamp;
    using J2_ACCESS = DeviceAccessJamp2;
    for( int icol = 0; icol < ncolor; icol++ )
      // NB: atomicAdd is needed after moving to cuda streams with one helicity per stream!
      atomicAdd( &J2_ACCESS::kernelAccessIcol( colAllJamp2s, icol ), cxabs2( J_ACCESS::kernelAccessIcolConst( allJamps, icol ) ) );
  }
#endif
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
      using J2_ACCESS = DeviceAccessJamp2;
      for( int icolC = 0; icolC < ncolor; icolC++ )
        jamp2_sv[icolC] = J2_ACCESS::kernelAccessIcolConst( allJamp2s, icolC );
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
      //printf( "select_col: ievt=%4d rndcol=%f\n", ievt, allrndcol[ievt] );
      for( int icolC = 0; icolC < ncolor; icolC++ )
      {
        if( allrndcol[ievt] < ( targetamp[icolC] / targetamp[ncolor - 1] ) )
        {
          allselcol[ievt] = icolC + 1; // NB Fortran [1,ncolor], cudacpp [0,ncolor-1]
          //printf( "select_col: ievt=%d icol=%d\n", ievt, icolC+1 );
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

#ifdef MGONGPUCPP_GPUIMPL /* clang-format off */
  void
  sigmaKin( const fptype* allmomenta,           // input: momenta[nevt*npar*4]
            const fptype* allcouplings,         // input: couplings[nevt*ndcoup*2]
            const fptype* allrndhel,            // input: random numbers[nevt] for helicity selection
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
            const fptype* allrndcol,            // input: random numbers[nevt] for color selection
            const unsigned int* allChannelIds,  // input: channelIds[nevt] (1 to #diagrams); nullptr to disable single-diagram enhancement (fix #899/#911)
#endif
            fptype* allMEs,                     // output: allMEs[nevt], |M|^2 final_avg_over_helicities
            int* allselhel,                     // output: helicity selection[nevt]
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
            int* allselcol,                     // output: helicity selection[nevt]
            fptype* colAllJamp2s,               // tmp: allJamp2s super-buffer for ncolor individual colors, running sum over colors and helicities
            fptype* ghelAllNumerators,          // tmp: allNumerators super-buffer for nGoodHel <= ncomb individual helicities (index is ighel)
            fptype* ghelAllDenominators,        // tmp: allDenominators super-buffer for nGoodHel <= ncomb individual helicities (index is ighel)
#endif
            fptype* ghelAllMEs,                 // tmp: allMEs super-buffer for nGoodHel <= ncomb individual helicities (index is ighel)
            fptype* ghelAllJamps,               // tmp: allJamps super-buffer for nGoodHel <= ncomb individual helicities (index is ighel)
            fptype* ghelAllWfs,                 // tmp: allWfs super-buffer for nGoodHel <= ncomb individual helicities (index is ighel)
            gpuStream_t* ghelStreams,           // input: cuda streams (index is ighel: only the first nGoodHel <= ncomb are non-null)
            const int gpublocks,                // input: cuda gpublocks
            const int gputhreads )              // input: cuda gputhreads
#else
  void
  sigmaKin( const fptype* allmomenta,           // input: momenta[nevt*npar*4]
            const fptype* allcouplings,         // input: couplings[nevt*ndcoup*2]
            const fptype* allrndhel,            // input: random numbers[nevt] for helicity selection
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
            const fptype* allrndcol,            // input: random numbers[nevt] for color selection
            const unsigned int* allChannelIds,  // input: channelIds[nevt] (1 to #diagrams); nullptr to disable single-diagram enhancement (fix #899/#911)
#endif
            fptype* allMEs,                     // output: allMEs[nevt], |M|^2 final_avg_over_helicities
            int* allselhel,                     // output: helicity selection[nevt]
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
            int* allselcol,                     // output: helicity selection[nevt]
            fptype* allNumerators,              // tmp: multichannel numerators[nevt], running_sum_over_helicities
            fptype* allDenominators,            // tmp: multichannel denominators[nevt], running_sum_over_helicities
#endif
            const int nevt                      // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
            )
#endif /* clang-format on */
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
    constexpr int helcolDenominators[1] = { 1536 }; // assume nprocesses == 1 (#272 and #343)

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

    // === PART 0 - INITIALISATION (before calculate_jamps) ===
    // Reset the "matrix elements" - running sums of |M|^2 over helicities for the given event
#ifdef MGONGPUCPP_GPUIMPL
    const int nevt = gpublocks * gputhreads;
    gpuMemset( allMEs, 0, nevt * sizeof( fptype ) );
    gpuMemset( ghelAllJamps, 0, cNGoodHel * ncolor * mgOnGpu::nx2 * nevt * sizeof( fptype ) );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    gpuMemset( colAllJamp2s, 0, ncolor * nevt * sizeof( fptype ) );
    gpuMemset( ghelAllNumerators, 0, cNGoodHel * nevt * sizeof( fptype ) );
    gpuMemset( ghelAllDenominators, 0, cNGoodHel * nevt * sizeof( fptype ) );
#endif
    gpuMemset( ghelAllMEs, 0, cNGoodHel * nevt * sizeof( fptype ) );
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
    // Use CUDA/HIP streams to process different helicities in parallel (one good helicity per stream)
    // (1a) First, within each helicity stream, compute the QCD partial amplitudes jamp's for each helicity
    // In multichannel mode, also compute the running sums over helicities of numerators, denominators and squared jamp2s
    for( int ighel = 0; ighel < cNGoodHel; ighel++ )
    {
      const int ihel = cGoodHel[ighel];
      fptype* hAllJamps = ghelAllJamps + ighel * nevt * ncolor * mgOnGpu::nx2;
      fptype* hAllWfs = ghelAllWfs + ighel * nwf * nevt * nw6 * mgOnGpu::nx2;
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      fptype* hAllNumerators = ghelAllNumerators + ighel * nevt;
      fptype* hAllDenominators = ghelAllDenominators + ighel * nevt;
      calculate_jamps( ihel, allmomenta, allcouplings, hAllJamps, hAllWfs, allChannelIds, hAllNumerators, hAllDenominators, ghelStreams[ighel], gpublocks, gputhreads );
#else
      calculate_jamps( ihel, allmomenta, allcouplings, hAllJamps, hAllWfs, ghelStreams[ighel], gpublocks, gputhreads );
#endif
    }
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // (1b) Then, In multichannel mode, also compute the running sums over helicities of squared jamp2s within each helicity stream
    for( int ighel = 0; ighel < cNGoodHel; ighel++ )
    {
      fptype* hAllJamps = ghelAllJamps + ighel * nevt * ncolor * mgOnGpu::nx2;
      gpuLaunchKernelStream( update_jamp2s, gpublocks, gputhreads, ghelStreams[ighel], hAllJamps, colAllJamp2s );
    }
#endif
    // (2) Then, within each helicity stream, compute the ME for that helicity from the color sum of QCD partial amplitudes jamps
    for( int ighel = 0; ighel < cNGoodHel; ighel++ )
    {
      fptype* hAllMEs = ghelAllMEs + ighel * nevt;
      fptype* hAllJamps = ghelAllJamps + ighel * nevt * ncolor * mgOnGpu::nx2;
      gpuLaunchKernelStream( color_sum, gpublocks, gputhreads, ghelStreams[ighel], hAllMEs, hAllJamps );
    }
    checkGpu( gpuDeviceSynchronize() ); // do not start helicity/color selection until the loop over helicities has completed
    // (3) Wait for all helicity streams to complete, then finally compute the ME sum over all helicities and choose one helicity and one color
    // Event-by-event random choice of helicity #403 and ME sum over helicities (defer this after the helicity loop to avoid breaking streams parallelism)
    gpuLaunchKernel( add_and_select_hel, gpublocks, gputhreads, allselhel, allrndhel, ghelAllMEs, allMEs, gpublocks * gputhreads );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Event-by-event random choice of color #402
    gpuLaunchKernel( select_col, gpublocks, gputhreads, allselcol, allrndcol, allChannelIds, colAllJamp2s, gpublocks * gputhreads );
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
      // **NB! in "mixed" precision, using SIMD, calculate_jamps computes MEs for TWO neppV pages with a single channelId! #924
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
        // **NB! in "mixed" precision, using SIMD, calculate_jamps computes MEs for TWO neppV pages with a single channelId! #924
        for( int i = 0; i < neppV; ++i )
        {
          assert( channelId == channelIds2_v[i] ); // SANITY CHECKS #898 #924: all events in the 2nd SIMD vector have the same channelId as that of the 1st SIMD vector
        }
#endif
      }
#endif
      // Running sum of partial amplitudes squared for event by event color selection (#402)
      // (jamp2[nParity][ncolor][neppV] for the SIMD vector - or the two SIMD vectors - of events processed in calculate_jamps)
      fptype_sv jamp2_sv[nParity * ncolor] = {};
      fptype_sv MEs_ighel[ncomb] = {};  // sum of MEs for all good helicities up to ighel (for the first - and/or only - neppV page)
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
      fptype_sv MEs_ighel2[ncomb] = {}; // sum of MEs for all good helicities up to ighel (for the second neppV page)
#endif
      for( int ighel = 0; ighel < cNGoodHel; ighel++ )
      {
        const int ihel = cGoodHel[ighel];
        cxtype_sv jamp_sv[nParity * ncolor] = {}; // fixed nasty bug (omitting 'nParity' caused memory corruptions after calling calculate_jamps)
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
        // **NB! in "mixed" precision, using SIMD, calculate_jamps computes MEs for TWO neppV pages with a single channelId! #924
        calculate_jamps( ihel, allmomenta, allcouplings, jamp_sv, channelId, allNumerators, allDenominators, ievt00 );
#else
        calculate_jamps( ihel, allmomenta, allcouplings, jamp_sv, ievt00 );
#endif
        color_sum( allMEs, jamp_sv, ievt00 );
        MEs_ighel[ighel] = E_ACCESS::kernelAccess( E_ACCESS::ieventAccessRecord( allMEs, ievt00 ) );
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
        MEs_ighel2[ighel] = E_ACCESS::kernelAccess( E_ACCESS::ieventAccessRecord( allMEs, ievt00 + neppV ) );
#endif
        using J_ACCESS = HostAccessJamp;
        for( int iParity = 0; iParity < nParity; ++iParity )
          for( int icol = 0; icol < ncolor; icol++ )
            jamp2_sv[ncolor * iParity + icol] += cxabs2( J_ACCESS::kernelAccessIcol( &( jamp_sv[ncolor * iParity] ), icol ) ); // may underflow #831
      }
      // Event-by-event random choice of helicity #403
      for( int ieppV = 0; ieppV < neppV; ++ieppV )
      {
        const int ievt = ievt00 + ieppV;
        //printf( "sigmaKin: ievt=%4d rndhel=%f\n", ievt, allrndhel[ievt] );
        for( int ighel = 0; ighel < cNGoodHel; ighel++ )
        {
#if defined MGONGPU_CPPSIMD
          //printf( "sigmaKin: ievt=%4d ighel=%d MEs_ighel=%f\n", ievt, ighel, MEs_ighel[ighel][ieppV] );
          const bool okhel = allrndhel[ievt] < ( MEs_ighel[ighel][ieppV] / MEs_ighel[cNGoodHel - 1][ieppV] );
#else
          //printf( "sigmaKin: ievt=%4d ighel=%d MEs_ighel=%f\n", ievt, ighel, MEs_ighel[ighel] );
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
          //printf( "sigmaKin: ievt=%4d ighel=%d MEs_ighel=%f\n", ievt2, ighel, MEs_ighel2[ighel][ieppV] );
          if( allrndhel[ievt2] < ( MEs_ighel2[ighel][ieppV] / MEs_ighel2[cNGoodHel - 1][ieppV] ) )
          {
            const int ihelF = cGoodHel[ighel] + 1; // NB Fortran [1,ncomb], cudacpp [0,ncomb-1]
            allselhel[ievt2] = ihelF;
            //printf( "sigmaKin: ievt=%4d ihel=%4d\n", ievt2, ihelF );
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

    // PART 2 - FINALISATION (after calculate_jamps)
    // Get the final |M|^2 as an average over helicities/colors of the running sum of |M|^2 over helicities for the given event
    // [NB 'sum over final spins, average over initial spins', eg see
    // https://www.uzh.ch/cmsssl/physik/dam/jcr:2e24b7b1-f4d7-4160-817e-47b13dbf1d7c/Handout_4_2016-UZH.pdf]
#ifdef MGONGPUCPP_GPUIMPL
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    gpuLaunchKernel( normalise_output, gpublocks, gputhreads, allMEs, ghelAllNumerators, ghelAllDenominators, allChannelIds, helcolDenominators[0] );
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
