// Copyright (C) 2010 The MadGraph5_aMC@NLO development team and contributors.
// Created by: J. Alwall (Oct 2010) for the MG5aMC CPP backend.
//==========================================================================
// Copyright (C) 2020-2025 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Modified by: S. Roiser (Feb 2020) for the MG5aMC CUDACPP plugin.
// Further modified by: S. Hageboeck, O. Mattelaer, S. Roiser, J. Teig, A. Valassi, Z. Wettersten (2020-2025) for the MG5aMC CUDACPP plugin.
//==========================================================================
// This file has been automatically generated for CUDA/C++ standalone by
// MadGraph5_aMC@NLO v. 3.6.4, 2025-09-13
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "CPPProcess.h"

#include "mgOnGpuConfig.h"

#include "GpuRuntime.h"
#include "HelAmps_sm.h"
#include "MemoryAccessAmplitudes.h"
#include "MemoryAccessChannelIds.h"
#include "MemoryAccessCouplings.h"
#include "MemoryAccessCouplingsFixed.h"
#include "MemoryAccessGs.h"
#include "MemoryAccessMatrixElements.h"
#include "MemoryAccessMomenta.h"
#include "MemoryAccessWavefunctions.h"
#include "color_sum.h"
#include "diagrams.h"

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
// Process: g g > t t~ g g g g WEIGHTED<=6 @1

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

#ifndef MGONGPU_RDC_DIAGRAMS
  constexpr int ndiagramgroups = CPPProcess::ndiagramgroups; // the number of Feynman diagram groups
#endif

  using Parameters_sm_dependentCouplings::ndcoup;   // #couplings that vary event by event (depend on running alphas QCD)
  using Parameters_sm_independentCouplings::nicoup; // #couplings that are fixed for all events (do not depend on running alphas QCD)

  constexpr int nIPD = CPPProcess::nIPD; // SM independent parameters
  constexpr int nIPC = CPPProcess::nIPC; // SM independent couplings

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
  static_assert( nIPC <= nicoup );
  static_assert( nIPD >= 0 ); // Hack to avoid build warnings when nIPD==0 is unused
  static_assert( nIPC >= 0 ); // Hack to avoid build warnings when nIPC==0 is unused
  // Hardcoded parameters (HRDCOD=1)
#ifdef MGONGPU_HARDCODE_PARAM
  __device__ const fptype dcIPD[nIPD] = { (fptype)Parameters_sm::mdl_MT, (fptype)Parameters_sm::mdl_WT };
  __device__ const fptype* dcIPC = nullptr; // unused as nIPC=0
#ifdef MGONGPUCPP_GPUIMPL
  static fptype* cIPD = nullptr; // symbol address
  static fptype* cIPC = nullptr; // symbol address
#else
  static const fptype* cIPD = dcIPD;
  static const fptype* cIPC = dcIPC;
#endif
  // Non-hardcoded parameters (HRDCOD=0)
#else
#ifdef MGONGPUCPP_GPUIMPL /* clang-format off */
  __device__ __constant__ fptype dcIPD[nIPD];
  __device__ __constant__ fptype* dcIPC = nullptr; // unused as nIPC=0
  static fptype* cIPD = nullptr; // symbol address
  static fptype* cIPC = nullptr; // symbol address
#else /* clang-format on */
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
  __device__ __constant__ short dcHel[ncomb][npar];
  __device__ __constant__ int dcNGoodHel;
  __device__ __constant__ int dcGoodHel[ncomb];
#ifndef MGONGPU_RDC_DIAGRAMS
  static short* cHelFlat = nullptr; // symbol address
#endif
#else
  static short cHel[ncomb][npar];
  static short* cHelFlat = (short*)cHel;
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
  // Launch a group of Feynman diagrams as a standalone kernel (sigmaKin_getGoodHel) or within a CUDA/HIP graph (sigmaKin)
  template<typename Func, typename... Args>
  void
  gpuDiagrams( bool useGraphs,
               gpuGraph_t* pGraph,
               gpuGraphExec_t* pGraphExec,
               gpuGraphNode_t* pNode,
               gpuGraphNode_t* pNodeDep,
               Func diagrams,
               int gpublocks,
               int gputhreads,
               gpuStream_t gpustream,
               Args... args )
  {
    // CASE 0: WITHOUT GRAPHS (graphs disabled)
    if( !useGraphs )
    {
      gpuLaunchKernelStream( diagrams, gpublocks, gputhreads, gpustream, args... );
    }
    // CASE 0: WITHOUT GRAPHS (graphs enabled - sigmaKin_getGoodHel)
    else if( gpustream == 0 )
    {
      gpuLaunchKernelStream( diagrams, gpublocks, gputhreads, gpustream, args... );
    }
    // CASE 1: WITH GRAPHS (graphs enabled - sigmaKin)
    else
    {
      // Define the parameters for the graph node for this Feynman diagram
      gpuKernelNodeParams params = {};
      void* kParams[] = { static_cast<void*>( &args )... };
      params.func = (void*)diagrams;
      params.gridDim = dim3( gpublocks );
      params.blockDim = dim3( gputhreads );
      params.kernelParams = kParams;
      // Create the graph node for this Feynman diagram if not yet done
      if( !( *pNode ) )
      {
        if( pNodeDep == nullptr )
        {
          checkGpu( gpuGraphAddKernelNode( pNode, *pGraph, nullptr, 0, &params ) );
          //std::cout << "Added graph node " << pNode << " with no dependencies" << std::endl;
        }
        else
        {
          checkGpu( gpuGraphAddKernelNode( pNode, *pGraph, pNodeDep, 1, &params ) );
          //std::cout << "Added graph node " << pNode << " with one dependency on " << pNodeDep << std::endl;
        }
      }
      // Update parameters if the graph node for this Feynman diagram already exists
      else
      {
        checkGpu( gpuGraphExecKernelNodeSetParams( *pGraphExec, *pNode, &params ) );
        //std::cout << "Updated parameters for graph node " << pNode << std::endl;
      }
    }
  }
#endif

  //--------------------------------------------------------------------------

  // Evaluate QCD partial amplitudes jamps for this given helicity from Feynman diagrams
  // Also compute running sums over helicities adding jamp2, numerator, denominator
  // (NB: this function no longer handles matrix elements as the color sum has now been moved to a separate function/kernel)
#ifdef MGONGPUCPP_GPUIMPL /* clang-format off */
  // In CUDA, this function processes a single event
  // ** NB1: NEW Nov2024! In CUDA this is now a kernel function (it used to be a device function)
  // ** NB2: NEW Nov2024! in CUDA this now takes a channelId array as input (it used to take a scalar channelId as input)
#ifndef MGONGPU_RDC_DIAGRAMS
  INLINE void
#else
  __global__ void
#endif
  calculate_jamps( const fptype* allmomenta,          // input: momenta[nevt*npar*4]
                   const fptype* allcouplings,        // input: couplings[nevt*ndcoup*2]
                   fptype* allJamps,                  // output: jamp[ncolor*2*nevt] for this helicity
                   fptype* allWfs,                    // output: wf[nwf*nw6*2*nevt]
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
                   const unsigned int* allChannelIds, // input: multichannel channelIds[nevt] (1 to #diagrams); nullptr to disable SDE (#899/#911)
                   fptype* allNumerators,             // input/output: multichannel numerators[nevt], add helicity ihel
                   fptype* allDenominators,           // input/output: multichannel denominators[nevt], add helicity ihel
#endif
#ifndef MGONGPU_RDC_DIAGRAMS
                   gpuStream_t gpustream,             // input: cuda stream for this helicity
                   const int gpublocks,               // input: cuda gpublocks
                   const int gputhreads,              // input: cuda gputhreads
#endif
                   int ihel )
#else
  // In C++, this function processes a single event "page" or SIMD vector (or for two in "mixed" precision mode, nParity=2)
  // *** NB: in C++, calculate_jamps accepts a SCALAR channelId because it is GUARANTEED that all events in a SIMD vector have the same channelId #898
  INLINE void
  calculate_jamps( const fptype* allmomenta,          // input: momenta[nevt*npar*4]
                   const fptype* allcouplings,        // input: couplings[nevt*ndcoup*2]
                   cxtype_sv* jamp_sv_1or2,           // output: jamp_sv[ncolor] (f/d) or [2*ncolor] (m) for SIMD event page(s) ievt00 and helicity ihel
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
                   const unsigned int channelId,      // input: SCALAR channelId (1 to #diagrams, 0 to disable SDE) for SIMD event page(s) ievt00
                   fptype* allNumerators,             // input/output: multichannel numerators[nevt], add helicity ihel
                   fptype* allDenominators,           // input/output: multichannel denominators[nevt], add helicity ihel
#endif
                   const int ievt00,                  // input: first event number in current C++ event page (for CUDA, ievt depends on threadid)
                   int ihel )
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
    // Local TEMPORARY variables for a subset of Feynman diagrams in the given C++ event page (ipagV)
    // [NB these variables are reused several times (and re-initialised each time) within the same event or event page]
    // ** NB: wavefunctions only need TRIVIAL ACCESS in C++ code
    cxtype_sv w_sv[nwf][nw6]; // particle wavefunctions within Feynman diagrams (nw6 is often 6, the dimension of spin 1/2 or spin 1 particles)
    fptype* wfs = reinterpret_cast<fptype*>( w_sv );
#else
#ifndef MGONGPU_RDC_DIAGRAMS
    // Global-memory variables for a subset of Feynman diagrams in the given CUDA event (ievt)
    // ** NB: wavefunctions need non-trivial access in CUDA code because of kernel splitting
    fptype* wfs = allWfs;
#else
    // Local TEMPORARY variables for a subset of Feynman diagrams in the given CUDA event (ievt)
    // [NB these variables are reused several times (and re-initialised each time) within the same event or event page]
    // ** NB: wavefunctions only need TRIVIAL ACCESS in C++ code
    assert( allWfs == nullptr ); // sanity check
    cxtype_sv w_sv[nwf][nw6];    // particle wavefunctions within Feynman diagrams (nw6 is often 6, the dimension of spin 1/2 or spin 1 particles)
    fptype* wfs = reinterpret_cast<fptype*>( w_sv );
#endif
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
      fptype* jamps = allJamps;
#else
      // In C++, write jamps to the output array [for one specific event or SIMD vector] passed as argument
      cxtype_sv* jamp_sv = ( iParity == 0 ? jamp_sv_1or2 : &( jamp_sv_1or2[ncolor] ) );
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

      // *** DIAGRAMS 1 TO 15495 ***
#ifdef MGONGPUCPP_GPUIMPL
#ifndef MGONGPU_RDC_DIAGRAMS
      // === GPU IMPLEMENTATION (DCDIAG=0): each diagram group is an individual kernel ===
      static bool useGraphs = false;
      static bool first = true;
      if( first )
      {
        first = false;
        // Analyse environment variable CUDACPP_RUNTIME_GPUGRAPHS
        const char* graphsEnv = getenv( "CUDACPP_RUNTIME_GPUGRAPHS" );
        if( graphsEnv && std::string( graphsEnv ) != "" )
        {
          useGraphs = true;
          std::cout << "INFO: Env variable CUDACPP_RUNTIME_GPUGRAPHS is set and non-empty: use GPU Graphs" << std::endl;
        }
        else
        {
          std::cout << "INFO: Env variable CUDACPP_RUNTIME_GPUGRAPHS is empty or not set: do not use GPU Graphs" << std::endl;
        }
      }
      static gpuGraph_t graphs[ncomb] = {};
      static gpuGraphExec_t graphExecs[ncomb] = {};
      static gpuGraphNode_t graphNodes[ncomb * ndiagramgroups] = {};
      gpuGraph_t& graph = graphs[ihel];
      gpuGraphExec_t& graphExec = graphExecs[ihel];
      // Case 1 with graphs (gpustream!=0, sigmaKin): create the graph if not yet done
      if( useGraphs && gpustream != 0 )
      {
        if( !graph )
        {
          checkGpu( gpuGraphCreate( &graph, 0 ) );
          //std::cout << "(ihel=" << ihel << ") Created graph " << graph << std::endl;
        }
      }
      // Case 0 without graphs (gpustream==0, sigmaKin_getGoodHel): launch all diagram kernels
      // Case 1 with graphs (gpustream!=0, sigmaKin): create graph nodes if not yet done, else update them with new parameters
      gpuGraphNode_t& node1 = graphNodes[ihel * ndiagramgroups + 0];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1, nullptr, diagramgroup1, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD, cHelFlat, momenta, ihel );
      gpuGraphNode_t& node2 = graphNodes[ihel * ndiagramgroups + 1];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node2, &node1, diagramgroup2, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node3 = graphNodes[ihel * ndiagramgroups + 2];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node3, &node2, diagramgroup3, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node4 = graphNodes[ihel * ndiagramgroups + 3];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node4, &node3, diagramgroup4, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node5 = graphNodes[ihel * ndiagramgroups + 4];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node5, &node4, diagramgroup5, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node6 = graphNodes[ihel * ndiagramgroups + 5];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node6, &node5, diagramgroup6, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node7 = graphNodes[ihel * ndiagramgroups + 6];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node7, &node6, diagramgroup7, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node8 = graphNodes[ihel * ndiagramgroups + 7];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node8, &node7, diagramgroup8, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node9 = graphNodes[ihel * ndiagramgroups + 8];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node9, &node8, diagramgroup9, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node10 = graphNodes[ihel * ndiagramgroups + 9];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node10, &node9, diagramgroup10, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node11 = graphNodes[ihel * ndiagramgroups + 10];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node11, &node10, diagramgroup11, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node12 = graphNodes[ihel * ndiagramgroups + 11];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node12, &node11, diagramgroup12, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node13 = graphNodes[ihel * ndiagramgroups + 12];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node13, &node12, diagramgroup13, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node14 = graphNodes[ihel * ndiagramgroups + 13];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node14, &node13, diagramgroup14, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node15 = graphNodes[ihel * ndiagramgroups + 14];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node15, &node14, diagramgroup15, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node16 = graphNodes[ihel * ndiagramgroups + 15];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node16, &node15, diagramgroup16, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node17 = graphNodes[ihel * ndiagramgroups + 16];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node17, &node16, diagramgroup17, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node18 = graphNodes[ihel * ndiagramgroups + 17];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node18, &node17, diagramgroup18, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node19 = graphNodes[ihel * ndiagramgroups + 18];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node19, &node18, diagramgroup19, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node20 = graphNodes[ihel * ndiagramgroups + 19];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node20, &node19, diagramgroup20, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node21 = graphNodes[ihel * ndiagramgroups + 20];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node21, &node20, diagramgroup21, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node22 = graphNodes[ihel * ndiagramgroups + 21];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node22, &node21, diagramgroup22, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node23 = graphNodes[ihel * ndiagramgroups + 22];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node23, &node22, diagramgroup23, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node24 = graphNodes[ihel * ndiagramgroups + 23];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node24, &node23, diagramgroup24, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node25 = graphNodes[ihel * ndiagramgroups + 24];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node25, &node24, diagramgroup25, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node26 = graphNodes[ihel * ndiagramgroups + 25];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node26, &node25, diagramgroup26, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node27 = graphNodes[ihel * ndiagramgroups + 26];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node27, &node26, diagramgroup27, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node28 = graphNodes[ihel * ndiagramgroups + 27];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node28, &node27, diagramgroup28, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node29 = graphNodes[ihel * ndiagramgroups + 28];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node29, &node28, diagramgroup29, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node30 = graphNodes[ihel * ndiagramgroups + 29];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node30, &node29, diagramgroup30, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node31 = graphNodes[ihel * ndiagramgroups + 30];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node31, &node30, diagramgroup31, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node32 = graphNodes[ihel * ndiagramgroups + 31];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node32, &node31, diagramgroup32, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node33 = graphNodes[ihel * ndiagramgroups + 32];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node33, &node32, diagramgroup33, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node34 = graphNodes[ihel * ndiagramgroups + 33];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node34, &node33, diagramgroup34, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node35 = graphNodes[ihel * ndiagramgroups + 34];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node35, &node34, diagramgroup35, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node36 = graphNodes[ihel * ndiagramgroups + 35];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node36, &node35, diagramgroup36, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node37 = graphNodes[ihel * ndiagramgroups + 36];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node37, &node36, diagramgroup37, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node38 = graphNodes[ihel * ndiagramgroups + 37];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node38, &node37, diagramgroup38, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node39 = graphNodes[ihel * ndiagramgroups + 38];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node39, &node38, diagramgroup39, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node40 = graphNodes[ihel * ndiagramgroups + 39];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node40, &node39, diagramgroup40, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node41 = graphNodes[ihel * ndiagramgroups + 40];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node41, &node40, diagramgroup41, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node42 = graphNodes[ihel * ndiagramgroups + 41];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node42, &node41, diagramgroup42, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node43 = graphNodes[ihel * ndiagramgroups + 42];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node43, &node42, diagramgroup43, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node44 = graphNodes[ihel * ndiagramgroups + 43];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node44, &node43, diagramgroup44, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node45 = graphNodes[ihel * ndiagramgroups + 44];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node45, &node44, diagramgroup45, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node46 = graphNodes[ihel * ndiagramgroups + 45];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node46, &node45, diagramgroup46, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node47 = graphNodes[ihel * ndiagramgroups + 46];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node47, &node46, diagramgroup47, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node48 = graphNodes[ihel * ndiagramgroups + 47];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node48, &node47, diagramgroup48, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node49 = graphNodes[ihel * ndiagramgroups + 48];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node49, &node48, diagramgroup49, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node50 = graphNodes[ihel * ndiagramgroups + 49];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node50, &node49, diagramgroup50, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node51 = graphNodes[ihel * ndiagramgroups + 50];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node51, &node50, diagramgroup51, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node52 = graphNodes[ihel * ndiagramgroups + 51];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node52, &node51, diagramgroup52, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node53 = graphNodes[ihel * ndiagramgroups + 52];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node53, &node52, diagramgroup53, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node54 = graphNodes[ihel * ndiagramgroups + 53];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node54, &node53, diagramgroup54, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node55 = graphNodes[ihel * ndiagramgroups + 54];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node55, &node54, diagramgroup55, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node56 = graphNodes[ihel * ndiagramgroups + 55];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node56, &node55, diagramgroup56, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node57 = graphNodes[ihel * ndiagramgroups + 56];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node57, &node56, diagramgroup57, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node58 = graphNodes[ihel * ndiagramgroups + 57];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node58, &node57, diagramgroup58, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node59 = graphNodes[ihel * ndiagramgroups + 58];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node59, &node58, diagramgroup59, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node60 = graphNodes[ihel * ndiagramgroups + 59];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node60, &node59, diagramgroup60, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node61 = graphNodes[ihel * ndiagramgroups + 60];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node61, &node60, diagramgroup61, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node62 = graphNodes[ihel * ndiagramgroups + 61];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node62, &node61, diagramgroup62, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node63 = graphNodes[ihel * ndiagramgroups + 62];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node63, &node62, diagramgroup63, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node64 = graphNodes[ihel * ndiagramgroups + 63];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node64, &node63, diagramgroup64, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node65 = graphNodes[ihel * ndiagramgroups + 64];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node65, &node64, diagramgroup65, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node66 = graphNodes[ihel * ndiagramgroups + 65];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node66, &node65, diagramgroup66, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node67 = graphNodes[ihel * ndiagramgroups + 66];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node67, &node66, diagramgroup67, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node68 = graphNodes[ihel * ndiagramgroups + 67];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node68, &node67, diagramgroup68, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node69 = graphNodes[ihel * ndiagramgroups + 68];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node69, &node68, diagramgroup69, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node70 = graphNodes[ihel * ndiagramgroups + 69];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node70, &node69, diagramgroup70, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node71 = graphNodes[ihel * ndiagramgroups + 70];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node71, &node70, diagramgroup71, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node72 = graphNodes[ihel * ndiagramgroups + 71];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node72, &node71, diagramgroup72, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node73 = graphNodes[ihel * ndiagramgroups + 72];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node73, &node72, diagramgroup73, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node74 = graphNodes[ihel * ndiagramgroups + 73];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node74, &node73, diagramgroup74, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node75 = graphNodes[ihel * ndiagramgroups + 74];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node75, &node74, diagramgroup75, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node76 = graphNodes[ihel * ndiagramgroups + 75];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node76, &node75, diagramgroup76, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node77 = graphNodes[ihel * ndiagramgroups + 76];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node77, &node76, diagramgroup77, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node78 = graphNodes[ihel * ndiagramgroups + 77];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node78, &node77, diagramgroup78, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node79 = graphNodes[ihel * ndiagramgroups + 78];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node79, &node78, diagramgroup79, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node80 = graphNodes[ihel * ndiagramgroups + 79];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node80, &node79, diagramgroup80, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node81 = graphNodes[ihel * ndiagramgroups + 80];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node81, &node80, diagramgroup81, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node82 = graphNodes[ihel * ndiagramgroups + 81];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node82, &node81, diagramgroup82, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node83 = graphNodes[ihel * ndiagramgroups + 82];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node83, &node82, diagramgroup83, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node84 = graphNodes[ihel * ndiagramgroups + 83];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node84, &node83, diagramgroup84, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node85 = graphNodes[ihel * ndiagramgroups + 84];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node85, &node84, diagramgroup85, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node86 = graphNodes[ihel * ndiagramgroups + 85];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node86, &node85, diagramgroup86, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node87 = graphNodes[ihel * ndiagramgroups + 86];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node87, &node86, diagramgroup87, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node88 = graphNodes[ihel * ndiagramgroups + 87];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node88, &node87, diagramgroup88, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node89 = graphNodes[ihel * ndiagramgroups + 88];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node89, &node88, diagramgroup89, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node90 = graphNodes[ihel * ndiagramgroups + 89];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node90, &node89, diagramgroup90, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node91 = graphNodes[ihel * ndiagramgroups + 90];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node91, &node90, diagramgroup91, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node92 = graphNodes[ihel * ndiagramgroups + 91];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node92, &node91, diagramgroup92, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node93 = graphNodes[ihel * ndiagramgroups + 92];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node93, &node92, diagramgroup93, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node94 = graphNodes[ihel * ndiagramgroups + 93];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node94, &node93, diagramgroup94, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node95 = graphNodes[ihel * ndiagramgroups + 94];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node95, &node94, diagramgroup95, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node96 = graphNodes[ihel * ndiagramgroups + 95];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node96, &node95, diagramgroup96, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node97 = graphNodes[ihel * ndiagramgroups + 96];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node97, &node96, diagramgroup97, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node98 = graphNodes[ihel * ndiagramgroups + 97];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node98, &node97, diagramgroup98, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node99 = graphNodes[ihel * ndiagramgroups + 98];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node99, &node98, diagramgroup99, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node100 = graphNodes[ihel * ndiagramgroups + 99];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node100, &node99, diagramgroup100, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node101 = graphNodes[ihel * ndiagramgroups + 100];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node101, &node100, diagramgroup101, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node102 = graphNodes[ihel * ndiagramgroups + 101];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node102, &node101, diagramgroup102, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node103 = graphNodes[ihel * ndiagramgroups + 102];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node103, &node102, diagramgroup103, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node104 = graphNodes[ihel * ndiagramgroups + 103];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node104, &node103, diagramgroup104, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node105 = graphNodes[ihel * ndiagramgroups + 104];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node105, &node104, diagramgroup105, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node106 = graphNodes[ihel * ndiagramgroups + 105];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node106, &node105, diagramgroup106, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node107 = graphNodes[ihel * ndiagramgroups + 106];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node107, &node106, diagramgroup107, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node108 = graphNodes[ihel * ndiagramgroups + 107];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node108, &node107, diagramgroup108, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node109 = graphNodes[ihel * ndiagramgroups + 108];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node109, &node108, diagramgroup109, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node110 = graphNodes[ihel * ndiagramgroups + 109];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node110, &node109, diagramgroup110, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node111 = graphNodes[ihel * ndiagramgroups + 110];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node111, &node110, diagramgroup111, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node112 = graphNodes[ihel * ndiagramgroups + 111];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node112, &node111, diagramgroup112, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node113 = graphNodes[ihel * ndiagramgroups + 112];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node113, &node112, diagramgroup113, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node114 = graphNodes[ihel * ndiagramgroups + 113];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node114, &node113, diagramgroup114, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node115 = graphNodes[ihel * ndiagramgroups + 114];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node115, &node114, diagramgroup115, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node116 = graphNodes[ihel * ndiagramgroups + 115];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node116, &node115, diagramgroup116, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node117 = graphNodes[ihel * ndiagramgroups + 116];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node117, &node116, diagramgroup117, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node118 = graphNodes[ihel * ndiagramgroups + 117];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node118, &node117, diagramgroup118, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node119 = graphNodes[ihel * ndiagramgroups + 118];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node119, &node118, diagramgroup119, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node120 = graphNodes[ihel * ndiagramgroups + 119];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node120, &node119, diagramgroup120, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node121 = graphNodes[ihel * ndiagramgroups + 120];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node121, &node120, diagramgroup121, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node122 = graphNodes[ihel * ndiagramgroups + 121];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node122, &node121, diagramgroup122, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node123 = graphNodes[ihel * ndiagramgroups + 122];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node123, &node122, diagramgroup123, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node124 = graphNodes[ihel * ndiagramgroups + 123];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node124, &node123, diagramgroup124, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node125 = graphNodes[ihel * ndiagramgroups + 124];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node125, &node124, diagramgroup125, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node126 = graphNodes[ihel * ndiagramgroups + 125];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node126, &node125, diagramgroup126, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node127 = graphNodes[ihel * ndiagramgroups + 126];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node127, &node126, diagramgroup127, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node128 = graphNodes[ihel * ndiagramgroups + 127];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node128, &node127, diagramgroup128, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node129 = graphNodes[ihel * ndiagramgroups + 128];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node129, &node128, diagramgroup129, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node130 = graphNodes[ihel * ndiagramgroups + 129];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node130, &node129, diagramgroup130, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node131 = graphNodes[ihel * ndiagramgroups + 130];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node131, &node130, diagramgroup131, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node132 = graphNodes[ihel * ndiagramgroups + 131];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node132, &node131, diagramgroup132, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node133 = graphNodes[ihel * ndiagramgroups + 132];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node133, &node132, diagramgroup133, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node134 = graphNodes[ihel * ndiagramgroups + 133];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node134, &node133, diagramgroup134, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node135 = graphNodes[ihel * ndiagramgroups + 134];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node135, &node134, diagramgroup135, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node136 = graphNodes[ihel * ndiagramgroups + 135];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node136, &node135, diagramgroup136, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node137 = graphNodes[ihel * ndiagramgroups + 136];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node137, &node136, diagramgroup137, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node138 = graphNodes[ihel * ndiagramgroups + 137];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node138, &node137, diagramgroup138, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node139 = graphNodes[ihel * ndiagramgroups + 138];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node139, &node138, diagramgroup139, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node140 = graphNodes[ihel * ndiagramgroups + 139];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node140, &node139, diagramgroup140, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node141 = graphNodes[ihel * ndiagramgroups + 140];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node141, &node140, diagramgroup141, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node142 = graphNodes[ihel * ndiagramgroups + 141];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node142, &node141, diagramgroup142, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node143 = graphNodes[ihel * ndiagramgroups + 142];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node143, &node142, diagramgroup143, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node144 = graphNodes[ihel * ndiagramgroups + 143];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node144, &node143, diagramgroup144, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node145 = graphNodes[ihel * ndiagramgroups + 144];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node145, &node144, diagramgroup145, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node146 = graphNodes[ihel * ndiagramgroups + 145];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node146, &node145, diagramgroup146, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node147 = graphNodes[ihel * ndiagramgroups + 146];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node147, &node146, diagramgroup147, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node148 = graphNodes[ihel * ndiagramgroups + 147];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node148, &node147, diagramgroup148, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node149 = graphNodes[ihel * ndiagramgroups + 148];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node149, &node148, diagramgroup149, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node150 = graphNodes[ihel * ndiagramgroups + 149];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node150, &node149, diagramgroup150, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node151 = graphNodes[ihel * ndiagramgroups + 150];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node151, &node150, diagramgroup151, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node152 = graphNodes[ihel * ndiagramgroups + 151];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node152, &node151, diagramgroup152, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node153 = graphNodes[ihel * ndiagramgroups + 152];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node153, &node152, diagramgroup153, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node154 = graphNodes[ihel * ndiagramgroups + 153];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node154, &node153, diagramgroup154, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node155 = graphNodes[ihel * ndiagramgroups + 154];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node155, &node154, diagramgroup155, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node156 = graphNodes[ihel * ndiagramgroups + 155];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node156, &node155, diagramgroup156, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node157 = graphNodes[ihel * ndiagramgroups + 156];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node157, &node156, diagramgroup157, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node158 = graphNodes[ihel * ndiagramgroups + 157];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node158, &node157, diagramgroup158, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node159 = graphNodes[ihel * ndiagramgroups + 158];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node159, &node158, diagramgroup159, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node160 = graphNodes[ihel * ndiagramgroups + 159];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node160, &node159, diagramgroup160, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node161 = graphNodes[ihel * ndiagramgroups + 160];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node161, &node160, diagramgroup161, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node162 = graphNodes[ihel * ndiagramgroups + 161];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node162, &node161, diagramgroup162, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node163 = graphNodes[ihel * ndiagramgroups + 162];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node163, &node162, diagramgroup163, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node164 = graphNodes[ihel * ndiagramgroups + 163];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node164, &node163, diagramgroup164, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node165 = graphNodes[ihel * ndiagramgroups + 164];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node165, &node164, diagramgroup165, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node166 = graphNodes[ihel * ndiagramgroups + 165];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node166, &node165, diagramgroup166, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node167 = graphNodes[ihel * ndiagramgroups + 166];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node167, &node166, diagramgroup167, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node168 = graphNodes[ihel * ndiagramgroups + 167];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node168, &node167, diagramgroup168, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node169 = graphNodes[ihel * ndiagramgroups + 168];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node169, &node168, diagramgroup169, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node170 = graphNodes[ihel * ndiagramgroups + 169];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node170, &node169, diagramgroup170, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node171 = graphNodes[ihel * ndiagramgroups + 170];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node171, &node170, diagramgroup171, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node172 = graphNodes[ihel * ndiagramgroups + 171];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node172, &node171, diagramgroup172, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node173 = graphNodes[ihel * ndiagramgroups + 172];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node173, &node172, diagramgroup173, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node174 = graphNodes[ihel * ndiagramgroups + 173];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node174, &node173, diagramgroup174, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node175 = graphNodes[ihel * ndiagramgroups + 174];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node175, &node174, diagramgroup175, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node176 = graphNodes[ihel * ndiagramgroups + 175];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node176, &node175, diagramgroup176, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node177 = graphNodes[ihel * ndiagramgroups + 176];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node177, &node176, diagramgroup177, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node178 = graphNodes[ihel * ndiagramgroups + 177];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node178, &node177, diagramgroup178, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node179 = graphNodes[ihel * ndiagramgroups + 178];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node179, &node178, diagramgroup179, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node180 = graphNodes[ihel * ndiagramgroups + 179];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node180, &node179, diagramgroup180, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node181 = graphNodes[ihel * ndiagramgroups + 180];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node181, &node180, diagramgroup181, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node182 = graphNodes[ihel * ndiagramgroups + 181];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node182, &node181, diagramgroup182, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node183 = graphNodes[ihel * ndiagramgroups + 182];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node183, &node182, diagramgroup183, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node184 = graphNodes[ihel * ndiagramgroups + 183];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node184, &node183, diagramgroup184, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node185 = graphNodes[ihel * ndiagramgroups + 184];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node185, &node184, diagramgroup185, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node186 = graphNodes[ihel * ndiagramgroups + 185];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node186, &node185, diagramgroup186, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node187 = graphNodes[ihel * ndiagramgroups + 186];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node187, &node186, diagramgroup187, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node188 = graphNodes[ihel * ndiagramgroups + 187];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node188, &node187, diagramgroup188, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node189 = graphNodes[ihel * ndiagramgroups + 188];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node189, &node188, diagramgroup189, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node190 = graphNodes[ihel * ndiagramgroups + 189];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node190, &node189, diagramgroup190, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node191 = graphNodes[ihel * ndiagramgroups + 190];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node191, &node190, diagramgroup191, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node192 = graphNodes[ihel * ndiagramgroups + 191];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node192, &node191, diagramgroup192, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node193 = graphNodes[ihel * ndiagramgroups + 192];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node193, &node192, diagramgroup193, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node194 = graphNodes[ihel * ndiagramgroups + 193];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node194, &node193, diagramgroup194, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node195 = graphNodes[ihel * ndiagramgroups + 194];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node195, &node194, diagramgroup195, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node196 = graphNodes[ihel * ndiagramgroups + 195];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node196, &node195, diagramgroup196, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node197 = graphNodes[ihel * ndiagramgroups + 196];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node197, &node196, diagramgroup197, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node198 = graphNodes[ihel * ndiagramgroups + 197];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node198, &node197, diagramgroup198, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node199 = graphNodes[ihel * ndiagramgroups + 198];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node199, &node198, diagramgroup199, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node200 = graphNodes[ihel * ndiagramgroups + 199];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node200, &node199, diagramgroup200, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node201 = graphNodes[ihel * ndiagramgroups + 200];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node201, &node200, diagramgroup201, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node202 = graphNodes[ihel * ndiagramgroups + 201];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node202, &node201, diagramgroup202, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node203 = graphNodes[ihel * ndiagramgroups + 202];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node203, &node202, diagramgroup203, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node204 = graphNodes[ihel * ndiagramgroups + 203];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node204, &node203, diagramgroup204, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node205 = graphNodes[ihel * ndiagramgroups + 204];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node205, &node204, diagramgroup205, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node206 = graphNodes[ihel * ndiagramgroups + 205];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node206, &node205, diagramgroup206, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node207 = graphNodes[ihel * ndiagramgroups + 206];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node207, &node206, diagramgroup207, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node208 = graphNodes[ihel * ndiagramgroups + 207];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node208, &node207, diagramgroup208, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node209 = graphNodes[ihel * ndiagramgroups + 208];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node209, &node208, diagramgroup209, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node210 = graphNodes[ihel * ndiagramgroups + 209];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node210, &node209, diagramgroup210, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node211 = graphNodes[ihel * ndiagramgroups + 210];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node211, &node210, diagramgroup211, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node212 = graphNodes[ihel * ndiagramgroups + 211];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node212, &node211, diagramgroup212, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node213 = graphNodes[ihel * ndiagramgroups + 212];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node213, &node212, diagramgroup213, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node214 = graphNodes[ihel * ndiagramgroups + 213];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node214, &node213, diagramgroup214, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node215 = graphNodes[ihel * ndiagramgroups + 214];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node215, &node214, diagramgroup215, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node216 = graphNodes[ihel * ndiagramgroups + 215];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node216, &node215, diagramgroup216, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node217 = graphNodes[ihel * ndiagramgroups + 216];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node217, &node216, diagramgroup217, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node218 = graphNodes[ihel * ndiagramgroups + 217];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node218, &node217, diagramgroup218, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node219 = graphNodes[ihel * ndiagramgroups + 218];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node219, &node218, diagramgroup219, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node220 = graphNodes[ihel * ndiagramgroups + 219];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node220, &node219, diagramgroup220, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node221 = graphNodes[ihel * ndiagramgroups + 220];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node221, &node220, diagramgroup221, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node222 = graphNodes[ihel * ndiagramgroups + 221];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node222, &node221, diagramgroup222, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node223 = graphNodes[ihel * ndiagramgroups + 222];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node223, &node222, diagramgroup223, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node224 = graphNodes[ihel * ndiagramgroups + 223];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node224, &node223, diagramgroup224, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node225 = graphNodes[ihel * ndiagramgroups + 224];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node225, &node224, diagramgroup225, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node226 = graphNodes[ihel * ndiagramgroups + 225];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node226, &node225, diagramgroup226, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node227 = graphNodes[ihel * ndiagramgroups + 226];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node227, &node226, diagramgroup227, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node228 = graphNodes[ihel * ndiagramgroups + 227];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node228, &node227, diagramgroup228, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node229 = graphNodes[ihel * ndiagramgroups + 228];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node229, &node228, diagramgroup229, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node230 = graphNodes[ihel * ndiagramgroups + 229];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node230, &node229, diagramgroup230, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node231 = graphNodes[ihel * ndiagramgroups + 230];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node231, &node230, diagramgroup231, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node232 = graphNodes[ihel * ndiagramgroups + 231];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node232, &node231, diagramgroup232, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node233 = graphNodes[ihel * ndiagramgroups + 232];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node233, &node232, diagramgroup233, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node234 = graphNodes[ihel * ndiagramgroups + 233];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node234, &node233, diagramgroup234, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node235 = graphNodes[ihel * ndiagramgroups + 234];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node235, &node234, diagramgroup235, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node236 = graphNodes[ihel * ndiagramgroups + 235];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node236, &node235, diagramgroup236, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node237 = graphNodes[ihel * ndiagramgroups + 236];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node237, &node236, diagramgroup237, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node238 = graphNodes[ihel * ndiagramgroups + 237];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node238, &node237, diagramgroup238, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node239 = graphNodes[ihel * ndiagramgroups + 238];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node239, &node238, diagramgroup239, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node240 = graphNodes[ihel * ndiagramgroups + 239];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node240, &node239, diagramgroup240, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node241 = graphNodes[ihel * ndiagramgroups + 240];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node241, &node240, diagramgroup241, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node242 = graphNodes[ihel * ndiagramgroups + 241];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node242, &node241, diagramgroup242, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node243 = graphNodes[ihel * ndiagramgroups + 242];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node243, &node242, diagramgroup243, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node244 = graphNodes[ihel * ndiagramgroups + 243];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node244, &node243, diagramgroup244, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node245 = graphNodes[ihel * ndiagramgroups + 244];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node245, &node244, diagramgroup245, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node246 = graphNodes[ihel * ndiagramgroups + 245];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node246, &node245, diagramgroup246, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node247 = graphNodes[ihel * ndiagramgroups + 246];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node247, &node246, diagramgroup247, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node248 = graphNodes[ihel * ndiagramgroups + 247];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node248, &node247, diagramgroup248, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node249 = graphNodes[ihel * ndiagramgroups + 248];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node249, &node248, diagramgroup249, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node250 = graphNodes[ihel * ndiagramgroups + 249];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node250, &node249, diagramgroup250, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node251 = graphNodes[ihel * ndiagramgroups + 250];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node251, &node250, diagramgroup251, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node252 = graphNodes[ihel * ndiagramgroups + 251];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node252, &node251, diagramgroup252, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node253 = graphNodes[ihel * ndiagramgroups + 252];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node253, &node252, diagramgroup253, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node254 = graphNodes[ihel * ndiagramgroups + 253];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node254, &node253, diagramgroup254, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node255 = graphNodes[ihel * ndiagramgroups + 254];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node255, &node254, diagramgroup255, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node256 = graphNodes[ihel * ndiagramgroups + 255];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node256, &node255, diagramgroup256, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node257 = graphNodes[ihel * ndiagramgroups + 256];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node257, &node256, diagramgroup257, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node258 = graphNodes[ihel * ndiagramgroups + 257];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node258, &node257, diagramgroup258, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node259 = graphNodes[ihel * ndiagramgroups + 258];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node259, &node258, diagramgroup259, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node260 = graphNodes[ihel * ndiagramgroups + 259];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node260, &node259, diagramgroup260, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node261 = graphNodes[ihel * ndiagramgroups + 260];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node261, &node260, diagramgroup261, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node262 = graphNodes[ihel * ndiagramgroups + 261];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node262, &node261, diagramgroup262, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node263 = graphNodes[ihel * ndiagramgroups + 262];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node263, &node262, diagramgroup263, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node264 = graphNodes[ihel * ndiagramgroups + 263];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node264, &node263, diagramgroup264, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node265 = graphNodes[ihel * ndiagramgroups + 264];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node265, &node264, diagramgroup265, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node266 = graphNodes[ihel * ndiagramgroups + 265];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node266, &node265, diagramgroup266, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node267 = graphNodes[ihel * ndiagramgroups + 266];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node267, &node266, diagramgroup267, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node268 = graphNodes[ihel * ndiagramgroups + 267];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node268, &node267, diagramgroup268, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node269 = graphNodes[ihel * ndiagramgroups + 268];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node269, &node268, diagramgroup269, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node270 = graphNodes[ihel * ndiagramgroups + 269];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node270, &node269, diagramgroup270, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node271 = graphNodes[ihel * ndiagramgroups + 270];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node271, &node270, diagramgroup271, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node272 = graphNodes[ihel * ndiagramgroups + 271];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node272, &node271, diagramgroup272, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node273 = graphNodes[ihel * ndiagramgroups + 272];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node273, &node272, diagramgroup273, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node274 = graphNodes[ihel * ndiagramgroups + 273];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node274, &node273, diagramgroup274, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node275 = graphNodes[ihel * ndiagramgroups + 274];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node275, &node274, diagramgroup275, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node276 = graphNodes[ihel * ndiagramgroups + 275];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node276, &node275, diagramgroup276, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node277 = graphNodes[ihel * ndiagramgroups + 276];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node277, &node276, diagramgroup277, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node278 = graphNodes[ihel * ndiagramgroups + 277];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node278, &node277, diagramgroup278, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node279 = graphNodes[ihel * ndiagramgroups + 278];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node279, &node278, diagramgroup279, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node280 = graphNodes[ihel * ndiagramgroups + 279];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node280, &node279, diagramgroup280, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node281 = graphNodes[ihel * ndiagramgroups + 280];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node281, &node280, diagramgroup281, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node282 = graphNodes[ihel * ndiagramgroups + 281];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node282, &node281, diagramgroup282, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node283 = graphNodes[ihel * ndiagramgroups + 282];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node283, &node282, diagramgroup283, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node284 = graphNodes[ihel * ndiagramgroups + 283];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node284, &node283, diagramgroup284, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node285 = graphNodes[ihel * ndiagramgroups + 284];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node285, &node284, diagramgroup285, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node286 = graphNodes[ihel * ndiagramgroups + 285];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node286, &node285, diagramgroup286, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node287 = graphNodes[ihel * ndiagramgroups + 286];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node287, &node286, diagramgroup287, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node288 = graphNodes[ihel * ndiagramgroups + 287];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node288, &node287, diagramgroup288, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node289 = graphNodes[ihel * ndiagramgroups + 288];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node289, &node288, diagramgroup289, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node290 = graphNodes[ihel * ndiagramgroups + 289];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node290, &node289, diagramgroup290, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node291 = graphNodes[ihel * ndiagramgroups + 290];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node291, &node290, diagramgroup291, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node292 = graphNodes[ihel * ndiagramgroups + 291];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node292, &node291, diagramgroup292, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node293 = graphNodes[ihel * ndiagramgroups + 292];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node293, &node292, diagramgroup293, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node294 = graphNodes[ihel * ndiagramgroups + 293];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node294, &node293, diagramgroup294, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node295 = graphNodes[ihel * ndiagramgroups + 294];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node295, &node294, diagramgroup295, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node296 = graphNodes[ihel * ndiagramgroups + 295];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node296, &node295, diagramgroup296, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node297 = graphNodes[ihel * ndiagramgroups + 296];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node297, &node296, diagramgroup297, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node298 = graphNodes[ihel * ndiagramgroups + 297];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node298, &node297, diagramgroup298, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node299 = graphNodes[ihel * ndiagramgroups + 298];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node299, &node298, diagramgroup299, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node300 = graphNodes[ihel * ndiagramgroups + 299];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node300, &node299, diagramgroup300, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node301 = graphNodes[ihel * ndiagramgroups + 300];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node301, &node300, diagramgroup301, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node302 = graphNodes[ihel * ndiagramgroups + 301];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node302, &node301, diagramgroup302, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node303 = graphNodes[ihel * ndiagramgroups + 302];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node303, &node302, diagramgroup303, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node304 = graphNodes[ihel * ndiagramgroups + 303];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node304, &node303, diagramgroup304, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node305 = graphNodes[ihel * ndiagramgroups + 304];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node305, &node304, diagramgroup305, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node306 = graphNodes[ihel * ndiagramgroups + 305];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node306, &node305, diagramgroup306, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node307 = graphNodes[ihel * ndiagramgroups + 306];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node307, &node306, diagramgroup307, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node308 = graphNodes[ihel * ndiagramgroups + 307];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node308, &node307, diagramgroup308, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node309 = graphNodes[ihel * ndiagramgroups + 308];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node309, &node308, diagramgroup309, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node310 = graphNodes[ihel * ndiagramgroups + 309];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node310, &node309, diagramgroup310, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node311 = graphNodes[ihel * ndiagramgroups + 310];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node311, &node310, diagramgroup311, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node312 = graphNodes[ihel * ndiagramgroups + 311];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node312, &node311, diagramgroup312, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node313 = graphNodes[ihel * ndiagramgroups + 312];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node313, &node312, diagramgroup313, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node314 = graphNodes[ihel * ndiagramgroups + 313];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node314, &node313, diagramgroup314, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node315 = graphNodes[ihel * ndiagramgroups + 314];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node315, &node314, diagramgroup315, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node316 = graphNodes[ihel * ndiagramgroups + 315];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node316, &node315, diagramgroup316, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node317 = graphNodes[ihel * ndiagramgroups + 316];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node317, &node316, diagramgroup317, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node318 = graphNodes[ihel * ndiagramgroups + 317];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node318, &node317, diagramgroup318, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node319 = graphNodes[ihel * ndiagramgroups + 318];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node319, &node318, diagramgroup319, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node320 = graphNodes[ihel * ndiagramgroups + 319];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node320, &node319, diagramgroup320, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node321 = graphNodes[ihel * ndiagramgroups + 320];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node321, &node320, diagramgroup321, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node322 = graphNodes[ihel * ndiagramgroups + 321];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node322, &node321, diagramgroup322, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node323 = graphNodes[ihel * ndiagramgroups + 322];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node323, &node322, diagramgroup323, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node324 = graphNodes[ihel * ndiagramgroups + 323];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node324, &node323, diagramgroup324, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node325 = graphNodes[ihel * ndiagramgroups + 324];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node325, &node324, diagramgroup325, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node326 = graphNodes[ihel * ndiagramgroups + 325];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node326, &node325, diagramgroup326, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node327 = graphNodes[ihel * ndiagramgroups + 326];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node327, &node326, diagramgroup327, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node328 = graphNodes[ihel * ndiagramgroups + 327];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node328, &node327, diagramgroup328, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node329 = graphNodes[ihel * ndiagramgroups + 328];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node329, &node328, diagramgroup329, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node330 = graphNodes[ihel * ndiagramgroups + 329];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node330, &node329, diagramgroup330, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node331 = graphNodes[ihel * ndiagramgroups + 330];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node331, &node330, diagramgroup331, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node332 = graphNodes[ihel * ndiagramgroups + 331];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node332, &node331, diagramgroup332, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node333 = graphNodes[ihel * ndiagramgroups + 332];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node333, &node332, diagramgroup333, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node334 = graphNodes[ihel * ndiagramgroups + 333];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node334, &node333, diagramgroup334, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node335 = graphNodes[ihel * ndiagramgroups + 334];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node335, &node334, diagramgroup335, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node336 = graphNodes[ihel * ndiagramgroups + 335];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node336, &node335, diagramgroup336, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node337 = graphNodes[ihel * ndiagramgroups + 336];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node337, &node336, diagramgroup337, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node338 = graphNodes[ihel * ndiagramgroups + 337];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node338, &node337, diagramgroup338, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node339 = graphNodes[ihel * ndiagramgroups + 338];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node339, &node338, diagramgroup339, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node340 = graphNodes[ihel * ndiagramgroups + 339];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node340, &node339, diagramgroup340, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node341 = graphNodes[ihel * ndiagramgroups + 340];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node341, &node340, diagramgroup341, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node342 = graphNodes[ihel * ndiagramgroups + 341];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node342, &node341, diagramgroup342, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node343 = graphNodes[ihel * ndiagramgroups + 342];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node343, &node342, diagramgroup343, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node344 = graphNodes[ihel * ndiagramgroups + 343];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node344, &node343, diagramgroup344, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node345 = graphNodes[ihel * ndiagramgroups + 344];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node345, &node344, diagramgroup345, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node346 = graphNodes[ihel * ndiagramgroups + 345];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node346, &node345, diagramgroup346, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node347 = graphNodes[ihel * ndiagramgroups + 346];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node347, &node346, diagramgroup347, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node348 = graphNodes[ihel * ndiagramgroups + 347];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node348, &node347, diagramgroup348, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node349 = graphNodes[ihel * ndiagramgroups + 348];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node349, &node348, diagramgroup349, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node350 = graphNodes[ihel * ndiagramgroups + 349];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node350, &node349, diagramgroup350, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node351 = graphNodes[ihel * ndiagramgroups + 350];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node351, &node350, diagramgroup351, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node352 = graphNodes[ihel * ndiagramgroups + 351];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node352, &node351, diagramgroup352, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node353 = graphNodes[ihel * ndiagramgroups + 352];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node353, &node352, diagramgroup353, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node354 = graphNodes[ihel * ndiagramgroups + 353];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node354, &node353, diagramgroup354, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node355 = graphNodes[ihel * ndiagramgroups + 354];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node355, &node354, diagramgroup355, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node356 = graphNodes[ihel * ndiagramgroups + 355];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node356, &node355, diagramgroup356, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node357 = graphNodes[ihel * ndiagramgroups + 356];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node357, &node356, diagramgroup357, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node358 = graphNodes[ihel * ndiagramgroups + 357];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node358, &node357, diagramgroup358, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node359 = graphNodes[ihel * ndiagramgroups + 358];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node359, &node358, diagramgroup359, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node360 = graphNodes[ihel * ndiagramgroups + 359];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node360, &node359, diagramgroup360, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node361 = graphNodes[ihel * ndiagramgroups + 360];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node361, &node360, diagramgroup361, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node362 = graphNodes[ihel * ndiagramgroups + 361];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node362, &node361, diagramgroup362, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node363 = graphNodes[ihel * ndiagramgroups + 362];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node363, &node362, diagramgroup363, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node364 = graphNodes[ihel * ndiagramgroups + 363];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node364, &node363, diagramgroup364, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node365 = graphNodes[ihel * ndiagramgroups + 364];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node365, &node364, diagramgroup365, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node366 = graphNodes[ihel * ndiagramgroups + 365];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node366, &node365, diagramgroup366, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node367 = graphNodes[ihel * ndiagramgroups + 366];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node367, &node366, diagramgroup367, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node368 = graphNodes[ihel * ndiagramgroups + 367];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node368, &node367, diagramgroup368, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node369 = graphNodes[ihel * ndiagramgroups + 368];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node369, &node368, diagramgroup369, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node370 = graphNodes[ihel * ndiagramgroups + 369];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node370, &node369, diagramgroup370, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node371 = graphNodes[ihel * ndiagramgroups + 370];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node371, &node370, diagramgroup371, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node372 = graphNodes[ihel * ndiagramgroups + 371];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node372, &node371, diagramgroup372, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node373 = graphNodes[ihel * ndiagramgroups + 372];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node373, &node372, diagramgroup373, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node374 = graphNodes[ihel * ndiagramgroups + 373];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node374, &node373, diagramgroup374, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node375 = graphNodes[ihel * ndiagramgroups + 374];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node375, &node374, diagramgroup375, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node376 = graphNodes[ihel * ndiagramgroups + 375];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node376, &node375, diagramgroup376, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node377 = graphNodes[ihel * ndiagramgroups + 376];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node377, &node376, diagramgroup377, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node378 = graphNodes[ihel * ndiagramgroups + 377];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node378, &node377, diagramgroup378, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node379 = graphNodes[ihel * ndiagramgroups + 378];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node379, &node378, diagramgroup379, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node380 = graphNodes[ihel * ndiagramgroups + 379];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node380, &node379, diagramgroup380, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node381 = graphNodes[ihel * ndiagramgroups + 380];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node381, &node380, diagramgroup381, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node382 = graphNodes[ihel * ndiagramgroups + 381];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node382, &node381, diagramgroup382, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node383 = graphNodes[ihel * ndiagramgroups + 382];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node383, &node382, diagramgroup383, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node384 = graphNodes[ihel * ndiagramgroups + 383];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node384, &node383, diagramgroup384, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node385 = graphNodes[ihel * ndiagramgroups + 384];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node385, &node384, diagramgroup385, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node386 = graphNodes[ihel * ndiagramgroups + 385];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node386, &node385, diagramgroup386, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node387 = graphNodes[ihel * ndiagramgroups + 386];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node387, &node386, diagramgroup387, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node388 = graphNodes[ihel * ndiagramgroups + 387];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node388, &node387, diagramgroup388, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node389 = graphNodes[ihel * ndiagramgroups + 388];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node389, &node388, diagramgroup389, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node390 = graphNodes[ihel * ndiagramgroups + 389];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node390, &node389, diagramgroup390, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node391 = graphNodes[ihel * ndiagramgroups + 390];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node391, &node390, diagramgroup391, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node392 = graphNodes[ihel * ndiagramgroups + 391];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node392, &node391, diagramgroup392, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node393 = graphNodes[ihel * ndiagramgroups + 392];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node393, &node392, diagramgroup393, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node394 = graphNodes[ihel * ndiagramgroups + 393];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node394, &node393, diagramgroup394, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node395 = graphNodes[ihel * ndiagramgroups + 394];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node395, &node394, diagramgroup395, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node396 = graphNodes[ihel * ndiagramgroups + 395];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node396, &node395, diagramgroup396, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node397 = graphNodes[ihel * ndiagramgroups + 396];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node397, &node396, diagramgroup397, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node398 = graphNodes[ihel * ndiagramgroups + 397];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node398, &node397, diagramgroup398, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node399 = graphNodes[ihel * ndiagramgroups + 398];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node399, &node398, diagramgroup399, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node400 = graphNodes[ihel * ndiagramgroups + 399];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node400, &node399, diagramgroup400, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node401 = graphNodes[ihel * ndiagramgroups + 400];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node401, &node400, diagramgroup401, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node402 = graphNodes[ihel * ndiagramgroups + 401];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node402, &node401, diagramgroup402, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node403 = graphNodes[ihel * ndiagramgroups + 402];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node403, &node402, diagramgroup403, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node404 = graphNodes[ihel * ndiagramgroups + 403];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node404, &node403, diagramgroup404, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node405 = graphNodes[ihel * ndiagramgroups + 404];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node405, &node404, diagramgroup405, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node406 = graphNodes[ihel * ndiagramgroups + 405];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node406, &node405, diagramgroup406, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node407 = graphNodes[ihel * ndiagramgroups + 406];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node407, &node406, diagramgroup407, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node408 = graphNodes[ihel * ndiagramgroups + 407];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node408, &node407, diagramgroup408, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node409 = graphNodes[ihel * ndiagramgroups + 408];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node409, &node408, diagramgroup409, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node410 = graphNodes[ihel * ndiagramgroups + 409];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node410, &node409, diagramgroup410, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node411 = graphNodes[ihel * ndiagramgroups + 410];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node411, &node410, diagramgroup411, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node412 = graphNodes[ihel * ndiagramgroups + 411];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node412, &node411, diagramgroup412, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node413 = graphNodes[ihel * ndiagramgroups + 412];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node413, &node412, diagramgroup413, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node414 = graphNodes[ihel * ndiagramgroups + 413];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node414, &node413, diagramgroup414, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node415 = graphNodes[ihel * ndiagramgroups + 414];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node415, &node414, diagramgroup415, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node416 = graphNodes[ihel * ndiagramgroups + 415];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node416, &node415, diagramgroup416, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node417 = graphNodes[ihel * ndiagramgroups + 416];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node417, &node416, diagramgroup417, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node418 = graphNodes[ihel * ndiagramgroups + 417];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node418, &node417, diagramgroup418, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node419 = graphNodes[ihel * ndiagramgroups + 418];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node419, &node418, diagramgroup419, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node420 = graphNodes[ihel * ndiagramgroups + 419];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node420, &node419, diagramgroup420, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node421 = graphNodes[ihel * ndiagramgroups + 420];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node421, &node420, diagramgroup421, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node422 = graphNodes[ihel * ndiagramgroups + 421];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node422, &node421, diagramgroup422, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node423 = graphNodes[ihel * ndiagramgroups + 422];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node423, &node422, diagramgroup423, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node424 = graphNodes[ihel * ndiagramgroups + 423];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node424, &node423, diagramgroup424, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node425 = graphNodes[ihel * ndiagramgroups + 424];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node425, &node424, diagramgroup425, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node426 = graphNodes[ihel * ndiagramgroups + 425];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node426, &node425, diagramgroup426, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node427 = graphNodes[ihel * ndiagramgroups + 426];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node427, &node426, diagramgroup427, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node428 = graphNodes[ihel * ndiagramgroups + 427];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node428, &node427, diagramgroup428, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node429 = graphNodes[ihel * ndiagramgroups + 428];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node429, &node428, diagramgroup429, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node430 = graphNodes[ihel * ndiagramgroups + 429];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node430, &node429, diagramgroup430, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node431 = graphNodes[ihel * ndiagramgroups + 430];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node431, &node430, diagramgroup431, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node432 = graphNodes[ihel * ndiagramgroups + 431];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node432, &node431, diagramgroup432, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node433 = graphNodes[ihel * ndiagramgroups + 432];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node433, &node432, diagramgroup433, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node434 = graphNodes[ihel * ndiagramgroups + 433];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node434, &node433, diagramgroup434, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node435 = graphNodes[ihel * ndiagramgroups + 434];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node435, &node434, diagramgroup435, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node436 = graphNodes[ihel * ndiagramgroups + 435];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node436, &node435, diagramgroup436, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node437 = graphNodes[ihel * ndiagramgroups + 436];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node437, &node436, diagramgroup437, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node438 = graphNodes[ihel * ndiagramgroups + 437];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node438, &node437, diagramgroup438, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node439 = graphNodes[ihel * ndiagramgroups + 438];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node439, &node438, diagramgroup439, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node440 = graphNodes[ihel * ndiagramgroups + 439];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node440, &node439, diagramgroup440, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node441 = graphNodes[ihel * ndiagramgroups + 440];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node441, &node440, diagramgroup441, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node442 = graphNodes[ihel * ndiagramgroups + 441];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node442, &node441, diagramgroup442, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node443 = graphNodes[ihel * ndiagramgroups + 442];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node443, &node442, diagramgroup443, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node444 = graphNodes[ihel * ndiagramgroups + 443];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node444, &node443, diagramgroup444, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node445 = graphNodes[ihel * ndiagramgroups + 444];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node445, &node444, diagramgroup445, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node446 = graphNodes[ihel * ndiagramgroups + 445];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node446, &node445, diagramgroup446, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node447 = graphNodes[ihel * ndiagramgroups + 446];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node447, &node446, diagramgroup447, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node448 = graphNodes[ihel * ndiagramgroups + 447];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node448, &node447, diagramgroup448, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node449 = graphNodes[ihel * ndiagramgroups + 448];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node449, &node448, diagramgroup449, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node450 = graphNodes[ihel * ndiagramgroups + 449];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node450, &node449, diagramgroup450, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node451 = graphNodes[ihel * ndiagramgroups + 450];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node451, &node450, diagramgroup451, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node452 = graphNodes[ihel * ndiagramgroups + 451];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node452, &node451, diagramgroup452, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node453 = graphNodes[ihel * ndiagramgroups + 452];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node453, &node452, diagramgroup453, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node454 = graphNodes[ihel * ndiagramgroups + 453];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node454, &node453, diagramgroup454, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node455 = graphNodes[ihel * ndiagramgroups + 454];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node455, &node454, diagramgroup455, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node456 = graphNodes[ihel * ndiagramgroups + 455];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node456, &node455, diagramgroup456, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node457 = graphNodes[ihel * ndiagramgroups + 456];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node457, &node456, diagramgroup457, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node458 = graphNodes[ihel * ndiagramgroups + 457];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node458, &node457, diagramgroup458, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node459 = graphNodes[ihel * ndiagramgroups + 458];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node459, &node458, diagramgroup459, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node460 = graphNodes[ihel * ndiagramgroups + 459];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node460, &node459, diagramgroup460, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node461 = graphNodes[ihel * ndiagramgroups + 460];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node461, &node460, diagramgroup461, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node462 = graphNodes[ihel * ndiagramgroups + 461];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node462, &node461, diagramgroup462, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node463 = graphNodes[ihel * ndiagramgroups + 462];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node463, &node462, diagramgroup463, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node464 = graphNodes[ihel * ndiagramgroups + 463];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node464, &node463, diagramgroup464, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node465 = graphNodes[ihel * ndiagramgroups + 464];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node465, &node464, diagramgroup465, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node466 = graphNodes[ihel * ndiagramgroups + 465];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node466, &node465, diagramgroup466, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node467 = graphNodes[ihel * ndiagramgroups + 466];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node467, &node466, diagramgroup467, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node468 = graphNodes[ihel * ndiagramgroups + 467];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node468, &node467, diagramgroup468, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node469 = graphNodes[ihel * ndiagramgroups + 468];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node469, &node468, diagramgroup469, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node470 = graphNodes[ihel * ndiagramgroups + 469];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node470, &node469, diagramgroup470, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node471 = graphNodes[ihel * ndiagramgroups + 470];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node471, &node470, diagramgroup471, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node472 = graphNodes[ihel * ndiagramgroups + 471];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node472, &node471, diagramgroup472, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node473 = graphNodes[ihel * ndiagramgroups + 472];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node473, &node472, diagramgroup473, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node474 = graphNodes[ihel * ndiagramgroups + 473];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node474, &node473, diagramgroup474, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node475 = graphNodes[ihel * ndiagramgroups + 474];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node475, &node474, diagramgroup475, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node476 = graphNodes[ihel * ndiagramgroups + 475];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node476, &node475, diagramgroup476, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node477 = graphNodes[ihel * ndiagramgroups + 476];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node477, &node476, diagramgroup477, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node478 = graphNodes[ihel * ndiagramgroups + 477];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node478, &node477, diagramgroup478, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node479 = graphNodes[ihel * ndiagramgroups + 478];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node479, &node478, diagramgroup479, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node480 = graphNodes[ihel * ndiagramgroups + 479];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node480, &node479, diagramgroup480, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node481 = graphNodes[ihel * ndiagramgroups + 480];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node481, &node480, diagramgroup481, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node482 = graphNodes[ihel * ndiagramgroups + 481];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node482, &node481, diagramgroup482, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node483 = graphNodes[ihel * ndiagramgroups + 482];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node483, &node482, diagramgroup483, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node484 = graphNodes[ihel * ndiagramgroups + 483];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node484, &node483, diagramgroup484, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node485 = graphNodes[ihel * ndiagramgroups + 484];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node485, &node484, diagramgroup485, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node486 = graphNodes[ihel * ndiagramgroups + 485];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node486, &node485, diagramgroup486, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node487 = graphNodes[ihel * ndiagramgroups + 486];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node487, &node486, diagramgroup487, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node488 = graphNodes[ihel * ndiagramgroups + 487];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node488, &node487, diagramgroup488, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node489 = graphNodes[ihel * ndiagramgroups + 488];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node489, &node488, diagramgroup489, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node490 = graphNodes[ihel * ndiagramgroups + 489];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node490, &node489, diagramgroup490, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node491 = graphNodes[ihel * ndiagramgroups + 490];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node491, &node490, diagramgroup491, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node492 = graphNodes[ihel * ndiagramgroups + 491];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node492, &node491, diagramgroup492, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node493 = graphNodes[ihel * ndiagramgroups + 492];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node493, &node492, diagramgroup493, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node494 = graphNodes[ihel * ndiagramgroups + 493];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node494, &node493, diagramgroup494, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node495 = graphNodes[ihel * ndiagramgroups + 494];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node495, &node494, diagramgroup495, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node496 = graphNodes[ihel * ndiagramgroups + 495];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node496, &node495, diagramgroup496, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node497 = graphNodes[ihel * ndiagramgroups + 496];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node497, &node496, diagramgroup497, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node498 = graphNodes[ihel * ndiagramgroups + 497];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node498, &node497, diagramgroup498, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node499 = graphNodes[ihel * ndiagramgroups + 498];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node499, &node498, diagramgroup499, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node500 = graphNodes[ihel * ndiagramgroups + 499];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node500, &node499, diagramgroup500, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node501 = graphNodes[ihel * ndiagramgroups + 500];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node501, &node500, diagramgroup501, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node502 = graphNodes[ihel * ndiagramgroups + 501];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node502, &node501, diagramgroup502, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node503 = graphNodes[ihel * ndiagramgroups + 502];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node503, &node502, diagramgroup503, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node504 = graphNodes[ihel * ndiagramgroups + 503];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node504, &node503, diagramgroup504, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node505 = graphNodes[ihel * ndiagramgroups + 504];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node505, &node504, diagramgroup505, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node506 = graphNodes[ihel * ndiagramgroups + 505];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node506, &node505, diagramgroup506, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node507 = graphNodes[ihel * ndiagramgroups + 506];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node507, &node506, diagramgroup507, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node508 = graphNodes[ihel * ndiagramgroups + 507];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node508, &node507, diagramgroup508, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node509 = graphNodes[ihel * ndiagramgroups + 508];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node509, &node508, diagramgroup509, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node510 = graphNodes[ihel * ndiagramgroups + 509];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node510, &node509, diagramgroup510, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node511 = graphNodes[ihel * ndiagramgroups + 510];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node511, &node510, diagramgroup511, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node512 = graphNodes[ihel * ndiagramgroups + 511];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node512, &node511, diagramgroup512, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node513 = graphNodes[ihel * ndiagramgroups + 512];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node513, &node512, diagramgroup513, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node514 = graphNodes[ihel * ndiagramgroups + 513];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node514, &node513, diagramgroup514, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node515 = graphNodes[ihel * ndiagramgroups + 514];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node515, &node514, diagramgroup515, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node516 = graphNodes[ihel * ndiagramgroups + 515];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node516, &node515, diagramgroup516, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node517 = graphNodes[ihel * ndiagramgroups + 516];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node517, &node516, diagramgroup517, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node518 = graphNodes[ihel * ndiagramgroups + 517];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node518, &node517, diagramgroup518, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node519 = graphNodes[ihel * ndiagramgroups + 518];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node519, &node518, diagramgroup519, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node520 = graphNodes[ihel * ndiagramgroups + 519];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node520, &node519, diagramgroup520, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node521 = graphNodes[ihel * ndiagramgroups + 520];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node521, &node520, diagramgroup521, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node522 = graphNodes[ihel * ndiagramgroups + 521];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node522, &node521, diagramgroup522, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node523 = graphNodes[ihel * ndiagramgroups + 522];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node523, &node522, diagramgroup523, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node524 = graphNodes[ihel * ndiagramgroups + 523];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node524, &node523, diagramgroup524, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node525 = graphNodes[ihel * ndiagramgroups + 524];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node525, &node524, diagramgroup525, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node526 = graphNodes[ihel * ndiagramgroups + 525];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node526, &node525, diagramgroup526, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node527 = graphNodes[ihel * ndiagramgroups + 526];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node527, &node526, diagramgroup527, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node528 = graphNodes[ihel * ndiagramgroups + 527];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node528, &node527, diagramgroup528, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node529 = graphNodes[ihel * ndiagramgroups + 528];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node529, &node528, diagramgroup529, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node530 = graphNodes[ihel * ndiagramgroups + 529];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node530, &node529, diagramgroup530, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node531 = graphNodes[ihel * ndiagramgroups + 530];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node531, &node530, diagramgroup531, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node532 = graphNodes[ihel * ndiagramgroups + 531];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node532, &node531, diagramgroup532, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node533 = graphNodes[ihel * ndiagramgroups + 532];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node533, &node532, diagramgroup533, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node534 = graphNodes[ihel * ndiagramgroups + 533];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node534, &node533, diagramgroup534, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node535 = graphNodes[ihel * ndiagramgroups + 534];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node535, &node534, diagramgroup535, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node536 = graphNodes[ihel * ndiagramgroups + 535];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node536, &node535, diagramgroup536, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node537 = graphNodes[ihel * ndiagramgroups + 536];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node537, &node536, diagramgroup537, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node538 = graphNodes[ihel * ndiagramgroups + 537];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node538, &node537, diagramgroup538, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node539 = graphNodes[ihel * ndiagramgroups + 538];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node539, &node538, diagramgroup539, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node540 = graphNodes[ihel * ndiagramgroups + 539];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node540, &node539, diagramgroup540, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node541 = graphNodes[ihel * ndiagramgroups + 540];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node541, &node540, diagramgroup541, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node542 = graphNodes[ihel * ndiagramgroups + 541];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node542, &node541, diagramgroup542, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node543 = graphNodes[ihel * ndiagramgroups + 542];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node543, &node542, diagramgroup543, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node544 = graphNodes[ihel * ndiagramgroups + 543];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node544, &node543, diagramgroup544, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node545 = graphNodes[ihel * ndiagramgroups + 544];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node545, &node544, diagramgroup545, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node546 = graphNodes[ihel * ndiagramgroups + 545];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node546, &node545, diagramgroup546, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node547 = graphNodes[ihel * ndiagramgroups + 546];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node547, &node546, diagramgroup547, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node548 = graphNodes[ihel * ndiagramgroups + 547];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node548, &node547, diagramgroup548, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node549 = graphNodes[ihel * ndiagramgroups + 548];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node549, &node548, diagramgroup549, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node550 = graphNodes[ihel * ndiagramgroups + 549];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node550, &node549, diagramgroup550, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node551 = graphNodes[ihel * ndiagramgroups + 550];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node551, &node550, diagramgroup551, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node552 = graphNodes[ihel * ndiagramgroups + 551];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node552, &node551, diagramgroup552, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node553 = graphNodes[ihel * ndiagramgroups + 552];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node553, &node552, diagramgroup553, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node554 = graphNodes[ihel * ndiagramgroups + 553];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node554, &node553, diagramgroup554, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node555 = graphNodes[ihel * ndiagramgroups + 554];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node555, &node554, diagramgroup555, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node556 = graphNodes[ihel * ndiagramgroups + 555];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node556, &node555, diagramgroup556, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node557 = graphNodes[ihel * ndiagramgroups + 556];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node557, &node556, diagramgroup557, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node558 = graphNodes[ihel * ndiagramgroups + 557];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node558, &node557, diagramgroup558, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node559 = graphNodes[ihel * ndiagramgroups + 558];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node559, &node558, diagramgroup559, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node560 = graphNodes[ihel * ndiagramgroups + 559];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node560, &node559, diagramgroup560, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node561 = graphNodes[ihel * ndiagramgroups + 560];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node561, &node560, diagramgroup561, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node562 = graphNodes[ihel * ndiagramgroups + 561];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node562, &node561, diagramgroup562, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node563 = graphNodes[ihel * ndiagramgroups + 562];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node563, &node562, diagramgroup563, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node564 = graphNodes[ihel * ndiagramgroups + 563];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node564, &node563, diagramgroup564, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node565 = graphNodes[ihel * ndiagramgroups + 564];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node565, &node564, diagramgroup565, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node566 = graphNodes[ihel * ndiagramgroups + 565];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node566, &node565, diagramgroup566, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node567 = graphNodes[ihel * ndiagramgroups + 566];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node567, &node566, diagramgroup567, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node568 = graphNodes[ihel * ndiagramgroups + 567];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node568, &node567, diagramgroup568, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node569 = graphNodes[ihel * ndiagramgroups + 568];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node569, &node568, diagramgroup569, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node570 = graphNodes[ihel * ndiagramgroups + 569];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node570, &node569, diagramgroup570, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node571 = graphNodes[ihel * ndiagramgroups + 570];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node571, &node570, diagramgroup571, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node572 = graphNodes[ihel * ndiagramgroups + 571];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node572, &node571, diagramgroup572, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node573 = graphNodes[ihel * ndiagramgroups + 572];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node573, &node572, diagramgroup573, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node574 = graphNodes[ihel * ndiagramgroups + 573];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node574, &node573, diagramgroup574, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node575 = graphNodes[ihel * ndiagramgroups + 574];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node575, &node574, diagramgroup575, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node576 = graphNodes[ihel * ndiagramgroups + 575];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node576, &node575, diagramgroup576, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node577 = graphNodes[ihel * ndiagramgroups + 576];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node577, &node576, diagramgroup577, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node578 = graphNodes[ihel * ndiagramgroups + 577];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node578, &node577, diagramgroup578, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node579 = graphNodes[ihel * ndiagramgroups + 578];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node579, &node578, diagramgroup579, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node580 = graphNodes[ihel * ndiagramgroups + 579];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node580, &node579, diagramgroup580, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node581 = graphNodes[ihel * ndiagramgroups + 580];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node581, &node580, diagramgroup581, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node582 = graphNodes[ihel * ndiagramgroups + 581];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node582, &node581, diagramgroup582, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node583 = graphNodes[ihel * ndiagramgroups + 582];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node583, &node582, diagramgroup583, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node584 = graphNodes[ihel * ndiagramgroups + 583];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node584, &node583, diagramgroup584, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node585 = graphNodes[ihel * ndiagramgroups + 584];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node585, &node584, diagramgroup585, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node586 = graphNodes[ihel * ndiagramgroups + 585];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node586, &node585, diagramgroup586, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node587 = graphNodes[ihel * ndiagramgroups + 586];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node587, &node586, diagramgroup587, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node588 = graphNodes[ihel * ndiagramgroups + 587];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node588, &node587, diagramgroup588, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node589 = graphNodes[ihel * ndiagramgroups + 588];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node589, &node588, diagramgroup589, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node590 = graphNodes[ihel * ndiagramgroups + 589];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node590, &node589, diagramgroup590, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node591 = graphNodes[ihel * ndiagramgroups + 590];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node591, &node590, diagramgroup591, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node592 = graphNodes[ihel * ndiagramgroups + 591];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node592, &node591, diagramgroup592, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node593 = graphNodes[ihel * ndiagramgroups + 592];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node593, &node592, diagramgroup593, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node594 = graphNodes[ihel * ndiagramgroups + 593];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node594, &node593, diagramgroup594, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node595 = graphNodes[ihel * ndiagramgroups + 594];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node595, &node594, diagramgroup595, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node596 = graphNodes[ihel * ndiagramgroups + 595];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node596, &node595, diagramgroup596, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node597 = graphNodes[ihel * ndiagramgroups + 596];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node597, &node596, diagramgroup597, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node598 = graphNodes[ihel * ndiagramgroups + 597];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node598, &node597, diagramgroup598, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node599 = graphNodes[ihel * ndiagramgroups + 598];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node599, &node598, diagramgroup599, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node600 = graphNodes[ihel * ndiagramgroups + 599];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node600, &node599, diagramgroup600, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node601 = graphNodes[ihel * ndiagramgroups + 600];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node601, &node600, diagramgroup601, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node602 = graphNodes[ihel * ndiagramgroups + 601];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node602, &node601, diagramgroup602, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node603 = graphNodes[ihel * ndiagramgroups + 602];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node603, &node602, diagramgroup603, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node604 = graphNodes[ihel * ndiagramgroups + 603];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node604, &node603, diagramgroup604, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node605 = graphNodes[ihel * ndiagramgroups + 604];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node605, &node604, diagramgroup605, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node606 = graphNodes[ihel * ndiagramgroups + 605];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node606, &node605, diagramgroup606, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node607 = graphNodes[ihel * ndiagramgroups + 606];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node607, &node606, diagramgroup607, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node608 = graphNodes[ihel * ndiagramgroups + 607];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node608, &node607, diagramgroup608, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node609 = graphNodes[ihel * ndiagramgroups + 608];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node609, &node608, diagramgroup609, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node610 = graphNodes[ihel * ndiagramgroups + 609];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node610, &node609, diagramgroup610, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node611 = graphNodes[ihel * ndiagramgroups + 610];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node611, &node610, diagramgroup611, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node612 = graphNodes[ihel * ndiagramgroups + 611];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node612, &node611, diagramgroup612, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node613 = graphNodes[ihel * ndiagramgroups + 612];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node613, &node612, diagramgroup613, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node614 = graphNodes[ihel * ndiagramgroups + 613];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node614, &node613, diagramgroup614, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node615 = graphNodes[ihel * ndiagramgroups + 614];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node615, &node614, diagramgroup615, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node616 = graphNodes[ihel * ndiagramgroups + 615];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node616, &node615, diagramgroup616, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node617 = graphNodes[ihel * ndiagramgroups + 616];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node617, &node616, diagramgroup617, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node618 = graphNodes[ihel * ndiagramgroups + 617];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node618, &node617, diagramgroup618, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node619 = graphNodes[ihel * ndiagramgroups + 618];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node619, &node618, diagramgroup619, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node620 = graphNodes[ihel * ndiagramgroups + 619];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node620, &node619, diagramgroup620, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node621 = graphNodes[ihel * ndiagramgroups + 620];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node621, &node620, diagramgroup621, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node622 = graphNodes[ihel * ndiagramgroups + 621];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node622, &node621, diagramgroup622, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node623 = graphNodes[ihel * ndiagramgroups + 622];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node623, &node622, diagramgroup623, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node624 = graphNodes[ihel * ndiagramgroups + 623];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node624, &node623, diagramgroup624, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node625 = graphNodes[ihel * ndiagramgroups + 624];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node625, &node624, diagramgroup625, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node626 = graphNodes[ihel * ndiagramgroups + 625];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node626, &node625, diagramgroup626, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node627 = graphNodes[ihel * ndiagramgroups + 626];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node627, &node626, diagramgroup627, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node628 = graphNodes[ihel * ndiagramgroups + 627];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node628, &node627, diagramgroup628, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node629 = graphNodes[ihel * ndiagramgroups + 628];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node629, &node628, diagramgroup629, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node630 = graphNodes[ihel * ndiagramgroups + 629];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node630, &node629, diagramgroup630, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node631 = graphNodes[ihel * ndiagramgroups + 630];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node631, &node630, diagramgroup631, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node632 = graphNodes[ihel * ndiagramgroups + 631];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node632, &node631, diagramgroup632, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node633 = graphNodes[ihel * ndiagramgroups + 632];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node633, &node632, diagramgroup633, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node634 = graphNodes[ihel * ndiagramgroups + 633];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node634, &node633, diagramgroup634, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node635 = graphNodes[ihel * ndiagramgroups + 634];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node635, &node634, diagramgroup635, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node636 = graphNodes[ihel * ndiagramgroups + 635];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node636, &node635, diagramgroup636, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node637 = graphNodes[ihel * ndiagramgroups + 636];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node637, &node636, diagramgroup637, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node638 = graphNodes[ihel * ndiagramgroups + 637];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node638, &node637, diagramgroup638, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node639 = graphNodes[ihel * ndiagramgroups + 638];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node639, &node638, diagramgroup639, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node640 = graphNodes[ihel * ndiagramgroups + 639];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node640, &node639, diagramgroup640, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node641 = graphNodes[ihel * ndiagramgroups + 640];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node641, &node640, diagramgroup641, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node642 = graphNodes[ihel * ndiagramgroups + 641];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node642, &node641, diagramgroup642, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node643 = graphNodes[ihel * ndiagramgroups + 642];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node643, &node642, diagramgroup643, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node644 = graphNodes[ihel * ndiagramgroups + 643];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node644, &node643, diagramgroup644, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node645 = graphNodes[ihel * ndiagramgroups + 644];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node645, &node644, diagramgroup645, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node646 = graphNodes[ihel * ndiagramgroups + 645];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node646, &node645, diagramgroup646, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node647 = graphNodes[ihel * ndiagramgroups + 646];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node647, &node646, diagramgroup647, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node648 = graphNodes[ihel * ndiagramgroups + 647];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node648, &node647, diagramgroup648, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node649 = graphNodes[ihel * ndiagramgroups + 648];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node649, &node648, diagramgroup649, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node650 = graphNodes[ihel * ndiagramgroups + 649];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node650, &node649, diagramgroup650, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node651 = graphNodes[ihel * ndiagramgroups + 650];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node651, &node650, diagramgroup651, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node652 = graphNodes[ihel * ndiagramgroups + 651];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node652, &node651, diagramgroup652, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node653 = graphNodes[ihel * ndiagramgroups + 652];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node653, &node652, diagramgroup653, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node654 = graphNodes[ihel * ndiagramgroups + 653];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node654, &node653, diagramgroup654, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node655 = graphNodes[ihel * ndiagramgroups + 654];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node655, &node654, diagramgroup655, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node656 = graphNodes[ihel * ndiagramgroups + 655];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node656, &node655, diagramgroup656, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node657 = graphNodes[ihel * ndiagramgroups + 656];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node657, &node656, diagramgroup657, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node658 = graphNodes[ihel * ndiagramgroups + 657];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node658, &node657, diagramgroup658, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node659 = graphNodes[ihel * ndiagramgroups + 658];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node659, &node658, diagramgroup659, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node660 = graphNodes[ihel * ndiagramgroups + 659];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node660, &node659, diagramgroup660, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node661 = graphNodes[ihel * ndiagramgroups + 660];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node661, &node660, diagramgroup661, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node662 = graphNodes[ihel * ndiagramgroups + 661];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node662, &node661, diagramgroup662, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node663 = graphNodes[ihel * ndiagramgroups + 662];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node663, &node662, diagramgroup663, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node664 = graphNodes[ihel * ndiagramgroups + 663];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node664, &node663, diagramgroup664, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node665 = graphNodes[ihel * ndiagramgroups + 664];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node665, &node664, diagramgroup665, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node666 = graphNodes[ihel * ndiagramgroups + 665];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node666, &node665, diagramgroup666, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node667 = graphNodes[ihel * ndiagramgroups + 666];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node667, &node666, diagramgroup667, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node668 = graphNodes[ihel * ndiagramgroups + 667];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node668, &node667, diagramgroup668, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node669 = graphNodes[ihel * ndiagramgroups + 668];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node669, &node668, diagramgroup669, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node670 = graphNodes[ihel * ndiagramgroups + 669];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node670, &node669, diagramgroup670, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node671 = graphNodes[ihel * ndiagramgroups + 670];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node671, &node670, diagramgroup671, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node672 = graphNodes[ihel * ndiagramgroups + 671];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node672, &node671, diagramgroup672, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node673 = graphNodes[ihel * ndiagramgroups + 672];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node673, &node672, diagramgroup673, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node674 = graphNodes[ihel * ndiagramgroups + 673];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node674, &node673, diagramgroup674, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node675 = graphNodes[ihel * ndiagramgroups + 674];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node675, &node674, diagramgroup675, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node676 = graphNodes[ihel * ndiagramgroups + 675];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node676, &node675, diagramgroup676, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node677 = graphNodes[ihel * ndiagramgroups + 676];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node677, &node676, diagramgroup677, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node678 = graphNodes[ihel * ndiagramgroups + 677];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node678, &node677, diagramgroup678, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node679 = graphNodes[ihel * ndiagramgroups + 678];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node679, &node678, diagramgroup679, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node680 = graphNodes[ihel * ndiagramgroups + 679];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node680, &node679, diagramgroup680, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node681 = graphNodes[ihel * ndiagramgroups + 680];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node681, &node680, diagramgroup681, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node682 = graphNodes[ihel * ndiagramgroups + 681];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node682, &node681, diagramgroup682, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node683 = graphNodes[ihel * ndiagramgroups + 682];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node683, &node682, diagramgroup683, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node684 = graphNodes[ihel * ndiagramgroups + 683];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node684, &node683, diagramgroup684, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node685 = graphNodes[ihel * ndiagramgroups + 684];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node685, &node684, diagramgroup685, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node686 = graphNodes[ihel * ndiagramgroups + 685];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node686, &node685, diagramgroup686, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node687 = graphNodes[ihel * ndiagramgroups + 686];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node687, &node686, diagramgroup687, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node688 = graphNodes[ihel * ndiagramgroups + 687];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node688, &node687, diagramgroup688, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node689 = graphNodes[ihel * ndiagramgroups + 688];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node689, &node688, diagramgroup689, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node690 = graphNodes[ihel * ndiagramgroups + 689];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node690, &node689, diagramgroup690, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node691 = graphNodes[ihel * ndiagramgroups + 690];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node691, &node690, diagramgroup691, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node692 = graphNodes[ihel * ndiagramgroups + 691];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node692, &node691, diagramgroup692, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node693 = graphNodes[ihel * ndiagramgroups + 692];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node693, &node692, diagramgroup693, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node694 = graphNodes[ihel * ndiagramgroups + 693];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node694, &node693, diagramgroup694, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node695 = graphNodes[ihel * ndiagramgroups + 694];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node695, &node694, diagramgroup695, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node696 = graphNodes[ihel * ndiagramgroups + 695];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node696, &node695, diagramgroup696, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node697 = graphNodes[ihel * ndiagramgroups + 696];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node697, &node696, diagramgroup697, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node698 = graphNodes[ihel * ndiagramgroups + 697];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node698, &node697, diagramgroup698, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node699 = graphNodes[ihel * ndiagramgroups + 698];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node699, &node698, diagramgroup699, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node700 = graphNodes[ihel * ndiagramgroups + 699];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node700, &node699, diagramgroup700, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node701 = graphNodes[ihel * ndiagramgroups + 700];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node701, &node700, diagramgroup701, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node702 = graphNodes[ihel * ndiagramgroups + 701];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node702, &node701, diagramgroup702, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node703 = graphNodes[ihel * ndiagramgroups + 702];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node703, &node702, diagramgroup703, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node704 = graphNodes[ihel * ndiagramgroups + 703];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node704, &node703, diagramgroup704, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node705 = graphNodes[ihel * ndiagramgroups + 704];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node705, &node704, diagramgroup705, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node706 = graphNodes[ihel * ndiagramgroups + 705];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node706, &node705, diagramgroup706, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node707 = graphNodes[ihel * ndiagramgroups + 706];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node707, &node706, diagramgroup707, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node708 = graphNodes[ihel * ndiagramgroups + 707];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node708, &node707, diagramgroup708, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node709 = graphNodes[ihel * ndiagramgroups + 708];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node709, &node708, diagramgroup709, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node710 = graphNodes[ihel * ndiagramgroups + 709];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node710, &node709, diagramgroup710, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node711 = graphNodes[ihel * ndiagramgroups + 710];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node711, &node710, diagramgroup711, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node712 = graphNodes[ihel * ndiagramgroups + 711];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node712, &node711, diagramgroup712, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node713 = graphNodes[ihel * ndiagramgroups + 712];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node713, &node712, diagramgroup713, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node714 = graphNodes[ihel * ndiagramgroups + 713];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node714, &node713, diagramgroup714, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node715 = graphNodes[ihel * ndiagramgroups + 714];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node715, &node714, diagramgroup715, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node716 = graphNodes[ihel * ndiagramgroups + 715];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node716, &node715, diagramgroup716, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node717 = graphNodes[ihel * ndiagramgroups + 716];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node717, &node716, diagramgroup717, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node718 = graphNodes[ihel * ndiagramgroups + 717];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node718, &node717, diagramgroup718, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node719 = graphNodes[ihel * ndiagramgroups + 718];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node719, &node718, diagramgroup719, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node720 = graphNodes[ihel * ndiagramgroups + 719];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node720, &node719, diagramgroup720, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node721 = graphNodes[ihel * ndiagramgroups + 720];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node721, &node720, diagramgroup721, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node722 = graphNodes[ihel * ndiagramgroups + 721];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node722, &node721, diagramgroup722, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node723 = graphNodes[ihel * ndiagramgroups + 722];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node723, &node722, diagramgroup723, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node724 = graphNodes[ihel * ndiagramgroups + 723];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node724, &node723, diagramgroup724, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node725 = graphNodes[ihel * ndiagramgroups + 724];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node725, &node724, diagramgroup725, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node726 = graphNodes[ihel * ndiagramgroups + 725];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node726, &node725, diagramgroup726, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node727 = graphNodes[ihel * ndiagramgroups + 726];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node727, &node726, diagramgroup727, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node728 = graphNodes[ihel * ndiagramgroups + 727];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node728, &node727, diagramgroup728, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node729 = graphNodes[ihel * ndiagramgroups + 728];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node729, &node728, diagramgroup729, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node730 = graphNodes[ihel * ndiagramgroups + 729];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node730, &node729, diagramgroup730, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node731 = graphNodes[ihel * ndiagramgroups + 730];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node731, &node730, diagramgroup731, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node732 = graphNodes[ihel * ndiagramgroups + 731];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node732, &node731, diagramgroup732, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node733 = graphNodes[ihel * ndiagramgroups + 732];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node733, &node732, diagramgroup733, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node734 = graphNodes[ihel * ndiagramgroups + 733];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node734, &node733, diagramgroup734, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node735 = graphNodes[ihel * ndiagramgroups + 734];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node735, &node734, diagramgroup735, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node736 = graphNodes[ihel * ndiagramgroups + 735];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node736, &node735, diagramgroup736, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node737 = graphNodes[ihel * ndiagramgroups + 736];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node737, &node736, diagramgroup737, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node738 = graphNodes[ihel * ndiagramgroups + 737];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node738, &node737, diagramgroup738, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node739 = graphNodes[ihel * ndiagramgroups + 738];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node739, &node738, diagramgroup739, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node740 = graphNodes[ihel * ndiagramgroups + 739];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node740, &node739, diagramgroup740, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node741 = graphNodes[ihel * ndiagramgroups + 740];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node741, &node740, diagramgroup741, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node742 = graphNodes[ihel * ndiagramgroups + 741];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node742, &node741, diagramgroup742, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node743 = graphNodes[ihel * ndiagramgroups + 742];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node743, &node742, diagramgroup743, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node744 = graphNodes[ihel * ndiagramgroups + 743];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node744, &node743, diagramgroup744, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node745 = graphNodes[ihel * ndiagramgroups + 744];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node745, &node744, diagramgroup745, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node746 = graphNodes[ihel * ndiagramgroups + 745];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node746, &node745, diagramgroup746, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node747 = graphNodes[ihel * ndiagramgroups + 746];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node747, &node746, diagramgroup747, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node748 = graphNodes[ihel * ndiagramgroups + 747];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node748, &node747, diagramgroup748, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node749 = graphNodes[ihel * ndiagramgroups + 748];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node749, &node748, diagramgroup749, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node750 = graphNodes[ihel * ndiagramgroups + 749];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node750, &node749, diagramgroup750, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node751 = graphNodes[ihel * ndiagramgroups + 750];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node751, &node750, diagramgroup751, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node752 = graphNodes[ihel * ndiagramgroups + 751];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node752, &node751, diagramgroup752, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node753 = graphNodes[ihel * ndiagramgroups + 752];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node753, &node752, diagramgroup753, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node754 = graphNodes[ihel * ndiagramgroups + 753];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node754, &node753, diagramgroup754, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node755 = graphNodes[ihel * ndiagramgroups + 754];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node755, &node754, diagramgroup755, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node756 = graphNodes[ihel * ndiagramgroups + 755];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node756, &node755, diagramgroup756, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node757 = graphNodes[ihel * ndiagramgroups + 756];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node757, &node756, diagramgroup757, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node758 = graphNodes[ihel * ndiagramgroups + 757];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node758, &node757, diagramgroup758, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node759 = graphNodes[ihel * ndiagramgroups + 758];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node759, &node758, diagramgroup759, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node760 = graphNodes[ihel * ndiagramgroups + 759];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node760, &node759, diagramgroup760, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node761 = graphNodes[ihel * ndiagramgroups + 760];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node761, &node760, diagramgroup761, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node762 = graphNodes[ihel * ndiagramgroups + 761];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node762, &node761, diagramgroup762, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node763 = graphNodes[ihel * ndiagramgroups + 762];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node763, &node762, diagramgroup763, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node764 = graphNodes[ihel * ndiagramgroups + 763];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node764, &node763, diagramgroup764, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node765 = graphNodes[ihel * ndiagramgroups + 764];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node765, &node764, diagramgroup765, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node766 = graphNodes[ihel * ndiagramgroups + 765];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node766, &node765, diagramgroup766, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node767 = graphNodes[ihel * ndiagramgroups + 766];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node767, &node766, diagramgroup767, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node768 = graphNodes[ihel * ndiagramgroups + 767];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node768, &node767, diagramgroup768, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node769 = graphNodes[ihel * ndiagramgroups + 768];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node769, &node768, diagramgroup769, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node770 = graphNodes[ihel * ndiagramgroups + 769];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node770, &node769, diagramgroup770, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node771 = graphNodes[ihel * ndiagramgroups + 770];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node771, &node770, diagramgroup771, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node772 = graphNodes[ihel * ndiagramgroups + 771];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node772, &node771, diagramgroup772, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node773 = graphNodes[ihel * ndiagramgroups + 772];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node773, &node772, diagramgroup773, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node774 = graphNodes[ihel * ndiagramgroups + 773];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node774, &node773, diagramgroup774, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node775 = graphNodes[ihel * ndiagramgroups + 774];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node775, &node774, diagramgroup775, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node776 = graphNodes[ihel * ndiagramgroups + 775];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node776, &node775, diagramgroup776, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node777 = graphNodes[ihel * ndiagramgroups + 776];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node777, &node776, diagramgroup777, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node778 = graphNodes[ihel * ndiagramgroups + 777];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node778, &node777, diagramgroup778, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node779 = graphNodes[ihel * ndiagramgroups + 778];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node779, &node778, diagramgroup779, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node780 = graphNodes[ihel * ndiagramgroups + 779];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node780, &node779, diagramgroup780, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node781 = graphNodes[ihel * ndiagramgroups + 780];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node781, &node780, diagramgroup781, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node782 = graphNodes[ihel * ndiagramgroups + 781];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node782, &node781, diagramgroup782, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node783 = graphNodes[ihel * ndiagramgroups + 782];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node783, &node782, diagramgroup783, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node784 = graphNodes[ihel * ndiagramgroups + 783];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node784, &node783, diagramgroup784, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node785 = graphNodes[ihel * ndiagramgroups + 784];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node785, &node784, diagramgroup785, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node786 = graphNodes[ihel * ndiagramgroups + 785];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node786, &node785, diagramgroup786, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node787 = graphNodes[ihel * ndiagramgroups + 786];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node787, &node786, diagramgroup787, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node788 = graphNodes[ihel * ndiagramgroups + 787];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node788, &node787, diagramgroup788, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node789 = graphNodes[ihel * ndiagramgroups + 788];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node789, &node788, diagramgroup789, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node790 = graphNodes[ihel * ndiagramgroups + 789];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node790, &node789, diagramgroup790, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node791 = graphNodes[ihel * ndiagramgroups + 790];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node791, &node790, diagramgroup791, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node792 = graphNodes[ihel * ndiagramgroups + 791];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node792, &node791, diagramgroup792, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node793 = graphNodes[ihel * ndiagramgroups + 792];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node793, &node792, diagramgroup793, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node794 = graphNodes[ihel * ndiagramgroups + 793];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node794, &node793, diagramgroup794, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node795 = graphNodes[ihel * ndiagramgroups + 794];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node795, &node794, diagramgroup795, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node796 = graphNodes[ihel * ndiagramgroups + 795];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node796, &node795, diagramgroup796, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node797 = graphNodes[ihel * ndiagramgroups + 796];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node797, &node796, diagramgroup797, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node798 = graphNodes[ihel * ndiagramgroups + 797];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node798, &node797, diagramgroup798, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node799 = graphNodes[ihel * ndiagramgroups + 798];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node799, &node798, diagramgroup799, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node800 = graphNodes[ihel * ndiagramgroups + 799];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node800, &node799, diagramgroup800, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node801 = graphNodes[ihel * ndiagramgroups + 800];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node801, &node800, diagramgroup801, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node802 = graphNodes[ihel * ndiagramgroups + 801];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node802, &node801, diagramgroup802, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node803 = graphNodes[ihel * ndiagramgroups + 802];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node803, &node802, diagramgroup803, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node804 = graphNodes[ihel * ndiagramgroups + 803];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node804, &node803, diagramgroup804, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node805 = graphNodes[ihel * ndiagramgroups + 804];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node805, &node804, diagramgroup805, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node806 = graphNodes[ihel * ndiagramgroups + 805];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node806, &node805, diagramgroup806, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node807 = graphNodes[ihel * ndiagramgroups + 806];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node807, &node806, diagramgroup807, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node808 = graphNodes[ihel * ndiagramgroups + 807];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node808, &node807, diagramgroup808, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node809 = graphNodes[ihel * ndiagramgroups + 808];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node809, &node808, diagramgroup809, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node810 = graphNodes[ihel * ndiagramgroups + 809];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node810, &node809, diagramgroup810, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node811 = graphNodes[ihel * ndiagramgroups + 810];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node811, &node810, diagramgroup811, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node812 = graphNodes[ihel * ndiagramgroups + 811];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node812, &node811, diagramgroup812, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node813 = graphNodes[ihel * ndiagramgroups + 812];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node813, &node812, diagramgroup813, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node814 = graphNodes[ihel * ndiagramgroups + 813];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node814, &node813, diagramgroup814, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node815 = graphNodes[ihel * ndiagramgroups + 814];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node815, &node814, diagramgroup815, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node816 = graphNodes[ihel * ndiagramgroups + 815];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node816, &node815, diagramgroup816, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node817 = graphNodes[ihel * ndiagramgroups + 816];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node817, &node816, diagramgroup817, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node818 = graphNodes[ihel * ndiagramgroups + 817];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node818, &node817, diagramgroup818, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node819 = graphNodes[ihel * ndiagramgroups + 818];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node819, &node818, diagramgroup819, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node820 = graphNodes[ihel * ndiagramgroups + 819];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node820, &node819, diagramgroup820, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node821 = graphNodes[ihel * ndiagramgroups + 820];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node821, &node820, diagramgroup821, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node822 = graphNodes[ihel * ndiagramgroups + 821];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node822, &node821, diagramgroup822, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node823 = graphNodes[ihel * ndiagramgroups + 822];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node823, &node822, diagramgroup823, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node824 = graphNodes[ihel * ndiagramgroups + 823];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node824, &node823, diagramgroup824, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node825 = graphNodes[ihel * ndiagramgroups + 824];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node825, &node824, diagramgroup825, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node826 = graphNodes[ihel * ndiagramgroups + 825];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node826, &node825, diagramgroup826, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node827 = graphNodes[ihel * ndiagramgroups + 826];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node827, &node826, diagramgroup827, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node828 = graphNodes[ihel * ndiagramgroups + 827];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node828, &node827, diagramgroup828, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node829 = graphNodes[ihel * ndiagramgroups + 828];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node829, &node828, diagramgroup829, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node830 = graphNodes[ihel * ndiagramgroups + 829];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node830, &node829, diagramgroup830, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node831 = graphNodes[ihel * ndiagramgroups + 830];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node831, &node830, diagramgroup831, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node832 = graphNodes[ihel * ndiagramgroups + 831];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node832, &node831, diagramgroup832, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node833 = graphNodes[ihel * ndiagramgroups + 832];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node833, &node832, diagramgroup833, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node834 = graphNodes[ihel * ndiagramgroups + 833];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node834, &node833, diagramgroup834, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node835 = graphNodes[ihel * ndiagramgroups + 834];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node835, &node834, diagramgroup835, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node836 = graphNodes[ihel * ndiagramgroups + 835];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node836, &node835, diagramgroup836, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node837 = graphNodes[ihel * ndiagramgroups + 836];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node837, &node836, diagramgroup837, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node838 = graphNodes[ihel * ndiagramgroups + 837];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node838, &node837, diagramgroup838, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node839 = graphNodes[ihel * ndiagramgroups + 838];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node839, &node838, diagramgroup839, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node840 = graphNodes[ihel * ndiagramgroups + 839];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node840, &node839, diagramgroup840, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node841 = graphNodes[ihel * ndiagramgroups + 840];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node841, &node840, diagramgroup841, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node842 = graphNodes[ihel * ndiagramgroups + 841];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node842, &node841, diagramgroup842, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node843 = graphNodes[ihel * ndiagramgroups + 842];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node843, &node842, diagramgroup843, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node844 = graphNodes[ihel * ndiagramgroups + 843];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node844, &node843, diagramgroup844, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node845 = graphNodes[ihel * ndiagramgroups + 844];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node845, &node844, diagramgroup845, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node846 = graphNodes[ihel * ndiagramgroups + 845];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node846, &node845, diagramgroup846, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node847 = graphNodes[ihel * ndiagramgroups + 846];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node847, &node846, diagramgroup847, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node848 = graphNodes[ihel * ndiagramgroups + 847];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node848, &node847, diagramgroup848, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node849 = graphNodes[ihel * ndiagramgroups + 848];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node849, &node848, diagramgroup849, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node850 = graphNodes[ihel * ndiagramgroups + 849];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node850, &node849, diagramgroup850, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node851 = graphNodes[ihel * ndiagramgroups + 850];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node851, &node850, diagramgroup851, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node852 = graphNodes[ihel * ndiagramgroups + 851];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node852, &node851, diagramgroup852, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node853 = graphNodes[ihel * ndiagramgroups + 852];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node853, &node852, diagramgroup853, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node854 = graphNodes[ihel * ndiagramgroups + 853];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node854, &node853, diagramgroup854, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node855 = graphNodes[ihel * ndiagramgroups + 854];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node855, &node854, diagramgroup855, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node856 = graphNodes[ihel * ndiagramgroups + 855];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node856, &node855, diagramgroup856, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node857 = graphNodes[ihel * ndiagramgroups + 856];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node857, &node856, diagramgroup857, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node858 = graphNodes[ihel * ndiagramgroups + 857];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node858, &node857, diagramgroup858, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node859 = graphNodes[ihel * ndiagramgroups + 858];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node859, &node858, diagramgroup859, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node860 = graphNodes[ihel * ndiagramgroups + 859];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node860, &node859, diagramgroup860, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node861 = graphNodes[ihel * ndiagramgroups + 860];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node861, &node860, diagramgroup861, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node862 = graphNodes[ihel * ndiagramgroups + 861];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node862, &node861, diagramgroup862, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node863 = graphNodes[ihel * ndiagramgroups + 862];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node863, &node862, diagramgroup863, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node864 = graphNodes[ihel * ndiagramgroups + 863];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node864, &node863, diagramgroup864, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node865 = graphNodes[ihel * ndiagramgroups + 864];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node865, &node864, diagramgroup865, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node866 = graphNodes[ihel * ndiagramgroups + 865];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node866, &node865, diagramgroup866, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node867 = graphNodes[ihel * ndiagramgroups + 866];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node867, &node866, diagramgroup867, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node868 = graphNodes[ihel * ndiagramgroups + 867];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node868, &node867, diagramgroup868, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node869 = graphNodes[ihel * ndiagramgroups + 868];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node869, &node868, diagramgroup869, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node870 = graphNodes[ihel * ndiagramgroups + 869];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node870, &node869, diagramgroup870, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node871 = graphNodes[ihel * ndiagramgroups + 870];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node871, &node870, diagramgroup871, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node872 = graphNodes[ihel * ndiagramgroups + 871];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node872, &node871, diagramgroup872, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node873 = graphNodes[ihel * ndiagramgroups + 872];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node873, &node872, diagramgroup873, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node874 = graphNodes[ihel * ndiagramgroups + 873];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node874, &node873, diagramgroup874, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node875 = graphNodes[ihel * ndiagramgroups + 874];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node875, &node874, diagramgroup875, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node876 = graphNodes[ihel * ndiagramgroups + 875];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node876, &node875, diagramgroup876, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node877 = graphNodes[ihel * ndiagramgroups + 876];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node877, &node876, diagramgroup877, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node878 = graphNodes[ihel * ndiagramgroups + 877];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node878, &node877, diagramgroup878, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node879 = graphNodes[ihel * ndiagramgroups + 878];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node879, &node878, diagramgroup879, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node880 = graphNodes[ihel * ndiagramgroups + 879];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node880, &node879, diagramgroup880, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node881 = graphNodes[ihel * ndiagramgroups + 880];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node881, &node880, diagramgroup881, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node882 = graphNodes[ihel * ndiagramgroups + 881];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node882, &node881, diagramgroup882, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node883 = graphNodes[ihel * ndiagramgroups + 882];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node883, &node882, diagramgroup883, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node884 = graphNodes[ihel * ndiagramgroups + 883];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node884, &node883, diagramgroup884, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node885 = graphNodes[ihel * ndiagramgroups + 884];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node885, &node884, diagramgroup885, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node886 = graphNodes[ihel * ndiagramgroups + 885];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node886, &node885, diagramgroup886, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node887 = graphNodes[ihel * ndiagramgroups + 886];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node887, &node886, diagramgroup887, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node888 = graphNodes[ihel * ndiagramgroups + 887];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node888, &node887, diagramgroup888, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node889 = graphNodes[ihel * ndiagramgroups + 888];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node889, &node888, diagramgroup889, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node890 = graphNodes[ihel * ndiagramgroups + 889];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node890, &node889, diagramgroup890, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node891 = graphNodes[ihel * ndiagramgroups + 890];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node891, &node890, diagramgroup891, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node892 = graphNodes[ihel * ndiagramgroups + 891];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node892, &node891, diagramgroup892, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node893 = graphNodes[ihel * ndiagramgroups + 892];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node893, &node892, diagramgroup893, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node894 = graphNodes[ihel * ndiagramgroups + 893];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node894, &node893, diagramgroup894, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node895 = graphNodes[ihel * ndiagramgroups + 894];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node895, &node894, diagramgroup895, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node896 = graphNodes[ihel * ndiagramgroups + 895];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node896, &node895, diagramgroup896, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node897 = graphNodes[ihel * ndiagramgroups + 896];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node897, &node896, diagramgroup897, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node898 = graphNodes[ihel * ndiagramgroups + 897];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node898, &node897, diagramgroup898, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node899 = graphNodes[ihel * ndiagramgroups + 898];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node899, &node898, diagramgroup899, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node900 = graphNodes[ihel * ndiagramgroups + 899];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node900, &node899, diagramgroup900, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node901 = graphNodes[ihel * ndiagramgroups + 900];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node901, &node900, diagramgroup901, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node902 = graphNodes[ihel * ndiagramgroups + 901];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node902, &node901, diagramgroup902, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node903 = graphNodes[ihel * ndiagramgroups + 902];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node903, &node902, diagramgroup903, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node904 = graphNodes[ihel * ndiagramgroups + 903];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node904, &node903, diagramgroup904, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node905 = graphNodes[ihel * ndiagramgroups + 904];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node905, &node904, diagramgroup905, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node906 = graphNodes[ihel * ndiagramgroups + 905];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node906, &node905, diagramgroup906, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node907 = graphNodes[ihel * ndiagramgroups + 906];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node907, &node906, diagramgroup907, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node908 = graphNodes[ihel * ndiagramgroups + 907];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node908, &node907, diagramgroup908, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node909 = graphNodes[ihel * ndiagramgroups + 908];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node909, &node908, diagramgroup909, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node910 = graphNodes[ihel * ndiagramgroups + 909];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node910, &node909, diagramgroup910, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node911 = graphNodes[ihel * ndiagramgroups + 910];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node911, &node910, diagramgroup911, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node912 = graphNodes[ihel * ndiagramgroups + 911];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node912, &node911, diagramgroup912, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node913 = graphNodes[ihel * ndiagramgroups + 912];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node913, &node912, diagramgroup913, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node914 = graphNodes[ihel * ndiagramgroups + 913];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node914, &node913, diagramgroup914, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node915 = graphNodes[ihel * ndiagramgroups + 914];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node915, &node914, diagramgroup915, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node916 = graphNodes[ihel * ndiagramgroups + 915];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node916, &node915, diagramgroup916, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node917 = graphNodes[ihel * ndiagramgroups + 916];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node917, &node916, diagramgroup917, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node918 = graphNodes[ihel * ndiagramgroups + 917];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node918, &node917, diagramgroup918, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node919 = graphNodes[ihel * ndiagramgroups + 918];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node919, &node918, diagramgroup919, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node920 = graphNodes[ihel * ndiagramgroups + 919];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node920, &node919, diagramgroup920, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node921 = graphNodes[ihel * ndiagramgroups + 920];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node921, &node920, diagramgroup921, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node922 = graphNodes[ihel * ndiagramgroups + 921];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node922, &node921, diagramgroup922, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node923 = graphNodes[ihel * ndiagramgroups + 922];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node923, &node922, diagramgroup923, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node924 = graphNodes[ihel * ndiagramgroups + 923];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node924, &node923, diagramgroup924, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node925 = graphNodes[ihel * ndiagramgroups + 924];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node925, &node924, diagramgroup925, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node926 = graphNodes[ihel * ndiagramgroups + 925];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node926, &node925, diagramgroup926, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node927 = graphNodes[ihel * ndiagramgroups + 926];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node927, &node926, diagramgroup927, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node928 = graphNodes[ihel * ndiagramgroups + 927];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node928, &node927, diagramgroup928, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node929 = graphNodes[ihel * ndiagramgroups + 928];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node929, &node928, diagramgroup929, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node930 = graphNodes[ihel * ndiagramgroups + 929];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node930, &node929, diagramgroup930, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node931 = graphNodes[ihel * ndiagramgroups + 930];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node931, &node930, diagramgroup931, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node932 = graphNodes[ihel * ndiagramgroups + 931];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node932, &node931, diagramgroup932, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node933 = graphNodes[ihel * ndiagramgroups + 932];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node933, &node932, diagramgroup933, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node934 = graphNodes[ihel * ndiagramgroups + 933];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node934, &node933, diagramgroup934, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node935 = graphNodes[ihel * ndiagramgroups + 934];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node935, &node934, diagramgroup935, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node936 = graphNodes[ihel * ndiagramgroups + 935];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node936, &node935, diagramgroup936, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node937 = graphNodes[ihel * ndiagramgroups + 936];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node937, &node936, diagramgroup937, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node938 = graphNodes[ihel * ndiagramgroups + 937];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node938, &node937, diagramgroup938, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node939 = graphNodes[ihel * ndiagramgroups + 938];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node939, &node938, diagramgroup939, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node940 = graphNodes[ihel * ndiagramgroups + 939];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node940, &node939, diagramgroup940, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node941 = graphNodes[ihel * ndiagramgroups + 940];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node941, &node940, diagramgroup941, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node942 = graphNodes[ihel * ndiagramgroups + 941];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node942, &node941, diagramgroup942, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node943 = graphNodes[ihel * ndiagramgroups + 942];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node943, &node942, diagramgroup943, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node944 = graphNodes[ihel * ndiagramgroups + 943];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node944, &node943, diagramgroup944, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node945 = graphNodes[ihel * ndiagramgroups + 944];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node945, &node944, diagramgroup945, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node946 = graphNodes[ihel * ndiagramgroups + 945];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node946, &node945, diagramgroup946, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node947 = graphNodes[ihel * ndiagramgroups + 946];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node947, &node946, diagramgroup947, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node948 = graphNodes[ihel * ndiagramgroups + 947];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node948, &node947, diagramgroup948, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node949 = graphNodes[ihel * ndiagramgroups + 948];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node949, &node948, diagramgroup949, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node950 = graphNodes[ihel * ndiagramgroups + 949];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node950, &node949, diagramgroup950, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node951 = graphNodes[ihel * ndiagramgroups + 950];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node951, &node950, diagramgroup951, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node952 = graphNodes[ihel * ndiagramgroups + 951];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node952, &node951, diagramgroup952, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node953 = graphNodes[ihel * ndiagramgroups + 952];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node953, &node952, diagramgroup953, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node954 = graphNodes[ihel * ndiagramgroups + 953];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node954, &node953, diagramgroup954, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node955 = graphNodes[ihel * ndiagramgroups + 954];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node955, &node954, diagramgroup955, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node956 = graphNodes[ihel * ndiagramgroups + 955];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node956, &node955, diagramgroup956, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node957 = graphNodes[ihel * ndiagramgroups + 956];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node957, &node956, diagramgroup957, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node958 = graphNodes[ihel * ndiagramgroups + 957];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node958, &node957, diagramgroup958, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node959 = graphNodes[ihel * ndiagramgroups + 958];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node959, &node958, diagramgroup959, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node960 = graphNodes[ihel * ndiagramgroups + 959];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node960, &node959, diagramgroup960, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node961 = graphNodes[ihel * ndiagramgroups + 960];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node961, &node960, diagramgroup961, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node962 = graphNodes[ihel * ndiagramgroups + 961];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node962, &node961, diagramgroup962, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node963 = graphNodes[ihel * ndiagramgroups + 962];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node963, &node962, diagramgroup963, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node964 = graphNodes[ihel * ndiagramgroups + 963];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node964, &node963, diagramgroup964, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node965 = graphNodes[ihel * ndiagramgroups + 964];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node965, &node964, diagramgroup965, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node966 = graphNodes[ihel * ndiagramgroups + 965];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node966, &node965, diagramgroup966, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node967 = graphNodes[ihel * ndiagramgroups + 966];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node967, &node966, diagramgroup967, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node968 = graphNodes[ihel * ndiagramgroups + 967];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node968, &node967, diagramgroup968, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node969 = graphNodes[ihel * ndiagramgroups + 968];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node969, &node968, diagramgroup969, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node970 = graphNodes[ihel * ndiagramgroups + 969];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node970, &node969, diagramgroup970, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node971 = graphNodes[ihel * ndiagramgroups + 970];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node971, &node970, diagramgroup971, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node972 = graphNodes[ihel * ndiagramgroups + 971];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node972, &node971, diagramgroup972, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node973 = graphNodes[ihel * ndiagramgroups + 972];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node973, &node972, diagramgroup973, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node974 = graphNodes[ihel * ndiagramgroups + 973];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node974, &node973, diagramgroup974, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node975 = graphNodes[ihel * ndiagramgroups + 974];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node975, &node974, diagramgroup975, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node976 = graphNodes[ihel * ndiagramgroups + 975];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node976, &node975, diagramgroup976, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node977 = graphNodes[ihel * ndiagramgroups + 976];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node977, &node976, diagramgroup977, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node978 = graphNodes[ihel * ndiagramgroups + 977];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node978, &node977, diagramgroup978, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node979 = graphNodes[ihel * ndiagramgroups + 978];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node979, &node978, diagramgroup979, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node980 = graphNodes[ihel * ndiagramgroups + 979];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node980, &node979, diagramgroup980, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node981 = graphNodes[ihel * ndiagramgroups + 980];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node981, &node980, diagramgroup981, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node982 = graphNodes[ihel * ndiagramgroups + 981];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node982, &node981, diagramgroup982, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node983 = graphNodes[ihel * ndiagramgroups + 982];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node983, &node982, diagramgroup983, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node984 = graphNodes[ihel * ndiagramgroups + 983];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node984, &node983, diagramgroup984, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node985 = graphNodes[ihel * ndiagramgroups + 984];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node985, &node984, diagramgroup985, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node986 = graphNodes[ihel * ndiagramgroups + 985];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node986, &node985, diagramgroup986, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node987 = graphNodes[ihel * ndiagramgroups + 986];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node987, &node986, diagramgroup987, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node988 = graphNodes[ihel * ndiagramgroups + 987];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node988, &node987, diagramgroup988, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node989 = graphNodes[ihel * ndiagramgroups + 988];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node989, &node988, diagramgroup989, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node990 = graphNodes[ihel * ndiagramgroups + 989];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node990, &node989, diagramgroup990, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node991 = graphNodes[ihel * ndiagramgroups + 990];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node991, &node990, diagramgroup991, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node992 = graphNodes[ihel * ndiagramgroups + 991];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node992, &node991, diagramgroup992, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node993 = graphNodes[ihel * ndiagramgroups + 992];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node993, &node992, diagramgroup993, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node994 = graphNodes[ihel * ndiagramgroups + 993];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node994, &node993, diagramgroup994, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node995 = graphNodes[ihel * ndiagramgroups + 994];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node995, &node994, diagramgroup995, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node996 = graphNodes[ihel * ndiagramgroups + 995];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node996, &node995, diagramgroup996, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node997 = graphNodes[ihel * ndiagramgroups + 996];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node997, &node996, diagramgroup997, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node998 = graphNodes[ihel * ndiagramgroups + 997];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node998, &node997, diagramgroup998, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node999 = graphNodes[ihel * ndiagramgroups + 998];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node999, &node998, diagramgroup999, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1000 = graphNodes[ihel * ndiagramgroups + 999];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1000, &node999, diagramgroup1000, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1001 = graphNodes[ihel * ndiagramgroups + 1000];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1001, &node1000, diagramgroup1001, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1002 = graphNodes[ihel * ndiagramgroups + 1001];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1002, &node1001, diagramgroup1002, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1003 = graphNodes[ihel * ndiagramgroups + 1002];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1003, &node1002, diagramgroup1003, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1004 = graphNodes[ihel * ndiagramgroups + 1003];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1004, &node1003, diagramgroup1004, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1005 = graphNodes[ihel * ndiagramgroups + 1004];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1005, &node1004, diagramgroup1005, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1006 = graphNodes[ihel * ndiagramgroups + 1005];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1006, &node1005, diagramgroup1006, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1007 = graphNodes[ihel * ndiagramgroups + 1006];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1007, &node1006, diagramgroup1007, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1008 = graphNodes[ihel * ndiagramgroups + 1007];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1008, &node1007, diagramgroup1008, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1009 = graphNodes[ihel * ndiagramgroups + 1008];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1009, &node1008, diagramgroup1009, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1010 = graphNodes[ihel * ndiagramgroups + 1009];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1010, &node1009, diagramgroup1010, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1011 = graphNodes[ihel * ndiagramgroups + 1010];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1011, &node1010, diagramgroup1011, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1012 = graphNodes[ihel * ndiagramgroups + 1011];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1012, &node1011, diagramgroup1012, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1013 = graphNodes[ihel * ndiagramgroups + 1012];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1013, &node1012, diagramgroup1013, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1014 = graphNodes[ihel * ndiagramgroups + 1013];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1014, &node1013, diagramgroup1014, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1015 = graphNodes[ihel * ndiagramgroups + 1014];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1015, &node1014, diagramgroup1015, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1016 = graphNodes[ihel * ndiagramgroups + 1015];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1016, &node1015, diagramgroup1016, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1017 = graphNodes[ihel * ndiagramgroups + 1016];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1017, &node1016, diagramgroup1017, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1018 = graphNodes[ihel * ndiagramgroups + 1017];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1018, &node1017, diagramgroup1018, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1019 = graphNodes[ihel * ndiagramgroups + 1018];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1019, &node1018, diagramgroup1019, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1020 = graphNodes[ihel * ndiagramgroups + 1019];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1020, &node1019, diagramgroup1020, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1021 = graphNodes[ihel * ndiagramgroups + 1020];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1021, &node1020, diagramgroup1021, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1022 = graphNodes[ihel * ndiagramgroups + 1021];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1022, &node1021, diagramgroup1022, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1023 = graphNodes[ihel * ndiagramgroups + 1022];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1023, &node1022, diagramgroup1023, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1024 = graphNodes[ihel * ndiagramgroups + 1023];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1024, &node1023, diagramgroup1024, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1025 = graphNodes[ihel * ndiagramgroups + 1024];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1025, &node1024, diagramgroup1025, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1026 = graphNodes[ihel * ndiagramgroups + 1025];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1026, &node1025, diagramgroup1026, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1027 = graphNodes[ihel * ndiagramgroups + 1026];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1027, &node1026, diagramgroup1027, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1028 = graphNodes[ihel * ndiagramgroups + 1027];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1028, &node1027, diagramgroup1028, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1029 = graphNodes[ihel * ndiagramgroups + 1028];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1029, &node1028, diagramgroup1029, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1030 = graphNodes[ihel * ndiagramgroups + 1029];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1030, &node1029, diagramgroup1030, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1031 = graphNodes[ihel * ndiagramgroups + 1030];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1031, &node1030, diagramgroup1031, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1032 = graphNodes[ihel * ndiagramgroups + 1031];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1032, &node1031, diagramgroup1032, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1033 = graphNodes[ihel * ndiagramgroups + 1032];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1033, &node1032, diagramgroup1033, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1034 = graphNodes[ihel * ndiagramgroups + 1033];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1034, &node1033, diagramgroup1034, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1035 = graphNodes[ihel * ndiagramgroups + 1034];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1035, &node1034, diagramgroup1035, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1036 = graphNodes[ihel * ndiagramgroups + 1035];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1036, &node1035, diagramgroup1036, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1037 = graphNodes[ihel * ndiagramgroups + 1036];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1037, &node1036, diagramgroup1037, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1038 = graphNodes[ihel * ndiagramgroups + 1037];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1038, &node1037, diagramgroup1038, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1039 = graphNodes[ihel * ndiagramgroups + 1038];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1039, &node1038, diagramgroup1039, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1040 = graphNodes[ihel * ndiagramgroups + 1039];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1040, &node1039, diagramgroup1040, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1041 = graphNodes[ihel * ndiagramgroups + 1040];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1041, &node1040, diagramgroup1041, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1042 = graphNodes[ihel * ndiagramgroups + 1041];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1042, &node1041, diagramgroup1042, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1043 = graphNodes[ihel * ndiagramgroups + 1042];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1043, &node1042, diagramgroup1043, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1044 = graphNodes[ihel * ndiagramgroups + 1043];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1044, &node1043, diagramgroup1044, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1045 = graphNodes[ihel * ndiagramgroups + 1044];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1045, &node1044, diagramgroup1045, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1046 = graphNodes[ihel * ndiagramgroups + 1045];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1046, &node1045, diagramgroup1046, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1047 = graphNodes[ihel * ndiagramgroups + 1046];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1047, &node1046, diagramgroup1047, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1048 = graphNodes[ihel * ndiagramgroups + 1047];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1048, &node1047, diagramgroup1048, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1049 = graphNodes[ihel * ndiagramgroups + 1048];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1049, &node1048, diagramgroup1049, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1050 = graphNodes[ihel * ndiagramgroups + 1049];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1050, &node1049, diagramgroup1050, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1051 = graphNodes[ihel * ndiagramgroups + 1050];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1051, &node1050, diagramgroup1051, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1052 = graphNodes[ihel * ndiagramgroups + 1051];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1052, &node1051, diagramgroup1052, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1053 = graphNodes[ihel * ndiagramgroups + 1052];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1053, &node1052, diagramgroup1053, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1054 = graphNodes[ihel * ndiagramgroups + 1053];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1054, &node1053, diagramgroup1054, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1055 = graphNodes[ihel * ndiagramgroups + 1054];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1055, &node1054, diagramgroup1055, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1056 = graphNodes[ihel * ndiagramgroups + 1055];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1056, &node1055, diagramgroup1056, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1057 = graphNodes[ihel * ndiagramgroups + 1056];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1057, &node1056, diagramgroup1057, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1058 = graphNodes[ihel * ndiagramgroups + 1057];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1058, &node1057, diagramgroup1058, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1059 = graphNodes[ihel * ndiagramgroups + 1058];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1059, &node1058, diagramgroup1059, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1060 = graphNodes[ihel * ndiagramgroups + 1059];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1060, &node1059, diagramgroup1060, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1061 = graphNodes[ihel * ndiagramgroups + 1060];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1061, &node1060, diagramgroup1061, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1062 = graphNodes[ihel * ndiagramgroups + 1061];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1062, &node1061, diagramgroup1062, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1063 = graphNodes[ihel * ndiagramgroups + 1062];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1063, &node1062, diagramgroup1063, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1064 = graphNodes[ihel * ndiagramgroups + 1063];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1064, &node1063, diagramgroup1064, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1065 = graphNodes[ihel * ndiagramgroups + 1064];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1065, &node1064, diagramgroup1065, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1066 = graphNodes[ihel * ndiagramgroups + 1065];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1066, &node1065, diagramgroup1066, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1067 = graphNodes[ihel * ndiagramgroups + 1066];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1067, &node1066, diagramgroup1067, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1068 = graphNodes[ihel * ndiagramgroups + 1067];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1068, &node1067, diagramgroup1068, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1069 = graphNodes[ihel * ndiagramgroups + 1068];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1069, &node1068, diagramgroup1069, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1070 = graphNodes[ihel * ndiagramgroups + 1069];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1070, &node1069, diagramgroup1070, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1071 = graphNodes[ihel * ndiagramgroups + 1070];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1071, &node1070, diagramgroup1071, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1072 = graphNodes[ihel * ndiagramgroups + 1071];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1072, &node1071, diagramgroup1072, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1073 = graphNodes[ihel * ndiagramgroups + 1072];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1073, &node1072, diagramgroup1073, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1074 = graphNodes[ihel * ndiagramgroups + 1073];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1074, &node1073, diagramgroup1074, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1075 = graphNodes[ihel * ndiagramgroups + 1074];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1075, &node1074, diagramgroup1075, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1076 = graphNodes[ihel * ndiagramgroups + 1075];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1076, &node1075, diagramgroup1076, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1077 = graphNodes[ihel * ndiagramgroups + 1076];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1077, &node1076, diagramgroup1077, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1078 = graphNodes[ihel * ndiagramgroups + 1077];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1078, &node1077, diagramgroup1078, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1079 = graphNodes[ihel * ndiagramgroups + 1078];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1079, &node1078, diagramgroup1079, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1080 = graphNodes[ihel * ndiagramgroups + 1079];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1080, &node1079, diagramgroup1080, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1081 = graphNodes[ihel * ndiagramgroups + 1080];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1081, &node1080, diagramgroup1081, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1082 = graphNodes[ihel * ndiagramgroups + 1081];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1082, &node1081, diagramgroup1082, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1083 = graphNodes[ihel * ndiagramgroups + 1082];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1083, &node1082, diagramgroup1083, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1084 = graphNodes[ihel * ndiagramgroups + 1083];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1084, &node1083, diagramgroup1084, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1085 = graphNodes[ihel * ndiagramgroups + 1084];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1085, &node1084, diagramgroup1085, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1086 = graphNodes[ihel * ndiagramgroups + 1085];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1086, &node1085, diagramgroup1086, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1087 = graphNodes[ihel * ndiagramgroups + 1086];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1087, &node1086, diagramgroup1087, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1088 = graphNodes[ihel * ndiagramgroups + 1087];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1088, &node1087, diagramgroup1088, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1089 = graphNodes[ihel * ndiagramgroups + 1088];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1089, &node1088, diagramgroup1089, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1090 = graphNodes[ihel * ndiagramgroups + 1089];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1090, &node1089, diagramgroup1090, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1091 = graphNodes[ihel * ndiagramgroups + 1090];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1091, &node1090, diagramgroup1091, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1092 = graphNodes[ihel * ndiagramgroups + 1091];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1092, &node1091, diagramgroup1092, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1093 = graphNodes[ihel * ndiagramgroups + 1092];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1093, &node1092, diagramgroup1093, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1094 = graphNodes[ihel * ndiagramgroups + 1093];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1094, &node1093, diagramgroup1094, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1095 = graphNodes[ihel * ndiagramgroups + 1094];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1095, &node1094, diagramgroup1095, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1096 = graphNodes[ihel * ndiagramgroups + 1095];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1096, &node1095, diagramgroup1096, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1097 = graphNodes[ihel * ndiagramgroups + 1096];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1097, &node1096, diagramgroup1097, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1098 = graphNodes[ihel * ndiagramgroups + 1097];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1098, &node1097, diagramgroup1098, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1099 = graphNodes[ihel * ndiagramgroups + 1098];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1099, &node1098, diagramgroup1099, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1100 = graphNodes[ihel * ndiagramgroups + 1099];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1100, &node1099, diagramgroup1100, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1101 = graphNodes[ihel * ndiagramgroups + 1100];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1101, &node1100, diagramgroup1101, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1102 = graphNodes[ihel * ndiagramgroups + 1101];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1102, &node1101, diagramgroup1102, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1103 = graphNodes[ihel * ndiagramgroups + 1102];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1103, &node1102, diagramgroup1103, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1104 = graphNodes[ihel * ndiagramgroups + 1103];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1104, &node1103, diagramgroup1104, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1105 = graphNodes[ihel * ndiagramgroups + 1104];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1105, &node1104, diagramgroup1105, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1106 = graphNodes[ihel * ndiagramgroups + 1105];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1106, &node1105, diagramgroup1106, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1107 = graphNodes[ihel * ndiagramgroups + 1106];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1107, &node1106, diagramgroup1107, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1108 = graphNodes[ihel * ndiagramgroups + 1107];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1108, &node1107, diagramgroup1108, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1109 = graphNodes[ihel * ndiagramgroups + 1108];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1109, &node1108, diagramgroup1109, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1110 = graphNodes[ihel * ndiagramgroups + 1109];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1110, &node1109, diagramgroup1110, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1111 = graphNodes[ihel * ndiagramgroups + 1110];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1111, &node1110, diagramgroup1111, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1112 = graphNodes[ihel * ndiagramgroups + 1111];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1112, &node1111, diagramgroup1112, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1113 = graphNodes[ihel * ndiagramgroups + 1112];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1113, &node1112, diagramgroup1113, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1114 = graphNodes[ihel * ndiagramgroups + 1113];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1114, &node1113, diagramgroup1114, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1115 = graphNodes[ihel * ndiagramgroups + 1114];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1115, &node1114, diagramgroup1115, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1116 = graphNodes[ihel * ndiagramgroups + 1115];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1116, &node1115, diagramgroup1116, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1117 = graphNodes[ihel * ndiagramgroups + 1116];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1117, &node1116, diagramgroup1117, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1118 = graphNodes[ihel * ndiagramgroups + 1117];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1118, &node1117, diagramgroup1118, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1119 = graphNodes[ihel * ndiagramgroups + 1118];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1119, &node1118, diagramgroup1119, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1120 = graphNodes[ihel * ndiagramgroups + 1119];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1120, &node1119, diagramgroup1120, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1121 = graphNodes[ihel * ndiagramgroups + 1120];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1121, &node1120, diagramgroup1121, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1122 = graphNodes[ihel * ndiagramgroups + 1121];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1122, &node1121, diagramgroup1122, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1123 = graphNodes[ihel * ndiagramgroups + 1122];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1123, &node1122, diagramgroup1123, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1124 = graphNodes[ihel * ndiagramgroups + 1123];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1124, &node1123, diagramgroup1124, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1125 = graphNodes[ihel * ndiagramgroups + 1124];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1125, &node1124, diagramgroup1125, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1126 = graphNodes[ihel * ndiagramgroups + 1125];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1126, &node1125, diagramgroup1126, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1127 = graphNodes[ihel * ndiagramgroups + 1126];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1127, &node1126, diagramgroup1127, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1128 = graphNodes[ihel * ndiagramgroups + 1127];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1128, &node1127, diagramgroup1128, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1129 = graphNodes[ihel * ndiagramgroups + 1128];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1129, &node1128, diagramgroup1129, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1130 = graphNodes[ihel * ndiagramgroups + 1129];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1130, &node1129, diagramgroup1130, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1131 = graphNodes[ihel * ndiagramgroups + 1130];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1131, &node1130, diagramgroup1131, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1132 = graphNodes[ihel * ndiagramgroups + 1131];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1132, &node1131, diagramgroup1132, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1133 = graphNodes[ihel * ndiagramgroups + 1132];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1133, &node1132, diagramgroup1133, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1134 = graphNodes[ihel * ndiagramgroups + 1133];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1134, &node1133, diagramgroup1134, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1135 = graphNodes[ihel * ndiagramgroups + 1134];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1135, &node1134, diagramgroup1135, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1136 = graphNodes[ihel * ndiagramgroups + 1135];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1136, &node1135, diagramgroup1136, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1137 = graphNodes[ihel * ndiagramgroups + 1136];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1137, &node1136, diagramgroup1137, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1138 = graphNodes[ihel * ndiagramgroups + 1137];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1138, &node1137, diagramgroup1138, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1139 = graphNodes[ihel * ndiagramgroups + 1138];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1139, &node1138, diagramgroup1139, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1140 = graphNodes[ihel * ndiagramgroups + 1139];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1140, &node1139, diagramgroup1140, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1141 = graphNodes[ihel * ndiagramgroups + 1140];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1141, &node1140, diagramgroup1141, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1142 = graphNodes[ihel * ndiagramgroups + 1141];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1142, &node1141, diagramgroup1142, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1143 = graphNodes[ihel * ndiagramgroups + 1142];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1143, &node1142, diagramgroup1143, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1144 = graphNodes[ihel * ndiagramgroups + 1143];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1144, &node1143, diagramgroup1144, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1145 = graphNodes[ihel * ndiagramgroups + 1144];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1145, &node1144, diagramgroup1145, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1146 = graphNodes[ihel * ndiagramgroups + 1145];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1146, &node1145, diagramgroup1146, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1147 = graphNodes[ihel * ndiagramgroups + 1146];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1147, &node1146, diagramgroup1147, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1148 = graphNodes[ihel * ndiagramgroups + 1147];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1148, &node1147, diagramgroup1148, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1149 = graphNodes[ihel * ndiagramgroups + 1148];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1149, &node1148, diagramgroup1149, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1150 = graphNodes[ihel * ndiagramgroups + 1149];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1150, &node1149, diagramgroup1150, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1151 = graphNodes[ihel * ndiagramgroups + 1150];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1151, &node1150, diagramgroup1151, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1152 = graphNodes[ihel * ndiagramgroups + 1151];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1152, &node1151, diagramgroup1152, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1153 = graphNodes[ihel * ndiagramgroups + 1152];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1153, &node1152, diagramgroup1153, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1154 = graphNodes[ihel * ndiagramgroups + 1153];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1154, &node1153, diagramgroup1154, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1155 = graphNodes[ihel * ndiagramgroups + 1154];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1155, &node1154, diagramgroup1155, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1156 = graphNodes[ihel * ndiagramgroups + 1155];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1156, &node1155, diagramgroup1156, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1157 = graphNodes[ihel * ndiagramgroups + 1156];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1157, &node1156, diagramgroup1157, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1158 = graphNodes[ihel * ndiagramgroups + 1157];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1158, &node1157, diagramgroup1158, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1159 = graphNodes[ihel * ndiagramgroups + 1158];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1159, &node1158, diagramgroup1159, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1160 = graphNodes[ihel * ndiagramgroups + 1159];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1160, &node1159, diagramgroup1160, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1161 = graphNodes[ihel * ndiagramgroups + 1160];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1161, &node1160, diagramgroup1161, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1162 = graphNodes[ihel * ndiagramgroups + 1161];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1162, &node1161, diagramgroup1162, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1163 = graphNodes[ihel * ndiagramgroups + 1162];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1163, &node1162, diagramgroup1163, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1164 = graphNodes[ihel * ndiagramgroups + 1163];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1164, &node1163, diagramgroup1164, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1165 = graphNodes[ihel * ndiagramgroups + 1164];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1165, &node1164, diagramgroup1165, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1166 = graphNodes[ihel * ndiagramgroups + 1165];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1166, &node1165, diagramgroup1166, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1167 = graphNodes[ihel * ndiagramgroups + 1166];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1167, &node1166, diagramgroup1167, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1168 = graphNodes[ihel * ndiagramgroups + 1167];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1168, &node1167, diagramgroup1168, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1169 = graphNodes[ihel * ndiagramgroups + 1168];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1169, &node1168, diagramgroup1169, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1170 = graphNodes[ihel * ndiagramgroups + 1169];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1170, &node1169, diagramgroup1170, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1171 = graphNodes[ihel * ndiagramgroups + 1170];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1171, &node1170, diagramgroup1171, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1172 = graphNodes[ihel * ndiagramgroups + 1171];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1172, &node1171, diagramgroup1172, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1173 = graphNodes[ihel * ndiagramgroups + 1172];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1173, &node1172, diagramgroup1173, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1174 = graphNodes[ihel * ndiagramgroups + 1173];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1174, &node1173, diagramgroup1174, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1175 = graphNodes[ihel * ndiagramgroups + 1174];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1175, &node1174, diagramgroup1175, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1176 = graphNodes[ihel * ndiagramgroups + 1175];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1176, &node1175, diagramgroup1176, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1177 = graphNodes[ihel * ndiagramgroups + 1176];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1177, &node1176, diagramgroup1177, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1178 = graphNodes[ihel * ndiagramgroups + 1177];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1178, &node1177, diagramgroup1178, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1179 = graphNodes[ihel * ndiagramgroups + 1178];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1179, &node1178, diagramgroup1179, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1180 = graphNodes[ihel * ndiagramgroups + 1179];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1180, &node1179, diagramgroup1180, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1181 = graphNodes[ihel * ndiagramgroups + 1180];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1181, &node1180, diagramgroup1181, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1182 = graphNodes[ihel * ndiagramgroups + 1181];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1182, &node1181, diagramgroup1182, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1183 = graphNodes[ihel * ndiagramgroups + 1182];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1183, &node1182, diagramgroup1183, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1184 = graphNodes[ihel * ndiagramgroups + 1183];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1184, &node1183, diagramgroup1184, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1185 = graphNodes[ihel * ndiagramgroups + 1184];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1185, &node1184, diagramgroup1185, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1186 = graphNodes[ihel * ndiagramgroups + 1185];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1186, &node1185, diagramgroup1186, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1187 = graphNodes[ihel * ndiagramgroups + 1186];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1187, &node1186, diagramgroup1187, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1188 = graphNodes[ihel * ndiagramgroups + 1187];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1188, &node1187, diagramgroup1188, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1189 = graphNodes[ihel * ndiagramgroups + 1188];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1189, &node1188, diagramgroup1189, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1190 = graphNodes[ihel * ndiagramgroups + 1189];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1190, &node1189, diagramgroup1190, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1191 = graphNodes[ihel * ndiagramgroups + 1190];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1191, &node1190, diagramgroup1191, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1192 = graphNodes[ihel * ndiagramgroups + 1191];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1192, &node1191, diagramgroup1192, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1193 = graphNodes[ihel * ndiagramgroups + 1192];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1193, &node1192, diagramgroup1193, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1194 = graphNodes[ihel * ndiagramgroups + 1193];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1194, &node1193, diagramgroup1194, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1195 = graphNodes[ihel * ndiagramgroups + 1194];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1195, &node1194, diagramgroup1195, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1196 = graphNodes[ihel * ndiagramgroups + 1195];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1196, &node1195, diagramgroup1196, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1197 = graphNodes[ihel * ndiagramgroups + 1196];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1197, &node1196, diagramgroup1197, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1198 = graphNodes[ihel * ndiagramgroups + 1197];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1198, &node1197, diagramgroup1198, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1199 = graphNodes[ihel * ndiagramgroups + 1198];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1199, &node1198, diagramgroup1199, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1200 = graphNodes[ihel * ndiagramgroups + 1199];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1200, &node1199, diagramgroup1200, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1201 = graphNodes[ihel * ndiagramgroups + 1200];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1201, &node1200, diagramgroup1201, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1202 = graphNodes[ihel * ndiagramgroups + 1201];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1202, &node1201, diagramgroup1202, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1203 = graphNodes[ihel * ndiagramgroups + 1202];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1203, &node1202, diagramgroup1203, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1204 = graphNodes[ihel * ndiagramgroups + 1203];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1204, &node1203, diagramgroup1204, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1205 = graphNodes[ihel * ndiagramgroups + 1204];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1205, &node1204, diagramgroup1205, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1206 = graphNodes[ihel * ndiagramgroups + 1205];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1206, &node1205, diagramgroup1206, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1207 = graphNodes[ihel * ndiagramgroups + 1206];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1207, &node1206, diagramgroup1207, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1208 = graphNodes[ihel * ndiagramgroups + 1207];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1208, &node1207, diagramgroup1208, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1209 = graphNodes[ihel * ndiagramgroups + 1208];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1209, &node1208, diagramgroup1209, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1210 = graphNodes[ihel * ndiagramgroups + 1209];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1210, &node1209, diagramgroup1210, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1211 = graphNodes[ihel * ndiagramgroups + 1210];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1211, &node1210, diagramgroup1211, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1212 = graphNodes[ihel * ndiagramgroups + 1211];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1212, &node1211, diagramgroup1212, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1213 = graphNodes[ihel * ndiagramgroups + 1212];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1213, &node1212, diagramgroup1213, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1214 = graphNodes[ihel * ndiagramgroups + 1213];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1214, &node1213, diagramgroup1214, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1215 = graphNodes[ihel * ndiagramgroups + 1214];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1215, &node1214, diagramgroup1215, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1216 = graphNodes[ihel * ndiagramgroups + 1215];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1216, &node1215, diagramgroup1216, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1217 = graphNodes[ihel * ndiagramgroups + 1216];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1217, &node1216, diagramgroup1217, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1218 = graphNodes[ihel * ndiagramgroups + 1217];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1218, &node1217, diagramgroup1218, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1219 = graphNodes[ihel * ndiagramgroups + 1218];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1219, &node1218, diagramgroup1219, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1220 = graphNodes[ihel * ndiagramgroups + 1219];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1220, &node1219, diagramgroup1220, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1221 = graphNodes[ihel * ndiagramgroups + 1220];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1221, &node1220, diagramgroup1221, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1222 = graphNodes[ihel * ndiagramgroups + 1221];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1222, &node1221, diagramgroup1222, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1223 = graphNodes[ihel * ndiagramgroups + 1222];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1223, &node1222, diagramgroup1223, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1224 = graphNodes[ihel * ndiagramgroups + 1223];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1224, &node1223, diagramgroup1224, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1225 = graphNodes[ihel * ndiagramgroups + 1224];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1225, &node1224, diagramgroup1225, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1226 = graphNodes[ihel * ndiagramgroups + 1225];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1226, &node1225, diagramgroup1226, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1227 = graphNodes[ihel * ndiagramgroups + 1226];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1227, &node1226, diagramgroup1227, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1228 = graphNodes[ihel * ndiagramgroups + 1227];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1228, &node1227, diagramgroup1228, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1229 = graphNodes[ihel * ndiagramgroups + 1228];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1229, &node1228, diagramgroup1229, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1230 = graphNodes[ihel * ndiagramgroups + 1229];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1230, &node1229, diagramgroup1230, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1231 = graphNodes[ihel * ndiagramgroups + 1230];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1231, &node1230, diagramgroup1231, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1232 = graphNodes[ihel * ndiagramgroups + 1231];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1232, &node1231, diagramgroup1232, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1233 = graphNodes[ihel * ndiagramgroups + 1232];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1233, &node1232, diagramgroup1233, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1234 = graphNodes[ihel * ndiagramgroups + 1233];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1234, &node1233, diagramgroup1234, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1235 = graphNodes[ihel * ndiagramgroups + 1234];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1235, &node1234, diagramgroup1235, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1236 = graphNodes[ihel * ndiagramgroups + 1235];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1236, &node1235, diagramgroup1236, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1237 = graphNodes[ihel * ndiagramgroups + 1236];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1237, &node1236, diagramgroup1237, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1238 = graphNodes[ihel * ndiagramgroups + 1237];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1238, &node1237, diagramgroup1238, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1239 = graphNodes[ihel * ndiagramgroups + 1238];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1239, &node1238, diagramgroup1239, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1240 = graphNodes[ihel * ndiagramgroups + 1239];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1240, &node1239, diagramgroup1240, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1241 = graphNodes[ihel * ndiagramgroups + 1240];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1241, &node1240, diagramgroup1241, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1242 = graphNodes[ihel * ndiagramgroups + 1241];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1242, &node1241, diagramgroup1242, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1243 = graphNodes[ihel * ndiagramgroups + 1242];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1243, &node1242, diagramgroup1243, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1244 = graphNodes[ihel * ndiagramgroups + 1243];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1244, &node1243, diagramgroup1244, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1245 = graphNodes[ihel * ndiagramgroups + 1244];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1245, &node1244, diagramgroup1245, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1246 = graphNodes[ihel * ndiagramgroups + 1245];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1246, &node1245, diagramgroup1246, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1247 = graphNodes[ihel * ndiagramgroups + 1246];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1247, &node1246, diagramgroup1247, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1248 = graphNodes[ihel * ndiagramgroups + 1247];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1248, &node1247, diagramgroup1248, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1249 = graphNodes[ihel * ndiagramgroups + 1248];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1249, &node1248, diagramgroup1249, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1250 = graphNodes[ihel * ndiagramgroups + 1249];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1250, &node1249, diagramgroup1250, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1251 = graphNodes[ihel * ndiagramgroups + 1250];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1251, &node1250, diagramgroup1251, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1252 = graphNodes[ihel * ndiagramgroups + 1251];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1252, &node1251, diagramgroup1252, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1253 = graphNodes[ihel * ndiagramgroups + 1252];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1253, &node1252, diagramgroup1253, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1254 = graphNodes[ihel * ndiagramgroups + 1253];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1254, &node1253, diagramgroup1254, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1255 = graphNodes[ihel * ndiagramgroups + 1254];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1255, &node1254, diagramgroup1255, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1256 = graphNodes[ihel * ndiagramgroups + 1255];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1256, &node1255, diagramgroup1256, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1257 = graphNodes[ihel * ndiagramgroups + 1256];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1257, &node1256, diagramgroup1257, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1258 = graphNodes[ihel * ndiagramgroups + 1257];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1258, &node1257, diagramgroup1258, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1259 = graphNodes[ihel * ndiagramgroups + 1258];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1259, &node1258, diagramgroup1259, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1260 = graphNodes[ihel * ndiagramgroups + 1259];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1260, &node1259, diagramgroup1260, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1261 = graphNodes[ihel * ndiagramgroups + 1260];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1261, &node1260, diagramgroup1261, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1262 = graphNodes[ihel * ndiagramgroups + 1261];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1262, &node1261, diagramgroup1262, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1263 = graphNodes[ihel * ndiagramgroups + 1262];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1263, &node1262, diagramgroup1263, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1264 = graphNodes[ihel * ndiagramgroups + 1263];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1264, &node1263, diagramgroup1264, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1265 = graphNodes[ihel * ndiagramgroups + 1264];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1265, &node1264, diagramgroup1265, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1266 = graphNodes[ihel * ndiagramgroups + 1265];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1266, &node1265, diagramgroup1266, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1267 = graphNodes[ihel * ndiagramgroups + 1266];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1267, &node1266, diagramgroup1267, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1268 = graphNodes[ihel * ndiagramgroups + 1267];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1268, &node1267, diagramgroup1268, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1269 = graphNodes[ihel * ndiagramgroups + 1268];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1269, &node1268, diagramgroup1269, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1270 = graphNodes[ihel * ndiagramgroups + 1269];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1270, &node1269, diagramgroup1270, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1271 = graphNodes[ihel * ndiagramgroups + 1270];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1271, &node1270, diagramgroup1271, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1272 = graphNodes[ihel * ndiagramgroups + 1271];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1272, &node1271, diagramgroup1272, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1273 = graphNodes[ihel * ndiagramgroups + 1272];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1273, &node1272, diagramgroup1273, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1274 = graphNodes[ihel * ndiagramgroups + 1273];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1274, &node1273, diagramgroup1274, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1275 = graphNodes[ihel * ndiagramgroups + 1274];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1275, &node1274, diagramgroup1275, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1276 = graphNodes[ihel * ndiagramgroups + 1275];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1276, &node1275, diagramgroup1276, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1277 = graphNodes[ihel * ndiagramgroups + 1276];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1277, &node1276, diagramgroup1277, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1278 = graphNodes[ihel * ndiagramgroups + 1277];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1278, &node1277, diagramgroup1278, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1279 = graphNodes[ihel * ndiagramgroups + 1278];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1279, &node1278, diagramgroup1279, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1280 = graphNodes[ihel * ndiagramgroups + 1279];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1280, &node1279, diagramgroup1280, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1281 = graphNodes[ihel * ndiagramgroups + 1280];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1281, &node1280, diagramgroup1281, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1282 = graphNodes[ihel * ndiagramgroups + 1281];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1282, &node1281, diagramgroup1282, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1283 = graphNodes[ihel * ndiagramgroups + 1282];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1283, &node1282, diagramgroup1283, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1284 = graphNodes[ihel * ndiagramgroups + 1283];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1284, &node1283, diagramgroup1284, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1285 = graphNodes[ihel * ndiagramgroups + 1284];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1285, &node1284, diagramgroup1285, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1286 = graphNodes[ihel * ndiagramgroups + 1285];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1286, &node1285, diagramgroup1286, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1287 = graphNodes[ihel * ndiagramgroups + 1286];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1287, &node1286, diagramgroup1287, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1288 = graphNodes[ihel * ndiagramgroups + 1287];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1288, &node1287, diagramgroup1288, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1289 = graphNodes[ihel * ndiagramgroups + 1288];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1289, &node1288, diagramgroup1289, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1290 = graphNodes[ihel * ndiagramgroups + 1289];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1290, &node1289, diagramgroup1290, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1291 = graphNodes[ihel * ndiagramgroups + 1290];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1291, &node1290, diagramgroup1291, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1292 = graphNodes[ihel * ndiagramgroups + 1291];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1292, &node1291, diagramgroup1292, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1293 = graphNodes[ihel * ndiagramgroups + 1292];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1293, &node1292, diagramgroup1293, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1294 = graphNodes[ihel * ndiagramgroups + 1293];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1294, &node1293, diagramgroup1294, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1295 = graphNodes[ihel * ndiagramgroups + 1294];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1295, &node1294, diagramgroup1295, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1296 = graphNodes[ihel * ndiagramgroups + 1295];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1296, &node1295, diagramgroup1296, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1297 = graphNodes[ihel * ndiagramgroups + 1296];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1297, &node1296, diagramgroup1297, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1298 = graphNodes[ihel * ndiagramgroups + 1297];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1298, &node1297, diagramgroup1298, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1299 = graphNodes[ihel * ndiagramgroups + 1298];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1299, &node1298, diagramgroup1299, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1300 = graphNodes[ihel * ndiagramgroups + 1299];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1300, &node1299, diagramgroup1300, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1301 = graphNodes[ihel * ndiagramgroups + 1300];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1301, &node1300, diagramgroup1301, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1302 = graphNodes[ihel * ndiagramgroups + 1301];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1302, &node1301, diagramgroup1302, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1303 = graphNodes[ihel * ndiagramgroups + 1302];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1303, &node1302, diagramgroup1303, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1304 = graphNodes[ihel * ndiagramgroups + 1303];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1304, &node1303, diagramgroup1304, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1305 = graphNodes[ihel * ndiagramgroups + 1304];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1305, &node1304, diagramgroup1305, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1306 = graphNodes[ihel * ndiagramgroups + 1305];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1306, &node1305, diagramgroup1306, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1307 = graphNodes[ihel * ndiagramgroups + 1306];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1307, &node1306, diagramgroup1307, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1308 = graphNodes[ihel * ndiagramgroups + 1307];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1308, &node1307, diagramgroup1308, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1309 = graphNodes[ihel * ndiagramgroups + 1308];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1309, &node1308, diagramgroup1309, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1310 = graphNodes[ihel * ndiagramgroups + 1309];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1310, &node1309, diagramgroup1310, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1311 = graphNodes[ihel * ndiagramgroups + 1310];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1311, &node1310, diagramgroup1311, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1312 = graphNodes[ihel * ndiagramgroups + 1311];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1312, &node1311, diagramgroup1312, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1313 = graphNodes[ihel * ndiagramgroups + 1312];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1313, &node1312, diagramgroup1313, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1314 = graphNodes[ihel * ndiagramgroups + 1313];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1314, &node1313, diagramgroup1314, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1315 = graphNodes[ihel * ndiagramgroups + 1314];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1315, &node1314, diagramgroup1315, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1316 = graphNodes[ihel * ndiagramgroups + 1315];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1316, &node1315, diagramgroup1316, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1317 = graphNodes[ihel * ndiagramgroups + 1316];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1317, &node1316, diagramgroup1317, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1318 = graphNodes[ihel * ndiagramgroups + 1317];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1318, &node1317, diagramgroup1318, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1319 = graphNodes[ihel * ndiagramgroups + 1318];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1319, &node1318, diagramgroup1319, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1320 = graphNodes[ihel * ndiagramgroups + 1319];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1320, &node1319, diagramgroup1320, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1321 = graphNodes[ihel * ndiagramgroups + 1320];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1321, &node1320, diagramgroup1321, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1322 = graphNodes[ihel * ndiagramgroups + 1321];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1322, &node1321, diagramgroup1322, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1323 = graphNodes[ihel * ndiagramgroups + 1322];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1323, &node1322, diagramgroup1323, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1324 = graphNodes[ihel * ndiagramgroups + 1323];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1324, &node1323, diagramgroup1324, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1325 = graphNodes[ihel * ndiagramgroups + 1324];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1325, &node1324, diagramgroup1325, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1326 = graphNodes[ihel * ndiagramgroups + 1325];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1326, &node1325, diagramgroup1326, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1327 = graphNodes[ihel * ndiagramgroups + 1326];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1327, &node1326, diagramgroup1327, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1328 = graphNodes[ihel * ndiagramgroups + 1327];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1328, &node1327, diagramgroup1328, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1329 = graphNodes[ihel * ndiagramgroups + 1328];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1329, &node1328, diagramgroup1329, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1330 = graphNodes[ihel * ndiagramgroups + 1329];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1330, &node1329, diagramgroup1330, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1331 = graphNodes[ihel * ndiagramgroups + 1330];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1331, &node1330, diagramgroup1331, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1332 = graphNodes[ihel * ndiagramgroups + 1331];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1332, &node1331, diagramgroup1332, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1333 = graphNodes[ihel * ndiagramgroups + 1332];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1333, &node1332, diagramgroup1333, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1334 = graphNodes[ihel * ndiagramgroups + 1333];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1334, &node1333, diagramgroup1334, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1335 = graphNodes[ihel * ndiagramgroups + 1334];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1335, &node1334, diagramgroup1335, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1336 = graphNodes[ihel * ndiagramgroups + 1335];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1336, &node1335, diagramgroup1336, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1337 = graphNodes[ihel * ndiagramgroups + 1336];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1337, &node1336, diagramgroup1337, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1338 = graphNodes[ihel * ndiagramgroups + 1337];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1338, &node1337, diagramgroup1338, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1339 = graphNodes[ihel * ndiagramgroups + 1338];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1339, &node1338, diagramgroup1339, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1340 = graphNodes[ihel * ndiagramgroups + 1339];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1340, &node1339, diagramgroup1340, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1341 = graphNodes[ihel * ndiagramgroups + 1340];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1341, &node1340, diagramgroup1341, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1342 = graphNodes[ihel * ndiagramgroups + 1341];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1342, &node1341, diagramgroup1342, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1343 = graphNodes[ihel * ndiagramgroups + 1342];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1343, &node1342, diagramgroup1343, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1344 = graphNodes[ihel * ndiagramgroups + 1343];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1344, &node1343, diagramgroup1344, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1345 = graphNodes[ihel * ndiagramgroups + 1344];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1345, &node1344, diagramgroup1345, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1346 = graphNodes[ihel * ndiagramgroups + 1345];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1346, &node1345, diagramgroup1346, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1347 = graphNodes[ihel * ndiagramgroups + 1346];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1347, &node1346, diagramgroup1347, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1348 = graphNodes[ihel * ndiagramgroups + 1347];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1348, &node1347, diagramgroup1348, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1349 = graphNodes[ihel * ndiagramgroups + 1348];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1349, &node1348, diagramgroup1349, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1350 = graphNodes[ihel * ndiagramgroups + 1349];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1350, &node1349, diagramgroup1350, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1351 = graphNodes[ihel * ndiagramgroups + 1350];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1351, &node1350, diagramgroup1351, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1352 = graphNodes[ihel * ndiagramgroups + 1351];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1352, &node1351, diagramgroup1352, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1353 = graphNodes[ihel * ndiagramgroups + 1352];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1353, &node1352, diagramgroup1353, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1354 = graphNodes[ihel * ndiagramgroups + 1353];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1354, &node1353, diagramgroup1354, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1355 = graphNodes[ihel * ndiagramgroups + 1354];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1355, &node1354, diagramgroup1355, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1356 = graphNodes[ihel * ndiagramgroups + 1355];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1356, &node1355, diagramgroup1356, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1357 = graphNodes[ihel * ndiagramgroups + 1356];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1357, &node1356, diagramgroup1357, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1358 = graphNodes[ihel * ndiagramgroups + 1357];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1358, &node1357, diagramgroup1358, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1359 = graphNodes[ihel * ndiagramgroups + 1358];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1359, &node1358, diagramgroup1359, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1360 = graphNodes[ihel * ndiagramgroups + 1359];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1360, &node1359, diagramgroup1360, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1361 = graphNodes[ihel * ndiagramgroups + 1360];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1361, &node1360, diagramgroup1361, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1362 = graphNodes[ihel * ndiagramgroups + 1361];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1362, &node1361, diagramgroup1362, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1363 = graphNodes[ihel * ndiagramgroups + 1362];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1363, &node1362, diagramgroup1363, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1364 = graphNodes[ihel * ndiagramgroups + 1363];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1364, &node1363, diagramgroup1364, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1365 = graphNodes[ihel * ndiagramgroups + 1364];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1365, &node1364, diagramgroup1365, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1366 = graphNodes[ihel * ndiagramgroups + 1365];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1366, &node1365, diagramgroup1366, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1367 = graphNodes[ihel * ndiagramgroups + 1366];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1367, &node1366, diagramgroup1367, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1368 = graphNodes[ihel * ndiagramgroups + 1367];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1368, &node1367, diagramgroup1368, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1369 = graphNodes[ihel * ndiagramgroups + 1368];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1369, &node1368, diagramgroup1369, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1370 = graphNodes[ihel * ndiagramgroups + 1369];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1370, &node1369, diagramgroup1370, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1371 = graphNodes[ihel * ndiagramgroups + 1370];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1371, &node1370, diagramgroup1371, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1372 = graphNodes[ihel * ndiagramgroups + 1371];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1372, &node1371, diagramgroup1372, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1373 = graphNodes[ihel * ndiagramgroups + 1372];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1373, &node1372, diagramgroup1373, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1374 = graphNodes[ihel * ndiagramgroups + 1373];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1374, &node1373, diagramgroup1374, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1375 = graphNodes[ihel * ndiagramgroups + 1374];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1375, &node1374, diagramgroup1375, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1376 = graphNodes[ihel * ndiagramgroups + 1375];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1376, &node1375, diagramgroup1376, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1377 = graphNodes[ihel * ndiagramgroups + 1376];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1377, &node1376, diagramgroup1377, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1378 = graphNodes[ihel * ndiagramgroups + 1377];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1378, &node1377, diagramgroup1378, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1379 = graphNodes[ihel * ndiagramgroups + 1378];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1379, &node1378, diagramgroup1379, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1380 = graphNodes[ihel * ndiagramgroups + 1379];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1380, &node1379, diagramgroup1380, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1381 = graphNodes[ihel * ndiagramgroups + 1380];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1381, &node1380, diagramgroup1381, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1382 = graphNodes[ihel * ndiagramgroups + 1381];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1382, &node1381, diagramgroup1382, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1383 = graphNodes[ihel * ndiagramgroups + 1382];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1383, &node1382, diagramgroup1383, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1384 = graphNodes[ihel * ndiagramgroups + 1383];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1384, &node1383, diagramgroup1384, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1385 = graphNodes[ihel * ndiagramgroups + 1384];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1385, &node1384, diagramgroup1385, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1386 = graphNodes[ihel * ndiagramgroups + 1385];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1386, &node1385, diagramgroup1386, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1387 = graphNodes[ihel * ndiagramgroups + 1386];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1387, &node1386, diagramgroup1387, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1388 = graphNodes[ihel * ndiagramgroups + 1387];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1388, &node1387, diagramgroup1388, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1389 = graphNodes[ihel * ndiagramgroups + 1388];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1389, &node1388, diagramgroup1389, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1390 = graphNodes[ihel * ndiagramgroups + 1389];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1390, &node1389, diagramgroup1390, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1391 = graphNodes[ihel * ndiagramgroups + 1390];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1391, &node1390, diagramgroup1391, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1392 = graphNodes[ihel * ndiagramgroups + 1391];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1392, &node1391, diagramgroup1392, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1393 = graphNodes[ihel * ndiagramgroups + 1392];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1393, &node1392, diagramgroup1393, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1394 = graphNodes[ihel * ndiagramgroups + 1393];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1394, &node1393, diagramgroup1394, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1395 = graphNodes[ihel * ndiagramgroups + 1394];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1395, &node1394, diagramgroup1395, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1396 = graphNodes[ihel * ndiagramgroups + 1395];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1396, &node1395, diagramgroup1396, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1397 = graphNodes[ihel * ndiagramgroups + 1396];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1397, &node1396, diagramgroup1397, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1398 = graphNodes[ihel * ndiagramgroups + 1397];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1398, &node1397, diagramgroup1398, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1399 = graphNodes[ihel * ndiagramgroups + 1398];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1399, &node1398, diagramgroup1399, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1400 = graphNodes[ihel * ndiagramgroups + 1399];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1400, &node1399, diagramgroup1400, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1401 = graphNodes[ihel * ndiagramgroups + 1400];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1401, &node1400, diagramgroup1401, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1402 = graphNodes[ihel * ndiagramgroups + 1401];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1402, &node1401, diagramgroup1402, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1403 = graphNodes[ihel * ndiagramgroups + 1402];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1403, &node1402, diagramgroup1403, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1404 = graphNodes[ihel * ndiagramgroups + 1403];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1404, &node1403, diagramgroup1404, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1405 = graphNodes[ihel * ndiagramgroups + 1404];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1405, &node1404, diagramgroup1405, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1406 = graphNodes[ihel * ndiagramgroups + 1405];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1406, &node1405, diagramgroup1406, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1407 = graphNodes[ihel * ndiagramgroups + 1406];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1407, &node1406, diagramgroup1407, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1408 = graphNodes[ihel * ndiagramgroups + 1407];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1408, &node1407, diagramgroup1408, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1409 = graphNodes[ihel * ndiagramgroups + 1408];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1409, &node1408, diagramgroup1409, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1410 = graphNodes[ihel * ndiagramgroups + 1409];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1410, &node1409, diagramgroup1410, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1411 = graphNodes[ihel * ndiagramgroups + 1410];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1411, &node1410, diagramgroup1411, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1412 = graphNodes[ihel * ndiagramgroups + 1411];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1412, &node1411, diagramgroup1412, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1413 = graphNodes[ihel * ndiagramgroups + 1412];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1413, &node1412, diagramgroup1413, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1414 = graphNodes[ihel * ndiagramgroups + 1413];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1414, &node1413, diagramgroup1414, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1415 = graphNodes[ihel * ndiagramgroups + 1414];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1415, &node1414, diagramgroup1415, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1416 = graphNodes[ihel * ndiagramgroups + 1415];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1416, &node1415, diagramgroup1416, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1417 = graphNodes[ihel * ndiagramgroups + 1416];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1417, &node1416, diagramgroup1417, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1418 = graphNodes[ihel * ndiagramgroups + 1417];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1418, &node1417, diagramgroup1418, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1419 = graphNodes[ihel * ndiagramgroups + 1418];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1419, &node1418, diagramgroup1419, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1420 = graphNodes[ihel * ndiagramgroups + 1419];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1420, &node1419, diagramgroup1420, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1421 = graphNodes[ihel * ndiagramgroups + 1420];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1421, &node1420, diagramgroup1421, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1422 = graphNodes[ihel * ndiagramgroups + 1421];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1422, &node1421, diagramgroup1422, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1423 = graphNodes[ihel * ndiagramgroups + 1422];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1423, &node1422, diagramgroup1423, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1424 = graphNodes[ihel * ndiagramgroups + 1423];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1424, &node1423, diagramgroup1424, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1425 = graphNodes[ihel * ndiagramgroups + 1424];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1425, &node1424, diagramgroup1425, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1426 = graphNodes[ihel * ndiagramgroups + 1425];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1426, &node1425, diagramgroup1426, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1427 = graphNodes[ihel * ndiagramgroups + 1426];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1427, &node1426, diagramgroup1427, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1428 = graphNodes[ihel * ndiagramgroups + 1427];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1428, &node1427, diagramgroup1428, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1429 = graphNodes[ihel * ndiagramgroups + 1428];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1429, &node1428, diagramgroup1429, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1430 = graphNodes[ihel * ndiagramgroups + 1429];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1430, &node1429, diagramgroup1430, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1431 = graphNodes[ihel * ndiagramgroups + 1430];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1431, &node1430, diagramgroup1431, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1432 = graphNodes[ihel * ndiagramgroups + 1431];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1432, &node1431, diagramgroup1432, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1433 = graphNodes[ihel * ndiagramgroups + 1432];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1433, &node1432, diagramgroup1433, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1434 = graphNodes[ihel * ndiagramgroups + 1433];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1434, &node1433, diagramgroup1434, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1435 = graphNodes[ihel * ndiagramgroups + 1434];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1435, &node1434, diagramgroup1435, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1436 = graphNodes[ihel * ndiagramgroups + 1435];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1436, &node1435, diagramgroup1436, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1437 = graphNodes[ihel * ndiagramgroups + 1436];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1437, &node1436, diagramgroup1437, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1438 = graphNodes[ihel * ndiagramgroups + 1437];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1438, &node1437, diagramgroup1438, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1439 = graphNodes[ihel * ndiagramgroups + 1438];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1439, &node1438, diagramgroup1439, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1440 = graphNodes[ihel * ndiagramgroups + 1439];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1440, &node1439, diagramgroup1440, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1441 = graphNodes[ihel * ndiagramgroups + 1440];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1441, &node1440, diagramgroup1441, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1442 = graphNodes[ihel * ndiagramgroups + 1441];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1442, &node1441, diagramgroup1442, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1443 = graphNodes[ihel * ndiagramgroups + 1442];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1443, &node1442, diagramgroup1443, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1444 = graphNodes[ihel * ndiagramgroups + 1443];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1444, &node1443, diagramgroup1444, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1445 = graphNodes[ihel * ndiagramgroups + 1444];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1445, &node1444, diagramgroup1445, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1446 = graphNodes[ihel * ndiagramgroups + 1445];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1446, &node1445, diagramgroup1446, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1447 = graphNodes[ihel * ndiagramgroups + 1446];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1447, &node1446, diagramgroup1447, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1448 = graphNodes[ihel * ndiagramgroups + 1447];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1448, &node1447, diagramgroup1448, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1449 = graphNodes[ihel * ndiagramgroups + 1448];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1449, &node1448, diagramgroup1449, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1450 = graphNodes[ihel * ndiagramgroups + 1449];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1450, &node1449, diagramgroup1450, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1451 = graphNodes[ihel * ndiagramgroups + 1450];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1451, &node1450, diagramgroup1451, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1452 = graphNodes[ihel * ndiagramgroups + 1451];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1452, &node1451, diagramgroup1452, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1453 = graphNodes[ihel * ndiagramgroups + 1452];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1453, &node1452, diagramgroup1453, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1454 = graphNodes[ihel * ndiagramgroups + 1453];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1454, &node1453, diagramgroup1454, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1455 = graphNodes[ihel * ndiagramgroups + 1454];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1455, &node1454, diagramgroup1455, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1456 = graphNodes[ihel * ndiagramgroups + 1455];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1456, &node1455, diagramgroup1456, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1457 = graphNodes[ihel * ndiagramgroups + 1456];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1457, &node1456, diagramgroup1457, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1458 = graphNodes[ihel * ndiagramgroups + 1457];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1458, &node1457, diagramgroup1458, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1459 = graphNodes[ihel * ndiagramgroups + 1458];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1459, &node1458, diagramgroup1459, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1460 = graphNodes[ihel * ndiagramgroups + 1459];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1460, &node1459, diagramgroup1460, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1461 = graphNodes[ihel * ndiagramgroups + 1460];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1461, &node1460, diagramgroup1461, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1462 = graphNodes[ihel * ndiagramgroups + 1461];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1462, &node1461, diagramgroup1462, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1463 = graphNodes[ihel * ndiagramgroups + 1462];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1463, &node1462, diagramgroup1463, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1464 = graphNodes[ihel * ndiagramgroups + 1463];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1464, &node1463, diagramgroup1464, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1465 = graphNodes[ihel * ndiagramgroups + 1464];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1465, &node1464, diagramgroup1465, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1466 = graphNodes[ihel * ndiagramgroups + 1465];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1466, &node1465, diagramgroup1466, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1467 = graphNodes[ihel * ndiagramgroups + 1466];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1467, &node1466, diagramgroup1467, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1468 = graphNodes[ihel * ndiagramgroups + 1467];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1468, &node1467, diagramgroup1468, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1469 = graphNodes[ihel * ndiagramgroups + 1468];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1469, &node1468, diagramgroup1469, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1470 = graphNodes[ihel * ndiagramgroups + 1469];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1470, &node1469, diagramgroup1470, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1471 = graphNodes[ihel * ndiagramgroups + 1470];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1471, &node1470, diagramgroup1471, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1472 = graphNodes[ihel * ndiagramgroups + 1471];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1472, &node1471, diagramgroup1472, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1473 = graphNodes[ihel * ndiagramgroups + 1472];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1473, &node1472, diagramgroup1473, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1474 = graphNodes[ihel * ndiagramgroups + 1473];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1474, &node1473, diagramgroup1474, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1475 = graphNodes[ihel * ndiagramgroups + 1474];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1475, &node1474, diagramgroup1475, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1476 = graphNodes[ihel * ndiagramgroups + 1475];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1476, &node1475, diagramgroup1476, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1477 = graphNodes[ihel * ndiagramgroups + 1476];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1477, &node1476, diagramgroup1477, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1478 = graphNodes[ihel * ndiagramgroups + 1477];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1478, &node1477, diagramgroup1478, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1479 = graphNodes[ihel * ndiagramgroups + 1478];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1479, &node1478, diagramgroup1479, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1480 = graphNodes[ihel * ndiagramgroups + 1479];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1480, &node1479, diagramgroup1480, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1481 = graphNodes[ihel * ndiagramgroups + 1480];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1481, &node1480, diagramgroup1481, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1482 = graphNodes[ihel * ndiagramgroups + 1481];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1482, &node1481, diagramgroup1482, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1483 = graphNodes[ihel * ndiagramgroups + 1482];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1483, &node1482, diagramgroup1483, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1484 = graphNodes[ihel * ndiagramgroups + 1483];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1484, &node1483, diagramgroup1484, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1485 = graphNodes[ihel * ndiagramgroups + 1484];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1485, &node1484, diagramgroup1485, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1486 = graphNodes[ihel * ndiagramgroups + 1485];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1486, &node1485, diagramgroup1486, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1487 = graphNodes[ihel * ndiagramgroups + 1486];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1487, &node1486, diagramgroup1487, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1488 = graphNodes[ihel * ndiagramgroups + 1487];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1488, &node1487, diagramgroup1488, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1489 = graphNodes[ihel * ndiagramgroups + 1488];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1489, &node1488, diagramgroup1489, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1490 = graphNodes[ihel * ndiagramgroups + 1489];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1490, &node1489, diagramgroup1490, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1491 = graphNodes[ihel * ndiagramgroups + 1490];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1491, &node1490, diagramgroup1491, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1492 = graphNodes[ihel * ndiagramgroups + 1491];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1492, &node1491, diagramgroup1492, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1493 = graphNodes[ihel * ndiagramgroups + 1492];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1493, &node1492, diagramgroup1493, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1494 = graphNodes[ihel * ndiagramgroups + 1493];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1494, &node1493, diagramgroup1494, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1495 = graphNodes[ihel * ndiagramgroups + 1494];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1495, &node1494, diagramgroup1495, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1496 = graphNodes[ihel * ndiagramgroups + 1495];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1496, &node1495, diagramgroup1496, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1497 = graphNodes[ihel * ndiagramgroups + 1496];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1497, &node1496, diagramgroup1497, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1498 = graphNodes[ihel * ndiagramgroups + 1497];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1498, &node1497, diagramgroup1498, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1499 = graphNodes[ihel * ndiagramgroups + 1498];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1499, &node1498, diagramgroup1499, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1500 = graphNodes[ihel * ndiagramgroups + 1499];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1500, &node1499, diagramgroup1500, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1501 = graphNodes[ihel * ndiagramgroups + 1500];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1501, &node1500, diagramgroup1501, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1502 = graphNodes[ihel * ndiagramgroups + 1501];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1502, &node1501, diagramgroup1502, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1503 = graphNodes[ihel * ndiagramgroups + 1502];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1503, &node1502, diagramgroup1503, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1504 = graphNodes[ihel * ndiagramgroups + 1503];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1504, &node1503, diagramgroup1504, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1505 = graphNodes[ihel * ndiagramgroups + 1504];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1505, &node1504, diagramgroup1505, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1506 = graphNodes[ihel * ndiagramgroups + 1505];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1506, &node1505, diagramgroup1506, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1507 = graphNodes[ihel * ndiagramgroups + 1506];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1507, &node1506, diagramgroup1507, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1508 = graphNodes[ihel * ndiagramgroups + 1507];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1508, &node1507, diagramgroup1508, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1509 = graphNodes[ihel * ndiagramgroups + 1508];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1509, &node1508, diagramgroup1509, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1510 = graphNodes[ihel * ndiagramgroups + 1509];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1510, &node1509, diagramgroup1510, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1511 = graphNodes[ihel * ndiagramgroups + 1510];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1511, &node1510, diagramgroup1511, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1512 = graphNodes[ihel * ndiagramgroups + 1511];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1512, &node1511, diagramgroup1512, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1513 = graphNodes[ihel * ndiagramgroups + 1512];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1513, &node1512, diagramgroup1513, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1514 = graphNodes[ihel * ndiagramgroups + 1513];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1514, &node1513, diagramgroup1514, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1515 = graphNodes[ihel * ndiagramgroups + 1514];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1515, &node1514, diagramgroup1515, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1516 = graphNodes[ihel * ndiagramgroups + 1515];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1516, &node1515, diagramgroup1516, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1517 = graphNodes[ihel * ndiagramgroups + 1516];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1517, &node1516, diagramgroup1517, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1518 = graphNodes[ihel * ndiagramgroups + 1517];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1518, &node1517, diagramgroup1518, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1519 = graphNodes[ihel * ndiagramgroups + 1518];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1519, &node1518, diagramgroup1519, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1520 = graphNodes[ihel * ndiagramgroups + 1519];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1520, &node1519, diagramgroup1520, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1521 = graphNodes[ihel * ndiagramgroups + 1520];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1521, &node1520, diagramgroup1521, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1522 = graphNodes[ihel * ndiagramgroups + 1521];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1522, &node1521, diagramgroup1522, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1523 = graphNodes[ihel * ndiagramgroups + 1522];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1523, &node1522, diagramgroup1523, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1524 = graphNodes[ihel * ndiagramgroups + 1523];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1524, &node1523, diagramgroup1524, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1525 = graphNodes[ihel * ndiagramgroups + 1524];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1525, &node1524, diagramgroup1525, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1526 = graphNodes[ihel * ndiagramgroups + 1525];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1526, &node1525, diagramgroup1526, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1527 = graphNodes[ihel * ndiagramgroups + 1526];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1527, &node1526, diagramgroup1527, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1528 = graphNodes[ihel * ndiagramgroups + 1527];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1528, &node1527, diagramgroup1528, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1529 = graphNodes[ihel * ndiagramgroups + 1528];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1529, &node1528, diagramgroup1529, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1530 = graphNodes[ihel * ndiagramgroups + 1529];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1530, &node1529, diagramgroup1530, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1531 = graphNodes[ihel * ndiagramgroups + 1530];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1531, &node1530, diagramgroup1531, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1532 = graphNodes[ihel * ndiagramgroups + 1531];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1532, &node1531, diagramgroup1532, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1533 = graphNodes[ihel * ndiagramgroups + 1532];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1533, &node1532, diagramgroup1533, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1534 = graphNodes[ihel * ndiagramgroups + 1533];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1534, &node1533, diagramgroup1534, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1535 = graphNodes[ihel * ndiagramgroups + 1534];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1535, &node1534, diagramgroup1535, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1536 = graphNodes[ihel * ndiagramgroups + 1535];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1536, &node1535, diagramgroup1536, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1537 = graphNodes[ihel * ndiagramgroups + 1536];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1537, &node1536, diagramgroup1537, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1538 = graphNodes[ihel * ndiagramgroups + 1537];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1538, &node1537, diagramgroup1538, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1539 = graphNodes[ihel * ndiagramgroups + 1538];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1539, &node1538, diagramgroup1539, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1540 = graphNodes[ihel * ndiagramgroups + 1539];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1540, &node1539, diagramgroup1540, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1541 = graphNodes[ihel * ndiagramgroups + 1540];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1541, &node1540, diagramgroup1541, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1542 = graphNodes[ihel * ndiagramgroups + 1541];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1542, &node1541, diagramgroup1542, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1543 = graphNodes[ihel * ndiagramgroups + 1542];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1543, &node1542, diagramgroup1543, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1544 = graphNodes[ihel * ndiagramgroups + 1543];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1544, &node1543, diagramgroup1544, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1545 = graphNodes[ihel * ndiagramgroups + 1544];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1545, &node1544, diagramgroup1545, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1546 = graphNodes[ihel * ndiagramgroups + 1545];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1546, &node1545, diagramgroup1546, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1547 = graphNodes[ihel * ndiagramgroups + 1546];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1547, &node1546, diagramgroup1547, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1548 = graphNodes[ihel * ndiagramgroups + 1547];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1548, &node1547, diagramgroup1548, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1549 = graphNodes[ihel * ndiagramgroups + 1548];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1549, &node1548, diagramgroup1549, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node1550 = graphNodes[ihel * ndiagramgroups + 1549];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1550, &node1549, diagramgroup1550, gpublocks, gputhreads, gpustream, wfs, jamps, cNGoodHel, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      // Case 1 with graphs (gpustream!=0, sigmaKin): create the graph executor if not yet done, then launch the graph executor
      if( useGraphs && gpustream != 0 )
      {
        if( !graphExec )
        {
          checkGpu( gpuGraphInstantiate( &graphExec, graph, nullptr, nullptr, 0 ) );
          //std::cout << "(ihel=" << ihel << ") Created graph executor " << &graphExec << " for graph " << graph << std::endl;
        }
        //std::cout << "(ihel=" << ihel << ") Launch graph executor " << &graphExec << " for graph " << graph << std::endl;
        checkGpu( gpuGraphLaunch( graphExec, gpustream ) );
      }
#else
      // === GPU IMPLEMENTATION (DCDIAG=1): merge all diagram groups into a single kernel ===
      diagramgroup1( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD, dcHel, momenta, ihel );
      diagramgroup2( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup3( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup4( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup5( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup6( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup7( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup8( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup9( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup10( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup11( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup12( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup13( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup14( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup15( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup16( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup17( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup18( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup19( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup20( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup21( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup22( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup23( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup24( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup25( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup26( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup27( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup28( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup29( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup30( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup31( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup32( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup33( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup34( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup35( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup36( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup37( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup38( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup39( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup40( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup41( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup42( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup43( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup44( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup45( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup46( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup47( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup48( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup49( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup50( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup51( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup52( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup53( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup54( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup55( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup56( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup57( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup58( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup59( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup60( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup61( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup62( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup63( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup64( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup65( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup66( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup67( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup68( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup69( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup70( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup71( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup72( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup73( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup74( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup75( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup76( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup77( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup78( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup79( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup80( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup81( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup82( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup83( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup84( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup85( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup86( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup87( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup88( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup89( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup90( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup91( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup92( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup93( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup94( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup95( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup96( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup97( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup98( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup99( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup100( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup101( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup102( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup103( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup104( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup105( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup106( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup107( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup108( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup109( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup110( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup111( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup112( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup113( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup114( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup115( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup116( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup117( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup118( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup119( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup120( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup121( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup122( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup123( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup124( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup125( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup126( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup127( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup128( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup129( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup130( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup131( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup132( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup133( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup134( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup135( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup136( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup137( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup138( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup139( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup140( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup141( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup142( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup143( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup144( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup145( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup146( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup147( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup148( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup149( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup150( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup151( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup152( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup153( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup154( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup155( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup156( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup157( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup158( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup159( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup160( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup161( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup162( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup163( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup164( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup165( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup166( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup167( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup168( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup169( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup170( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup171( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup172( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup173( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup174( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup175( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup176( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup177( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup178( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup179( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup180( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup181( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup182( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup183( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup184( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup185( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup186( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup187( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup188( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup189( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup190( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup191( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup192( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup193( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup194( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup195( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup196( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup197( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup198( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup199( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup200( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup201( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup202( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup203( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup204( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup205( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup206( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup207( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup208( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup209( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup210( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup211( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup212( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup213( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup214( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup215( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup216( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup217( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup218( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup219( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup220( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup221( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup222( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup223( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup224( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup225( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup226( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup227( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup228( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup229( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup230( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup231( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup232( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup233( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup234( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup235( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup236( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup237( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup238( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup239( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup240( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup241( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup242( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup243( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup244( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup245( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup246( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup247( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup248( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup249( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup250( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup251( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup252( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup253( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup254( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup255( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup256( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup257( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup258( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup259( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup260( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup261( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup262( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup263( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup264( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup265( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup266( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup267( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup268( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup269( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup270( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup271( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup272( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup273( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup274( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup275( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup276( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup277( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup278( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup279( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup280( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup281( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup282( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup283( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup284( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup285( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup286( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup287( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup288( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup289( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup290( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup291( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup292( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup293( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup294( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup295( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup296( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup297( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup298( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup299( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup300( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup301( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup302( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup303( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup304( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup305( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup306( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup307( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup308( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup309( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup310( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup311( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup312( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup313( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup314( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup315( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup316( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup317( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup318( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup319( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup320( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup321( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup322( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup323( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup324( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup325( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup326( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup327( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup328( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup329( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup330( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup331( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup332( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup333( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup334( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup335( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup336( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup337( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup338( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup339( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup340( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup341( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup342( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup343( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup344( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup345( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup346( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup347( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup348( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup349( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup350( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup351( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup352( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup353( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup354( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup355( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup356( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup357( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup358( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup359( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup360( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup361( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup362( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup363( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup364( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup365( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup366( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup367( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup368( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup369( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup370( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup371( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup372( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup373( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup374( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup375( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup376( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup377( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup378( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup379( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup380( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup381( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup382( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup383( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup384( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup385( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup386( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup387( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup388( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup389( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup390( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup391( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup392( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup393( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup394( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup395( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup396( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup397( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup398( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup399( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup400( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup401( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup402( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup403( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup404( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup405( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup406( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup407( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup408( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup409( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup410( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup411( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup412( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup413( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup414( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup415( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup416( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup417( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup418( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup419( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup420( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup421( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup422( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup423( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup424( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup425( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup426( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup427( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup428( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup429( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup430( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup431( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup432( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup433( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup434( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup435( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup436( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup437( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup438( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup439( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup440( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup441( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup442( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup443( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup444( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup445( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup446( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup447( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup448( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup449( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup450( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup451( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup452( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup453( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup454( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup455( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup456( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup457( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup458( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup459( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup460( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup461( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup462( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup463( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup464( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup465( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup466( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup467( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup468( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup469( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup470( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup471( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup472( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup473( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup474( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup475( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup476( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup477( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup478( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup479( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup480( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup481( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup482( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup483( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup484( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup485( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup486( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup487( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup488( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup489( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup490( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup491( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup492( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup493( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup494( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup495( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup496( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup497( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup498( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup499( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup500( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup501( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup502( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup503( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup504( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup505( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup506( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup507( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup508( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup509( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup510( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup511( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup512( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup513( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup514( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup515( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup516( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup517( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup518( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup519( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup520( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup521( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup522( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup523( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup524( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup525( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup526( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup527( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup528( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup529( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup530( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup531( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup532( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup533( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup534( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup535( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup536( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup537( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup538( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup539( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup540( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup541( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup542( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup543( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup544( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup545( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup546( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup547( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup548( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup549( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup550( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup551( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup552( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup553( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup554( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup555( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup556( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup557( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup558( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup559( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup560( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup561( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup562( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup563( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup564( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup565( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup566( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup567( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup568( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup569( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup570( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup571( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup572( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup573( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup574( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup575( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup576( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup577( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup578( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup579( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup580( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup581( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup582( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup583( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup584( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup585( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup586( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup587( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup588( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup589( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup590( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup591( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup592( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup593( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup594( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup595( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup596( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup597( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup598( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup599( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup600( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup601( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup602( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup603( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup604( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup605( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup606( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup607( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup608( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup609( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup610( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup611( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup612( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup613( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup614( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup615( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup616( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup617( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup618( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup619( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup620( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup621( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup622( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup623( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup624( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup625( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup626( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup627( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup628( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup629( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup630( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup631( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup632( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup633( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup634( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup635( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup636( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup637( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup638( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup639( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup640( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup641( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup642( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup643( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup644( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup645( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup646( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup647( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup648( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup649( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup650( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup651( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup652( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup653( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup654( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup655( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup656( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup657( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup658( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup659( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup660( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup661( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup662( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup663( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup664( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup665( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup666( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup667( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup668( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup669( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup670( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup671( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup672( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup673( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup674( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup675( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup676( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup677( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup678( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup679( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup680( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup681( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup682( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup683( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup684( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup685( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup686( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup687( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup688( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup689( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup690( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup691( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup692( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup693( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup694( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup695( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup696( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup697( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup698( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup699( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup700( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup701( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup702( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup703( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup704( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup705( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup706( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup707( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup708( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup709( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup710( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup711( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup712( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup713( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup714( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup715( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup716( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup717( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup718( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup719( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup720( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup721( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup722( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup723( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup724( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup725( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup726( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup727( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup728( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup729( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup730( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup731( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup732( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup733( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup734( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup735( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup736( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup737( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup738( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup739( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup740( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup741( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup742( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup743( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup744( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup745( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup746( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup747( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup748( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup749( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup750( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup751( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup752( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup753( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup754( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup755( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup756( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup757( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup758( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup759( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup760( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup761( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup762( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup763( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup764( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup765( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup766( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup767( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup768( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup769( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup770( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup771( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup772( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup773( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup774( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup775( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup776( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup777( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup778( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup779( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup780( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup781( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup782( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup783( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup784( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup785( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup786( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup787( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup788( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup789( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup790( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup791( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup792( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup793( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup794( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup795( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup796( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup797( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup798( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup799( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup800( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup801( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup802( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup803( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup804( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup805( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup806( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup807( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup808( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup809( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup810( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup811( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup812( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup813( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup814( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup815( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup816( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup817( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup818( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup819( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup820( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup821( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup822( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup823( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup824( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup825( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup826( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup827( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup828( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup829( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup830( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup831( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup832( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup833( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup834( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup835( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup836( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup837( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup838( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup839( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup840( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup841( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup842( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup843( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup844( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup845( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup846( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup847( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup848( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup849( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup850( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup851( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup852( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup853( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup854( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup855( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup856( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup857( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup858( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup859( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup860( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup861( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup862( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup863( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup864( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup865( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup866( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup867( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup868( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup869( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup870( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup871( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup872( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup873( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup874( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup875( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup876( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup877( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup878( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup879( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup880( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup881( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup882( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup883( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup884( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup885( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup886( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup887( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup888( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup889( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup890( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup891( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup892( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup893( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup894( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup895( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup896( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup897( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup898( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup899( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup900( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup901( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup902( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup903( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup904( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup905( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup906( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup907( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup908( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup909( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup910( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup911( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup912( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup913( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup914( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup915( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup916( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup917( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup918( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup919( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup920( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup921( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup922( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup923( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup924( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup925( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup926( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup927( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup928( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup929( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup930( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup931( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup932( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup933( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup934( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup935( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup936( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup937( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup938( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup939( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup940( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup941( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup942( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup943( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup944( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup945( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup946( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup947( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup948( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup949( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup950( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup951( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup952( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup953( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup954( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup955( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup956( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup957( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup958( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup959( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup960( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup961( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup962( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup963( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup964( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup965( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup966( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup967( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup968( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup969( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup970( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup971( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup972( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup973( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup974( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup975( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup976( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup977( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup978( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup979( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup980( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup981( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup982( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup983( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup984( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup985( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup986( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup987( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup988( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup989( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup990( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup991( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup992( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup993( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup994( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup995( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup996( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup997( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup998( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup999( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1000( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1001( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1002( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1003( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1004( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1005( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1006( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1007( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1008( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1009( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1010( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1011( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1012( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1013( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1014( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1015( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1016( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1017( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1018( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1019( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1020( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1021( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1022( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1023( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1024( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1025( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1026( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1027( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1028( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1029( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1030( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1031( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1032( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1033( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1034( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1035( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1036( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1037( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1038( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1039( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1040( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1041( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1042( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1043( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1044( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1045( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1046( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1047( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1048( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1049( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1050( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1051( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1052( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1053( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1054( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1055( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1056( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1057( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1058( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1059( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1060( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1061( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1062( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1063( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1064( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1065( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1066( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1067( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1068( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1069( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1070( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1071( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1072( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1073( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1074( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1075( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1076( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1077( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1078( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1079( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1080( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1081( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1082( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1083( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1084( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1085( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1086( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1087( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1088( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1089( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1090( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1091( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1092( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1093( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1094( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1095( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1096( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1097( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1098( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1099( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1100( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1101( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1102( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1103( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1104( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1105( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1106( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1107( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1108( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1109( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1110( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1111( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1112( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1113( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1114( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1115( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1116( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1117( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1118( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1119( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1120( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1121( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1122( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1123( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1124( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1125( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1126( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1127( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1128( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1129( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1130( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1131( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1132( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1133( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1134( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1135( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1136( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1137( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1138( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1139( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1140( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1141( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1142( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1143( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1144( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1145( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1146( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1147( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1148( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1149( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1150( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1151( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1152( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1153( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1154( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1155( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1156( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1157( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1158( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1159( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1160( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1161( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1162( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1163( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1164( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1165( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1166( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1167( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1168( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1169( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1170( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1171( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1172( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1173( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1174( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1175( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1176( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1177( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1178( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1179( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1180( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1181( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1182( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1183( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1184( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1185( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1186( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1187( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1188( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1189( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1190( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1191( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1192( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1193( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1194( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1195( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1196( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1197( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1198( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1199( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1200( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1201( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1202( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1203( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1204( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1205( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1206( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1207( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1208( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1209( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1210( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1211( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1212( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1213( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1214( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1215( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1216( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1217( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1218( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1219( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1220( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1221( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1222( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1223( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1224( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1225( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1226( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1227( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1228( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1229( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1230( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1231( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1232( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1233( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1234( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1235( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1236( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1237( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1238( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1239( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1240( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1241( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1242( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1243( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1244( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1245( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1246( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1247( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1248( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1249( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1250( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1251( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1252( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1253( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1254( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1255( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1256( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1257( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1258( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1259( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1260( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1261( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1262( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1263( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1264( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1265( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1266( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1267( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1268( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1269( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1270( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1271( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1272( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1273( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1274( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1275( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1276( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1277( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1278( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1279( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1280( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1281( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1282( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1283( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1284( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1285( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1286( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1287( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1288( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1289( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1290( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1291( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1292( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1293( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1294( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1295( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1296( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1297( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1298( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1299( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1300( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1301( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1302( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1303( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1304( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1305( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1306( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1307( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1308( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1309( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1310( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1311( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1312( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1313( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1314( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1315( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1316( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1317( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1318( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1319( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1320( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1321( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1322( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1323( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1324( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1325( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1326( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1327( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1328( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1329( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1330( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1331( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1332( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1333( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1334( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1335( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1336( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1337( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1338( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1339( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1340( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1341( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1342( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1343( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1344( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1345( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1346( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1347( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1348( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1349( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1350( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1351( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1352( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1353( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1354( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1355( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1356( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1357( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1358( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1359( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1360( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1361( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1362( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1363( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1364( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1365( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1366( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1367( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1368( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1369( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1370( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1371( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1372( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1373( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1374( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1375( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1376( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1377( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1378( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1379( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1380( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1381( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1382( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1383( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1384( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1385( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1386( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1387( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1388( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1389( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1390( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1391( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1392( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1393( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1394( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1395( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1396( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1397( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1398( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1399( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1400( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1401( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1402( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1403( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1404( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1405( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1406( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1407( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1408( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1409( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1410( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1411( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1412( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1413( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1414( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1415( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1416( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1417( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1418( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1419( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1420( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1421( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1422( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1423( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1424( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1425( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1426( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1427( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1428( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1429( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1430( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1431( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1432( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1433( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1434( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1435( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1436( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1437( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1438( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1439( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1440( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1441( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1442( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1443( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1444( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1445( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1446( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1447( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1448( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1449( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1450( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1451( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1452( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1453( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1454( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1455( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1456( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1457( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1458( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1459( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1460( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1461( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1462( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1463( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1464( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1465( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1466( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1467( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1468( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1469( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1470( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1471( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1472( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1473( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1474( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1475( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1476( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1477( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1478( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1479( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1480( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1481( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1482( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1483( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1484( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1485( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1486( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1487( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1488( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1489( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1490( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1491( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1492( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1493( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1494( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1495( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1496( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1497( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1498( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1499( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1500( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1501( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1502( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1503( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1504( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1505( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1506( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1507( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1508( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1509( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1510( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1511( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1512( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1513( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1514( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1515( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1516( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1517( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1518( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1519( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1520( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1521( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1522( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1523( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1524( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1525( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1526( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1527( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1528( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1529( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1530( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1531( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1532( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1533( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1534( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1535( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1536( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1537( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1538( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1539( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1540( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1541( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1542( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1543( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1544( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1545( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1546( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1547( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1548( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1549( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
      diagramgroup1550( wfs, jamps, dcNGoodHel, couplings, channelIds, numerators, denominators, dcIPC, dcIPD );
#endif
#else
      // === C++ IMPLEMENTATION ===
      diagramgroup1( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD, cHelFlat, momenta, ihel );
      diagramgroup2( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup3( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup4( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup5( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup6( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup7( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup8( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup9( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup10( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup11( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup12( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup13( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup14( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup15( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup16( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup17( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup18( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup19( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup20( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup21( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup22( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup23( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup24( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup25( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup26( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup27( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup28( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup29( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup30( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup31( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup32( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup33( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup34( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup35( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup36( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup37( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup38( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup39( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup40( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup41( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup42( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup43( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup44( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup45( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup46( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup47( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup48( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup49( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup50( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup51( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup52( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup53( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup54( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup55( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup56( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup57( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup58( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup59( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup60( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup61( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup62( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup63( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup64( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup65( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup66( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup67( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup68( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup69( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup70( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup71( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup72( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup73( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup74( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup75( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup76( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup77( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup78( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup79( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup80( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup81( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup82( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup83( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup84( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup85( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup86( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup87( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup88( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup89( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup90( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup91( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup92( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup93( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup94( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup95( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup96( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup97( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup98( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup99( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup100( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup101( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup102( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup103( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup104( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup105( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup106( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup107( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup108( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup109( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup110( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup111( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup112( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup113( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup114( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup115( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup116( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup117( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup118( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup119( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup120( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup121( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup122( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup123( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup124( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup125( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup126( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup127( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup128( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup129( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup130( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup131( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup132( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup133( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup134( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup135( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup136( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup137( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup138( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup139( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup140( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup141( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup142( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup143( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup144( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup145( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup146( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup147( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup148( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup149( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup150( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup151( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup152( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup153( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup154( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup155( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup156( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup157( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup158( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup159( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup160( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup161( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup162( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup163( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup164( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup165( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup166( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup167( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup168( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup169( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup170( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup171( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup172( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup173( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup174( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup175( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup176( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup177( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup178( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup179( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup180( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup181( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup182( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup183( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup184( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup185( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup186( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup187( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup188( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup189( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup190( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup191( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup192( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup193( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup194( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup195( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup196( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup197( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup198( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup199( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup200( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup201( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup202( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup203( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup204( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup205( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup206( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup207( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup208( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup209( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup210( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup211( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup212( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup213( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup214( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup215( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup216( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup217( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup218( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup219( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup220( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup221( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup222( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup223( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup224( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup225( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup226( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup227( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup228( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup229( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup230( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup231( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup232( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup233( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup234( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup235( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup236( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup237( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup238( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup239( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup240( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup241( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup242( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup243( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup244( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup245( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup246( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup247( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup248( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup249( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup250( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup251( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup252( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup253( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup254( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup255( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup256( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup257( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup258( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup259( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup260( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup261( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup262( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup263( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup264( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup265( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup266( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup267( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup268( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup269( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup270( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup271( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup272( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup273( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup274( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup275( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup276( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup277( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup278( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup279( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup280( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup281( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup282( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup283( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup284( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup285( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup286( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup287( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup288( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup289( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup290( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup291( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup292( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup293( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup294( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup295( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup296( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup297( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup298( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup299( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup300( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup301( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup302( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup303( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup304( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup305( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup306( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup307( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup308( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup309( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup310( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup311( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup312( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup313( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup314( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup315( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup316( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup317( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup318( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup319( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup320( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup321( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup322( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup323( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup324( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup325( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup326( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup327( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup328( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup329( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup330( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup331( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup332( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup333( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup334( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup335( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup336( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup337( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup338( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup339( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup340( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup341( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup342( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup343( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup344( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup345( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup346( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup347( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup348( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup349( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup350( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup351( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup352( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup353( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup354( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup355( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup356( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup357( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup358( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup359( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup360( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup361( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup362( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup363( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup364( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup365( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup366( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup367( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup368( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup369( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup370( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup371( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup372( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup373( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup374( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup375( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup376( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup377( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup378( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup379( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup380( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup381( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup382( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup383( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup384( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup385( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup386( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup387( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup388( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup389( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup390( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup391( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup392( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup393( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup394( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup395( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup396( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup397( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup398( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup399( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup400( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup401( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup402( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup403( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup404( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup405( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup406( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup407( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup408( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup409( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup410( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup411( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup412( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup413( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup414( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup415( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup416( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup417( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup418( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup419( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup420( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup421( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup422( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup423( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup424( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup425( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup426( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup427( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup428( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup429( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup430( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup431( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup432( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup433( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup434( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup435( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup436( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup437( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup438( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup439( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup440( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup441( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup442( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup443( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup444( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup445( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup446( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup447( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup448( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup449( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup450( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup451( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup452( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup453( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup454( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup455( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup456( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup457( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup458( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup459( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup460( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup461( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup462( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup463( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup464( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup465( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup466( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup467( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup468( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup469( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup470( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup471( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup472( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup473( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup474( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup475( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup476( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup477( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup478( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup479( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup480( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup481( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup482( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup483( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup484( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup485( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup486( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup487( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup488( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup489( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup490( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup491( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup492( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup493( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup494( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup495( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup496( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup497( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup498( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup499( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup500( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup501( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup502( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup503( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup504( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup505( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup506( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup507( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup508( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup509( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup510( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup511( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup512( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup513( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup514( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup515( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup516( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup517( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup518( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup519( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup520( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup521( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup522( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup523( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup524( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup525( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup526( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup527( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup528( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup529( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup530( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup531( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup532( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup533( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup534( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup535( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup536( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup537( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup538( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup539( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup540( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup541( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup542( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup543( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup544( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup545( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup546( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup547( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup548( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup549( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup550( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup551( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup552( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup553( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup554( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup555( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup556( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup557( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup558( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup559( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup560( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup561( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup562( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup563( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup564( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup565( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup566( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup567( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup568( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup569( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup570( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup571( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup572( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup573( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup574( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup575( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup576( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup577( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup578( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup579( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup580( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup581( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup582( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup583( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup584( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup585( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup586( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup587( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup588( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup589( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup590( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup591( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup592( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup593( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup594( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup595( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup596( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup597( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup598( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup599( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup600( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup601( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup602( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup603( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup604( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup605( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup606( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup607( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup608( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup609( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup610( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup611( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup612( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup613( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup614( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup615( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup616( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup617( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup618( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup619( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup620( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup621( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup622( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup623( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup624( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup625( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup626( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup627( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup628( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup629( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup630( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup631( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup632( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup633( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup634( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup635( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup636( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup637( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup638( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup639( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup640( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup641( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup642( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup643( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup644( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup645( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup646( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup647( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup648( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup649( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup650( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup651( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup652( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup653( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup654( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup655( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup656( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup657( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup658( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup659( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup660( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup661( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup662( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup663( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup664( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup665( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup666( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup667( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup668( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup669( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup670( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup671( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup672( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup673( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup674( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup675( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup676( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup677( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup678( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup679( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup680( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup681( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup682( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup683( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup684( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup685( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup686( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup687( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup688( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup689( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup690( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup691( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup692( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup693( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup694( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup695( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup696( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup697( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup698( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup699( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup700( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup701( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup702( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup703( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup704( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup705( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup706( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup707( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup708( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup709( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup710( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup711( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup712( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup713( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup714( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup715( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup716( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup717( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup718( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup719( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup720( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup721( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup722( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup723( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup724( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup725( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup726( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup727( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup728( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup729( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup730( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup731( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup732( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup733( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup734( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup735( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup736( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup737( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup738( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup739( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup740( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup741( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup742( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup743( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup744( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup745( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup746( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup747( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup748( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup749( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup750( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup751( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup752( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup753( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup754( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup755( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup756( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup757( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup758( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup759( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup760( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup761( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup762( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup763( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup764( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup765( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup766( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup767( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup768( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup769( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup770( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup771( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup772( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup773( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup774( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup775( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup776( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup777( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup778( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup779( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup780( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup781( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup782( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup783( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup784( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup785( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup786( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup787( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup788( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup789( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup790( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup791( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup792( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup793( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup794( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup795( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup796( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup797( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup798( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup799( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup800( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup801( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup802( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup803( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup804( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup805( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup806( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup807( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup808( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup809( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup810( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup811( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup812( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup813( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup814( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup815( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup816( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup817( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup818( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup819( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup820( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup821( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup822( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup823( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup824( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup825( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup826( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup827( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup828( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup829( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup830( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup831( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup832( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup833( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup834( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup835( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup836( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup837( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup838( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup839( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup840( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup841( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup842( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup843( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup844( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup845( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup846( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup847( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup848( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup849( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup850( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup851( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup852( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup853( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup854( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup855( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup856( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup857( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup858( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup859( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup860( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup861( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup862( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup863( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup864( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup865( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup866( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup867( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup868( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup869( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup870( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup871( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup872( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup873( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup874( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup875( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup876( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup877( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup878( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup879( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup880( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup881( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup882( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup883( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup884( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup885( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup886( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup887( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup888( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup889( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup890( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup891( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup892( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup893( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup894( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup895( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup896( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup897( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup898( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup899( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup900( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup901( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup902( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup903( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup904( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup905( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup906( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup907( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup908( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup909( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup910( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup911( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup912( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup913( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup914( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup915( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup916( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup917( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup918( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup919( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup920( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup921( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup922( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup923( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup924( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup925( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup926( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup927( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup928( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup929( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup930( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup931( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup932( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup933( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup934( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup935( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup936( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup937( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup938( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup939( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup940( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup941( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup942( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup943( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup944( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup945( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup946( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup947( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup948( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup949( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup950( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup951( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup952( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup953( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup954( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup955( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup956( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup957( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup958( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup959( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup960( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup961( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup962( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup963( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup964( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup965( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup966( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup967( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup968( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup969( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup970( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup971( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup972( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup973( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup974( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup975( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup976( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup977( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup978( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup979( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup980( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup981( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup982( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup983( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup984( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup985( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup986( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup987( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup988( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup989( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup990( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup991( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup992( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup993( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup994( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup995( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup996( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup997( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup998( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup999( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1000( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1001( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1002( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1003( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1004( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1005( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1006( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1007( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1008( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1009( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1010( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1011( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1012( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1013( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1014( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1015( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1016( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1017( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1018( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1019( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1020( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1021( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1022( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1023( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1024( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1025( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1026( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1027( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1028( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1029( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1030( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1031( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1032( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1033( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1034( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1035( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1036( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1037( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1038( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1039( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1040( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1041( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1042( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1043( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1044( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1045( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1046( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1047( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1048( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1049( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1050( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1051( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1052( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1053( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1054( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1055( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1056( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1057( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1058( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1059( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1060( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1061( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1062( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1063( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1064( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1065( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1066( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1067( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1068( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1069( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1070( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1071( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1072( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1073( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1074( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1075( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1076( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1077( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1078( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1079( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1080( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1081( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1082( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1083( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1084( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1085( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1086( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1087( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1088( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1089( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1090( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1091( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1092( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1093( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1094( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1095( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1096( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1097( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1098( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1099( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1100( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1101( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1102( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1103( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1104( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1105( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1106( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1107( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1108( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1109( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1110( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1111( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1112( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1113( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1114( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1115( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1116( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1117( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1118( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1119( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1120( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1121( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1122( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1123( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1124( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1125( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1126( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1127( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1128( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1129( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1130( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1131( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1132( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1133( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1134( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1135( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1136( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1137( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1138( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1139( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1140( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1141( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1142( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1143( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1144( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1145( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1146( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1147( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1148( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1149( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1150( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1151( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1152( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1153( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1154( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1155( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1156( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1157( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1158( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1159( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1160( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1161( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1162( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1163( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1164( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1165( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1166( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1167( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1168( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1169( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1170( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1171( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1172( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1173( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1174( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1175( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1176( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1177( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1178( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1179( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1180( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1181( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1182( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1183( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1184( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1185( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1186( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1187( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1188( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1189( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1190( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1191( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1192( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1193( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1194( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1195( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1196( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1197( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1198( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1199( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1200( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1201( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1202( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1203( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1204( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1205( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1206( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1207( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1208( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1209( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1210( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1211( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1212( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1213( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1214( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1215( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1216( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1217( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1218( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1219( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1220( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1221( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1222( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1223( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1224( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1225( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1226( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1227( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1228( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1229( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1230( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1231( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1232( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1233( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1234( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1235( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1236( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1237( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1238( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1239( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1240( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1241( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1242( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1243( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1244( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1245( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1246( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1247( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1248( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1249( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1250( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1251( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1252( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1253( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1254( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1255( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1256( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1257( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1258( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1259( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1260( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1261( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1262( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1263( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1264( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1265( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1266( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1267( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1268( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1269( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1270( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1271( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1272( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1273( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1274( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1275( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1276( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1277( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1278( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1279( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1280( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1281( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1282( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1283( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1284( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1285( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1286( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1287( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1288( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1289( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1290( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1291( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1292( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1293( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1294( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1295( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1296( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1297( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1298( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1299( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1300( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1301( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1302( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1303( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1304( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1305( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1306( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1307( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1308( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1309( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1310( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1311( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1312( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1313( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1314( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1315( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1316( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1317( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1318( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1319( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1320( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1321( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1322( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1323( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1324( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1325( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1326( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1327( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1328( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1329( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1330( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1331( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1332( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1333( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1334( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1335( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1336( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1337( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1338( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1339( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1340( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1341( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1342( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1343( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1344( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1345( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1346( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1347( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1348( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1349( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1350( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1351( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1352( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1353( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1354( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1355( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1356( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1357( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1358( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1359( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1360( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1361( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1362( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1363( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1364( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1365( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1366( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1367( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1368( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1369( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1370( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1371( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1372( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1373( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1374( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1375( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1376( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1377( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1378( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1379( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1380( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1381( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1382( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1383( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1384( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1385( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1386( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1387( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1388( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1389( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1390( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1391( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1392( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1393( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1394( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1395( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1396( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1397( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1398( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1399( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1400( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1401( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1402( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1403( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1404( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1405( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1406( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1407( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1408( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1409( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1410( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1411( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1412( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1413( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1414( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1415( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1416( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1417( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1418( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1419( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1420( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1421( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1422( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1423( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1424( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1425( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1426( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1427( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1428( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1429( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1430( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1431( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1432( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1433( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1434( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1435( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1436( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1437( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1438( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1439( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1440( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1441( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1442( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1443( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1444( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1445( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1446( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1447( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1448( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1449( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1450( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1451( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1452( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1453( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1454( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1455( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1456( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1457( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1458( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1459( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1460( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1461( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1462( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1463( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1464( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1465( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1466( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1467( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1468( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1469( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1470( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1471( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1472( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1473( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1474( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1475( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1476( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1477( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1478( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1479( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1480( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1481( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1482( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1483( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1484( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1485( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1486( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1487( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1488( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1489( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1490( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1491( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1492( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1493( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1494( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1495( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1496( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1497( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1498( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1499( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1500( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1501( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1502( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1503( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1504( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1505( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1506( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1507( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1508( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1509( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1510( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1511( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1512( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1513( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1514( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1515( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1516( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1517( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1518( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1519( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1520( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1521( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1522( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1523( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1524( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1525( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1526( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1527( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1528( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1529( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1530( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1531( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1532( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1533( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1534( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1535( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1536( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1537( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1538( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1539( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1540( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1541( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1542( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1543( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1544( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1545( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1546( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1547( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1548( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1549( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
      diagramgroup1550( wfs, jamp_sv, COUPs, channelIds, numerators, denominators, cIPC, cIPD );
#endif
    }
    // *****************************
    // *** END LOOP ON IPARITY ***
    // *****************************

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
      { -1, -1, -1, 1, -1, -1, -1, -1 },
      { -1, -1, -1, 1, -1, -1, -1, 1 },
      { -1, -1, -1, 1, -1, -1, 1, -1 },
      { -1, -1, -1, 1, -1, -1, 1, 1 },
      { -1, -1, -1, 1, -1, 1, -1, -1 },
      { -1, -1, -1, 1, -1, 1, -1, 1 },
      { -1, -1, -1, 1, -1, 1, 1, -1 },
      { -1, -1, -1, 1, -1, 1, 1, 1 },
      { -1, -1, -1, 1, 1, -1, -1, -1 },
      { -1, -1, -1, 1, 1, -1, -1, 1 },
      { -1, -1, -1, 1, 1, -1, 1, -1 },
      { -1, -1, -1, 1, 1, -1, 1, 1 },
      { -1, -1, -1, 1, 1, 1, -1, -1 },
      { -1, -1, -1, 1, 1, 1, -1, 1 },
      { -1, -1, -1, 1, 1, 1, 1, -1 },
      { -1, -1, -1, 1, 1, 1, 1, 1 },
      { -1, -1, -1, -1, -1, -1, -1, -1 },
      { -1, -1, -1, -1, -1, -1, -1, 1 },
      { -1, -1, -1, -1, -1, -1, 1, -1 },
      { -1, -1, -1, -1, -1, -1, 1, 1 },
      { -1, -1, -1, -1, -1, 1, -1, -1 },
      { -1, -1, -1, -1, -1, 1, -1, 1 },
      { -1, -1, -1, -1, -1, 1, 1, -1 },
      { -1, -1, -1, -1, -1, 1, 1, 1 },
      { -1, -1, -1, -1, 1, -1, -1, -1 },
      { -1, -1, -1, -1, 1, -1, -1, 1 },
      { -1, -1, -1, -1, 1, -1, 1, -1 },
      { -1, -1, -1, -1, 1, -1, 1, 1 },
      { -1, -1, -1, -1, 1, 1, -1, -1 },
      { -1, -1, -1, -1, 1, 1, -1, 1 },
      { -1, -1, -1, -1, 1, 1, 1, -1 },
      { -1, -1, -1, -1, 1, 1, 1, 1 },
      { -1, -1, 1, 1, -1, -1, -1, -1 },
      { -1, -1, 1, 1, -1, -1, -1, 1 },
      { -1, -1, 1, 1, -1, -1, 1, -1 },
      { -1, -1, 1, 1, -1, -1, 1, 1 },
      { -1, -1, 1, 1, -1, 1, -1, -1 },
      { -1, -1, 1, 1, -1, 1, -1, 1 },
      { -1, -1, 1, 1, -1, 1, 1, -1 },
      { -1, -1, 1, 1, -1, 1, 1, 1 },
      { -1, -1, 1, 1, 1, -1, -1, -1 },
      { -1, -1, 1, 1, 1, -1, -1, 1 },
      { -1, -1, 1, 1, 1, -1, 1, -1 },
      { -1, -1, 1, 1, 1, -1, 1, 1 },
      { -1, -1, 1, 1, 1, 1, -1, -1 },
      { -1, -1, 1, 1, 1, 1, -1, 1 },
      { -1, -1, 1, 1, 1, 1, 1, -1 },
      { -1, -1, 1, 1, 1, 1, 1, 1 },
      { -1, -1, 1, -1, -1, -1, -1, -1 },
      { -1, -1, 1, -1, -1, -1, -1, 1 },
      { -1, -1, 1, -1, -1, -1, 1, -1 },
      { -1, -1, 1, -1, -1, -1, 1, 1 },
      { -1, -1, 1, -1, -1, 1, -1, -1 },
      { -1, -1, 1, -1, -1, 1, -1, 1 },
      { -1, -1, 1, -1, -1, 1, 1, -1 },
      { -1, -1, 1, -1, -1, 1, 1, 1 },
      { -1, -1, 1, -1, 1, -1, -1, -1 },
      { -1, -1, 1, -1, 1, -1, -1, 1 },
      { -1, -1, 1, -1, 1, -1, 1, -1 },
      { -1, -1, 1, -1, 1, -1, 1, 1 },
      { -1, -1, 1, -1, 1, 1, -1, -1 },
      { -1, -1, 1, -1, 1, 1, -1, 1 },
      { -1, -1, 1, -1, 1, 1, 1, -1 },
      { -1, -1, 1, -1, 1, 1, 1, 1 },
      { -1, 1, -1, 1, -1, -1, -1, -1 },
      { -1, 1, -1, 1, -1, -1, -1, 1 },
      { -1, 1, -1, 1, -1, -1, 1, -1 },
      { -1, 1, -1, 1, -1, -1, 1, 1 },
      { -1, 1, -1, 1, -1, 1, -1, -1 },
      { -1, 1, -1, 1, -1, 1, -1, 1 },
      { -1, 1, -1, 1, -1, 1, 1, -1 },
      { -1, 1, -1, 1, -1, 1, 1, 1 },
      { -1, 1, -1, 1, 1, -1, -1, -1 },
      { -1, 1, -1, 1, 1, -1, -1, 1 },
      { -1, 1, -1, 1, 1, -1, 1, -1 },
      { -1, 1, -1, 1, 1, -1, 1, 1 },
      { -1, 1, -1, 1, 1, 1, -1, -1 },
      { -1, 1, -1, 1, 1, 1, -1, 1 },
      { -1, 1, -1, 1, 1, 1, 1, -1 },
      { -1, 1, -1, 1, 1, 1, 1, 1 },
      { -1, 1, -1, -1, -1, -1, -1, -1 },
      { -1, 1, -1, -1, -1, -1, -1, 1 },
      { -1, 1, -1, -1, -1, -1, 1, -1 },
      { -1, 1, -1, -1, -1, -1, 1, 1 },
      { -1, 1, -1, -1, -1, 1, -1, -1 },
      { -1, 1, -1, -1, -1, 1, -1, 1 },
      { -1, 1, -1, -1, -1, 1, 1, -1 },
      { -1, 1, -1, -1, -1, 1, 1, 1 },
      { -1, 1, -1, -1, 1, -1, -1, -1 },
      { -1, 1, -1, -1, 1, -1, -1, 1 },
      { -1, 1, -1, -1, 1, -1, 1, -1 },
      { -1, 1, -1, -1, 1, -1, 1, 1 },
      { -1, 1, -1, -1, 1, 1, -1, -1 },
      { -1, 1, -1, -1, 1, 1, -1, 1 },
      { -1, 1, -1, -1, 1, 1, 1, -1 },
      { -1, 1, -1, -1, 1, 1, 1, 1 },
      { -1, 1, 1, 1, -1, -1, -1, -1 },
      { -1, 1, 1, 1, -1, -1, -1, 1 },
      { -1, 1, 1, 1, -1, -1, 1, -1 },
      { -1, 1, 1, 1, -1, -1, 1, 1 },
      { -1, 1, 1, 1, -1, 1, -1, -1 },
      { -1, 1, 1, 1, -1, 1, -1, 1 },
      { -1, 1, 1, 1, -1, 1, 1, -1 },
      { -1, 1, 1, 1, -1, 1, 1, 1 },
      { -1, 1, 1, 1, 1, -1, -1, -1 },
      { -1, 1, 1, 1, 1, -1, -1, 1 },
      { -1, 1, 1, 1, 1, -1, 1, -1 },
      { -1, 1, 1, 1, 1, -1, 1, 1 },
      { -1, 1, 1, 1, 1, 1, -1, -1 },
      { -1, 1, 1, 1, 1, 1, -1, 1 },
      { -1, 1, 1, 1, 1, 1, 1, -1 },
      { -1, 1, 1, 1, 1, 1, 1, 1 },
      { -1, 1, 1, -1, -1, -1, -1, -1 },
      { -1, 1, 1, -1, -1, -1, -1, 1 },
      { -1, 1, 1, -1, -1, -1, 1, -1 },
      { -1, 1, 1, -1, -1, -1, 1, 1 },
      { -1, 1, 1, -1, -1, 1, -1, -1 },
      { -1, 1, 1, -1, -1, 1, -1, 1 },
      { -1, 1, 1, -1, -1, 1, 1, -1 },
      { -1, 1, 1, -1, -1, 1, 1, 1 },
      { -1, 1, 1, -1, 1, -1, -1, -1 },
      { -1, 1, 1, -1, 1, -1, -1, 1 },
      { -1, 1, 1, -1, 1, -1, 1, -1 },
      { -1, 1, 1, -1, 1, -1, 1, 1 },
      { -1, 1, 1, -1, 1, 1, -1, -1 },
      { -1, 1, 1, -1, 1, 1, -1, 1 },
      { -1, 1, 1, -1, 1, 1, 1, -1 },
      { -1, 1, 1, -1, 1, 1, 1, 1 },
      { 1, -1, -1, 1, -1, -1, -1, -1 },
      { 1, -1, -1, 1, -1, -1, -1, 1 },
      { 1, -1, -1, 1, -1, -1, 1, -1 },
      { 1, -1, -1, 1, -1, -1, 1, 1 },
      { 1, -1, -1, 1, -1, 1, -1, -1 },
      { 1, -1, -1, 1, -1, 1, -1, 1 },
      { 1, -1, -1, 1, -1, 1, 1, -1 },
      { 1, -1, -1, 1, -1, 1, 1, 1 },
      { 1, -1, -1, 1, 1, -1, -1, -1 },
      { 1, -1, -1, 1, 1, -1, -1, 1 },
      { 1, -1, -1, 1, 1, -1, 1, -1 },
      { 1, -1, -1, 1, 1, -1, 1, 1 },
      { 1, -1, -1, 1, 1, 1, -1, -1 },
      { 1, -1, -1, 1, 1, 1, -1, 1 },
      { 1, -1, -1, 1, 1, 1, 1, -1 },
      { 1, -1, -1, 1, 1, 1, 1, 1 },
      { 1, -1, -1, -1, -1, -1, -1, -1 },
      { 1, -1, -1, -1, -1, -1, -1, 1 },
      { 1, -1, -1, -1, -1, -1, 1, -1 },
      { 1, -1, -1, -1, -1, -1, 1, 1 },
      { 1, -1, -1, -1, -1, 1, -1, -1 },
      { 1, -1, -1, -1, -1, 1, -1, 1 },
      { 1, -1, -1, -1, -1, 1, 1, -1 },
      { 1, -1, -1, -1, -1, 1, 1, 1 },
      { 1, -1, -1, -1, 1, -1, -1, -1 },
      { 1, -1, -1, -1, 1, -1, -1, 1 },
      { 1, -1, -1, -1, 1, -1, 1, -1 },
      { 1, -1, -1, -1, 1, -1, 1, 1 },
      { 1, -1, -1, -1, 1, 1, -1, -1 },
      { 1, -1, -1, -1, 1, 1, -1, 1 },
      { 1, -1, -1, -1, 1, 1, 1, -1 },
      { 1, -1, -1, -1, 1, 1, 1, 1 },
      { 1, -1, 1, 1, -1, -1, -1, -1 },
      { 1, -1, 1, 1, -1, -1, -1, 1 },
      { 1, -1, 1, 1, -1, -1, 1, -1 },
      { 1, -1, 1, 1, -1, -1, 1, 1 },
      { 1, -1, 1, 1, -1, 1, -1, -1 },
      { 1, -1, 1, 1, -1, 1, -1, 1 },
      { 1, -1, 1, 1, -1, 1, 1, -1 },
      { 1, -1, 1, 1, -1, 1, 1, 1 },
      { 1, -1, 1, 1, 1, -1, -1, -1 },
      { 1, -1, 1, 1, 1, -1, -1, 1 },
      { 1, -1, 1, 1, 1, -1, 1, -1 },
      { 1, -1, 1, 1, 1, -1, 1, 1 },
      { 1, -1, 1, 1, 1, 1, -1, -1 },
      { 1, -1, 1, 1, 1, 1, -1, 1 },
      { 1, -1, 1, 1, 1, 1, 1, -1 },
      { 1, -1, 1, 1, 1, 1, 1, 1 },
      { 1, -1, 1, -1, -1, -1, -1, -1 },
      { 1, -1, 1, -1, -1, -1, -1, 1 },
      { 1, -1, 1, -1, -1, -1, 1, -1 },
      { 1, -1, 1, -1, -1, -1, 1, 1 },
      { 1, -1, 1, -1, -1, 1, -1, -1 },
      { 1, -1, 1, -1, -1, 1, -1, 1 },
      { 1, -1, 1, -1, -1, 1, 1, -1 },
      { 1, -1, 1, -1, -1, 1, 1, 1 },
      { 1, -1, 1, -1, 1, -1, -1, -1 },
      { 1, -1, 1, -1, 1, -1, -1, 1 },
      { 1, -1, 1, -1, 1, -1, 1, -1 },
      { 1, -1, 1, -1, 1, -1, 1, 1 },
      { 1, -1, 1, -1, 1, 1, -1, -1 },
      { 1, -1, 1, -1, 1, 1, -1, 1 },
      { 1, -1, 1, -1, 1, 1, 1, -1 },
      { 1, -1, 1, -1, 1, 1, 1, 1 },
      { 1, 1, -1, 1, -1, -1, -1, -1 },
      { 1, 1, -1, 1, -1, -1, -1, 1 },
      { 1, 1, -1, 1, -1, -1, 1, -1 },
      { 1, 1, -1, 1, -1, -1, 1, 1 },
      { 1, 1, -1, 1, -1, 1, -1, -1 },
      { 1, 1, -1, 1, -1, 1, -1, 1 },
      { 1, 1, -1, 1, -1, 1, 1, -1 },
      { 1, 1, -1, 1, -1, 1, 1, 1 },
      { 1, 1, -1, 1, 1, -1, -1, -1 },
      { 1, 1, -1, 1, 1, -1, -1, 1 },
      { 1, 1, -1, 1, 1, -1, 1, -1 },
      { 1, 1, -1, 1, 1, -1, 1, 1 },
      { 1, 1, -1, 1, 1, 1, -1, -1 },
      { 1, 1, -1, 1, 1, 1, -1, 1 },
      { 1, 1, -1, 1, 1, 1, 1, -1 },
      { 1, 1, -1, 1, 1, 1, 1, 1 },
      { 1, 1, -1, -1, -1, -1, -1, -1 },
      { 1, 1, -1, -1, -1, -1, -1, 1 },
      { 1, 1, -1, -1, -1, -1, 1, -1 },
      { 1, 1, -1, -1, -1, -1, 1, 1 },
      { 1, 1, -1, -1, -1, 1, -1, -1 },
      { 1, 1, -1, -1, -1, 1, -1, 1 },
      { 1, 1, -1, -1, -1, 1, 1, -1 },
      { 1, 1, -1, -1, -1, 1, 1, 1 },
      { 1, 1, -1, -1, 1, -1, -1, -1 },
      { 1, 1, -1, -1, 1, -1, -1, 1 },
      { 1, 1, -1, -1, 1, -1, 1, -1 },
      { 1, 1, -1, -1, 1, -1, 1, 1 },
      { 1, 1, -1, -1, 1, 1, -1, -1 },
      { 1, 1, -1, -1, 1, 1, -1, 1 },
      { 1, 1, -1, -1, 1, 1, 1, -1 },
      { 1, 1, -1, -1, 1, 1, 1, 1 },
      { 1, 1, 1, 1, -1, -1, -1, -1 },
      { 1, 1, 1, 1, -1, -1, -1, 1 },
      { 1, 1, 1, 1, -1, -1, 1, -1 },
      { 1, 1, 1, 1, -1, -1, 1, 1 },
      { 1, 1, 1, 1, -1, 1, -1, -1 },
      { 1, 1, 1, 1, -1, 1, -1, 1 },
      { 1, 1, 1, 1, -1, 1, 1, -1 },
      { 1, 1, 1, 1, -1, 1, 1, 1 },
      { 1, 1, 1, 1, 1, -1, -1, -1 },
      { 1, 1, 1, 1, 1, -1, -1, 1 },
      { 1, 1, 1, 1, 1, -1, 1, -1 },
      { 1, 1, 1, 1, 1, -1, 1, 1 },
      { 1, 1, 1, 1, 1, 1, -1, -1 },
      { 1, 1, 1, 1, 1, 1, -1, 1 },
      { 1, 1, 1, 1, 1, 1, 1, -1 },
      { 1, 1, 1, 1, 1, 1, 1, 1 },
      { 1, 1, 1, -1, -1, -1, -1, -1 },
      { 1, 1, 1, -1, -1, -1, -1, 1 },
      { 1, 1, 1, -1, -1, -1, 1, -1 },
      { 1, 1, 1, -1, -1, -1, 1, 1 },
      { 1, 1, 1, -1, -1, 1, -1, -1 },
      { 1, 1, 1, -1, -1, 1, -1, 1 },
      { 1, 1, 1, -1, -1, 1, 1, -1 },
      { 1, 1, 1, -1, -1, 1, 1, 1 },
      { 1, 1, 1, -1, 1, -1, -1, -1 },
      { 1, 1, 1, -1, 1, -1, -1, 1 },
      { 1, 1, 1, -1, 1, -1, 1, -1 },
      { 1, 1, 1, -1, 1, -1, 1, 1 },
      { 1, 1, 1, -1, 1, 1, -1, -1 },
      { 1, 1, 1, -1, 1, 1, -1, 1 },
      { 1, 1, 1, -1, 1, 1, 1, -1 },
      { 1, 1, 1, -1, 1, 1, 1, 1 } };
#ifdef MGONGPUCPP_GPUIMPL
    gpuMemcpyToSymbol( dcHel, tHel, ncomb * npar * sizeof( short ) );
#ifndef MGONGPU_RDC_DIAGRAMS
    gpuGetSymbolAddress( (void**)( &cHelFlat ), dcHel );
#endif
#else
    memcpy( cHel, tHel, ncomb * npar * sizeof( short ) );
#endif

    // Enable SIGFPE traps for Floating Point Exceptions
#ifdef MGONGPUCPP_DEBUG
    fpeEnable();
#endif
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
    m_masses.push_back( m_pars->ZERO );
#ifdef MGONGPUCPP_GPUIMPL
    // Create the normalized color matrix in device memory
    createNormalizedColorMatrix();
#endif
    // Read physics parameters like masses and couplings from user configuration files (static: initialize once)
    // Then copy them to CUDA constant memory (issue #39) or its C++ emulation in file-scope static memory
    const fptype tIPD[nIPD] = { (fptype)m_pars->mdl_MT, (fptype)m_pars->mdl_WT };
    //const cxtype tIPC[0] = { ... }; // nIPC=0
#ifdef MGONGPUCPP_GPUIMPL
    gpuMemcpyToSymbol( dcIPD, tIPD, nIPD * sizeof( fptype ) );
    //gpuMemcpyToSymbol( dcIPC, tIPC, 0 * sizeof( cxtype ) ); // nIPC=0
    if constexpr( nIPD > 0 ) gpuGetSymbolAddress( (void**)( &cIPD ), dcIPD );
    if constexpr( nIPC > 0 ) gpuGetSymbolAddress( (void**)( &cIPC ), dcIPC );
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
    m_masses.push_back( Parameters_sm::ZERO );
#ifdef MGONGPUCPP_GPUIMPL
#ifdef __HIPCC__
#warning HRDCOD=1 in CUDACPP is no longer supported on HIP
#warning This code builds but fails at runtime "Cannot create GlobalVar Obj for symbol: _ZN9mg5amcGpuL5dcIPDE"
#endif
    if constexpr( nIPD > 0 ) gpuGetSymbolAddress( (void**)( &cIPD ), dcIPD );
    if constexpr( nIPC > 0 ) gpuGetSymbolAddress( (void**)( &cIPC ), dcIPC );
    // Create the normalized color matrix in device memory
    createNormalizedColorMatrix();
#endif
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
                       fptype* allJamps,           // tmp: jamp[ncolor*2*nevt] _for one helicity_ (reused in the getGoodHel helicity loop)
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
      constexpr int nOneHel = 1; // use a jamp buffer for a single helicity
      gpuMemcpyToSymbol( dcNGoodHel, &nOneHel, sizeof( int ) );
      cNGoodHel = nOneHel; // fix nasty bug (which was causing failures only in heftggbb)
      // NEW IMPLEMENTATION OF GETGOODHEL (#630): RESET THE RUNNING SUM OVER HELICITIES TO 0 BEFORE ADDING A NEW HELICITY
      gpuMemset( allMEs, 0, maxtry * sizeof( fptype ) );
      gpuMemset( allJamps, 0, maxtry * ncolor * mgOnGpu::nx2 * sizeof( fptype ) );
      // NB: color_sum ADDS |M|^2 for one helicity to the running sum of |M|^2 over helicities for the given event(s)
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      constexpr unsigned int* allChannelIds = nullptr; // disable multichannel single-diagram enhancement
#ifndef MGONGPU_RDC_DIAGRAMS
      calculate_jamps( allmomenta, allcouplings, allJamps, allWfs, allChannelIds, allNumerators, allDenominators, 0, gpublocks, gputhreads, ihel );
#else
      gpuLaunchKernelStream( calculate_jamps, gpublocks, gputhreads, 0, allmomenta, allcouplings, allJamps, allWfs, allChannelIds, allNumerators, allDenominators, ihel );
#endif
#else
#ifndef MGONGPU_RDC_DIAGRAMS
      calculate_jamps( allmomenta, allcouplings, allJamps, allWfs, 0, gpublocks, gputhreads, ihel );
#else
      gpuLaunchKernelStream( calculate_jamps, gpublocks, gputhreads, 0, allmomenta, allcouplings, allJamps, allWfs, ihel );
#endif
#endif
      gpuLaunchKernel( color_sum_kernel, gpublocks, gputhreads, allMEs, allJamps, nOneHel );
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
        cxtype_sv jamp_sv_1or2[2 * ncolor] = {}; // all zeros
#else
        cxtype_sv jamp_sv_1or2[ncolor] = {}; // all zeros
#endif
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL /* clang-format off */
        constexpr unsigned int channelId = 0; // disable multichannel single-diagram enhancement
        calculate_jamps( allmomenta, allcouplings, jamp_sv_1or2, channelId, allNumerators, allDenominators, ievt00, ihel ); //maxtry?
#else
        calculate_jamps( allmomenta, allcouplings, jamp_sv_1or2, ievt00, ihel ); //maxtry?
#endif /* clang-format on */
        color_sum_cpu( allMEs, jamp_sv_1or2, ievt00 );
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
                    const unsigned int* allChannelIds, // input: multichannel channelIds[nevt] (1 to #diagrams); nullptr to disable SDE (#899/#911)
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
                 fptype* colAllJamp2s,      // output: allJamp2s[ncolor][nevt] super-buffer, sum over col/hel (nullptr to disable)
                 const int nGoodHel )       // input: number of good helicities
  {
    using J_ACCESS = DeviceAccessJamp;
    using J2_ACCESS = DeviceAccessJamp2;
    constexpr int ihel0 = 0; // the allJamps buffer already points to a specific helicity _within a super-buffer for dcNGoodHel helicities_
    for( int icol = 0; icol < ncolor; icol++ )
      // NB: atomicAdd is needed after moving to cuda streams with one helicity per stream!
      atomicAdd( &J2_ACCESS::kernelAccessIcol( colAllJamp2s, icol ),
                 cxabs2( J_ACCESS::kernelAccessIcolIhelNhelConst( allJamps, icol, ihel0, nGoodHel ) ) );
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
            const unsigned int* allChannelIds,  // input: multichannel channelIds[nevt] (1 to #diagrams); nullptr to disable SDE (#899/#911)
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
            fptype2* ghelAllBlasTmp,            // tmp: allBlasTmp super-buffer for nGoodHel <= ncomb individual helicities (index is ighel)
            gpuBlasHandle_t* pBlasHandle,       // input: cuBLAS/hipBLAS handle
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
            const unsigned int* allChannelIds,  // input: multichannel channelIds[nevt] (1 to #diagrams); nullptr to disable SDE (#899/#911)
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
    constexpr int helcolDenominators[1] = { 6144 }; // assume nprocesses == 1 (#272 and #343)

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
    // *** PART 0a - CUDA ***
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
    // *** PART 0b - C++ ***
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
      fptype* hAllJamps = ghelAllJamps + ighel * nevt; // HACK: bypass DeviceAccessJamp (consistent with layout defined there)
      fptype* hAllWfs = ( ghelAllWfs ? ghelAllWfs + ighel * nwf * nevt * nw6 * mgOnGpu::nx2 : nullptr );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      fptype* hAllNumerators = ghelAllNumerators + ighel * nevt;
      fptype* hAllDenominators = ghelAllDenominators + ighel * nevt;
#ifndef MGONGPU_RDC_DIAGRAMS
      calculate_jamps( allmomenta, allcouplings, hAllJamps, hAllWfs, allChannelIds, hAllNumerators, hAllDenominators, ghelStreams[ighel], gpublocks, gputhreads, ihel );
#else
      gpuLaunchKernelStream( calculate_jamps, gpublocks, gputhreads, ghelStreams[ighel], allmomenta, allcouplings, hAllJamps, hAllWfs, allChannelIds, hAllNumerators, hAllDenominators, ihel );
#endif
#else
#ifndef MGONGPU_RDC_DIAGRAMS
      calculate_jamps( allmomenta, allcouplings, hAllJamps, hAllWfs, ghelStreams[ighel], gpublocks, gputhreads, ihel );
#else
      gpuLaunchKernelStream( calculate_jamps, gpublocks, gputhreads, ghelStreams[ighel], allmomenta, allcouplings, hAllJamps, hAllWfs, ihel );
#endif
#endif
    }
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // (1b) Then, in multichannel mode, also compute the running sums over helicities of squared jamp2s within each helicity stream
    for( int ighel = 0; ighel < cNGoodHel; ighel++ )
    {
      fptype* hAllJamps = ghelAllJamps + ighel * nevt; // HACK: bypass DeviceAccessJamp (consistent with layout defined there)
      gpuLaunchKernelStream( update_jamp2s, gpublocks, gputhreads, ghelStreams[ighel], hAllJamps, colAllJamp2s, cNGoodHel );
    }
#endif
    // (2) Then compute the ME for that helicity from the color sum of QCD partial amplitudes jamps
    color_sum_gpu( ghelAllMEs, ghelAllJamps, ghelAllBlasTmp, pBlasHandle, ghelStreams, cNGoodHel, gpublocks, gputhreads );
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
#else /* clang-format off */
    const int npagV2 = npagV; // loop on one SIMD page (neppV events) at a time
#endif /* clang-format on */
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
#else /* clang-format off */
      const int ievt00 = ipagV2 * neppV; // loop on one SIMD page (neppV events) at a time
#endif /* clang-format on */
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
        cxtype_sv jamp_sv_1or2[nParity * ncolor] = {}; // fixed nasty bug (omitting 'nParity' caused memory corruptions after calling calculate_jamps)
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
        // **NB! in "mixed" precision, using SIMD, calculate_jamps computes MEs for TWO neppV pages with a single channelId! #924
        calculate_jamps( allmomenta, allcouplings, jamp_sv_1or2, channelId, allNumerators, allDenominators, ievt00, ihel );
#else
        calculate_jamps( allmomenta, allcouplings, jamp_sv_1or2, ievt00, ihel );
#endif
        color_sum_cpu( allMEs, jamp_sv_1or2, ievt00 );
        MEs_ighel[ighel] = E_ACCESS::kernelAccess( E_ACCESS::ieventAccessRecord( allMEs, ievt00 ) );
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
        MEs_ighel2[ighel] = E_ACCESS::kernelAccess( E_ACCESS::ieventAccessRecord( allMEs, ievt00 + neppV ) );
#endif
        for( int iParity = 0; iParity < nParity; ++iParity )
          for( int icol = 0; icol < ncolor; icol++ )
            jamp2_sv[ncolor * iParity + icol] += cxabs2( jamp_sv_1or2[ncolor * iParity + icol] ); // may underflow #831
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
