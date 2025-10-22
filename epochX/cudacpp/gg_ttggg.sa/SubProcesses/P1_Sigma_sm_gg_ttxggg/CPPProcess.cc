// Copyright (C) 2010 The MadGraph5_aMC@NLO development team and contributors.
// Created by: J. Alwall (Oct 2010) for the MG5aMC CPP backend.
//==========================================================================
// Copyright (C) 2020-2025 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Modified by: S. Roiser (Feb 2020) for the MG5aMC CUDACPP plugin.
// Further modified by: S. Hageboeck, O. Mattelaer, S. Roiser, J. Teig, A. Valassi, Z. Wettersten (2020-2025) for the MG5aMC CUDACPP plugin.
//==========================================================================
// This file has been automatically generated for CUDA/C++ standalone by
// MadGraph5_aMC@NLO v. 3.6.3, 2025-06-12
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

  constexpr int ndiagramgroups = CPPProcess::ndiagramgroups; // the number of Feynman diagram groups

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
  static short* cHelFlat = nullptr; // symbol address
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
                   cxtype_sv* jamp_sv_1or2,           // output: jamp_sv[ncolor] (f/d) or [2*ncolor] (m) for SIMD event page(s) ievt00 and helicity ihel
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

      // *** DIAGRAMS 1 TO 1240 ***
#ifdef MGONGPUCPP_GPUIMPL
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
      gpuDiagrams( useGraphs, &graph, &graphExec, &node1, nullptr, diagramgroup1, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD, cHelFlat, momenta, ihel );
      gpuGraphNode_t& node2 = graphNodes[ihel * ndiagramgroups + 1];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node2, &node1, diagramgroup2, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node3 = graphNodes[ihel * ndiagramgroups + 2];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node3, &node2, diagramgroup3, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node4 = graphNodes[ihel * ndiagramgroups + 3];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node4, &node3, diagramgroup4, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node5 = graphNodes[ihel * ndiagramgroups + 4];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node5, &node4, diagramgroup5, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node6 = graphNodes[ihel * ndiagramgroups + 5];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node6, &node5, diagramgroup6, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node7 = graphNodes[ihel * ndiagramgroups + 6];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node7, &node6, diagramgroup7, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node8 = graphNodes[ihel * ndiagramgroups + 7];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node8, &node7, diagramgroup8, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node9 = graphNodes[ihel * ndiagramgroups + 8];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node9, &node8, diagramgroup9, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node10 = graphNodes[ihel * ndiagramgroups + 9];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node10, &node9, diagramgroup10, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node11 = graphNodes[ihel * ndiagramgroups + 10];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node11, &node10, diagramgroup11, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node12 = graphNodes[ihel * ndiagramgroups + 11];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node12, &node11, diagramgroup12, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node13 = graphNodes[ihel * ndiagramgroups + 12];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node13, &node12, diagramgroup13, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node14 = graphNodes[ihel * ndiagramgroups + 13];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node14, &node13, diagramgroup14, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node15 = graphNodes[ihel * ndiagramgroups + 14];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node15, &node14, diagramgroup15, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node16 = graphNodes[ihel * ndiagramgroups + 15];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node16, &node15, diagramgroup16, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node17 = graphNodes[ihel * ndiagramgroups + 16];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node17, &node16, diagramgroup17, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node18 = graphNodes[ihel * ndiagramgroups + 17];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node18, &node17, diagramgroup18, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node19 = graphNodes[ihel * ndiagramgroups + 18];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node19, &node18, diagramgroup19, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node20 = graphNodes[ihel * ndiagramgroups + 19];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node20, &node19, diagramgroup20, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node21 = graphNodes[ihel * ndiagramgroups + 20];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node21, &node20, diagramgroup21, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node22 = graphNodes[ihel * ndiagramgroups + 21];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node22, &node21, diagramgroup22, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node23 = graphNodes[ihel * ndiagramgroups + 22];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node23, &node22, diagramgroup23, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node24 = graphNodes[ihel * ndiagramgroups + 23];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node24, &node23, diagramgroup24, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node25 = graphNodes[ihel * ndiagramgroups + 24];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node25, &node24, diagramgroup25, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node26 = graphNodes[ihel * ndiagramgroups + 25];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node26, &node25, diagramgroup26, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node27 = graphNodes[ihel * ndiagramgroups + 26];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node27, &node26, diagramgroup27, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node28 = graphNodes[ihel * ndiagramgroups + 27];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node28, &node27, diagramgroup28, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node29 = graphNodes[ihel * ndiagramgroups + 28];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node29, &node28, diagramgroup29, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node30 = graphNodes[ihel * ndiagramgroups + 29];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node30, &node29, diagramgroup30, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node31 = graphNodes[ihel * ndiagramgroups + 30];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node31, &node30, diagramgroup31, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node32 = graphNodes[ihel * ndiagramgroups + 31];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node32, &node31, diagramgroup32, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node33 = graphNodes[ihel * ndiagramgroups + 32];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node33, &node32, diagramgroup33, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node34 = graphNodes[ihel * ndiagramgroups + 33];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node34, &node33, diagramgroup34, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node35 = graphNodes[ihel * ndiagramgroups + 34];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node35, &node34, diagramgroup35, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node36 = graphNodes[ihel * ndiagramgroups + 35];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node36, &node35, diagramgroup36, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node37 = graphNodes[ihel * ndiagramgroups + 36];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node37, &node36, diagramgroup37, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node38 = graphNodes[ihel * ndiagramgroups + 37];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node38, &node37, diagramgroup38, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node39 = graphNodes[ihel * ndiagramgroups + 38];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node39, &node38, diagramgroup39, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node40 = graphNodes[ihel * ndiagramgroups + 39];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node40, &node39, diagramgroup40, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node41 = graphNodes[ihel * ndiagramgroups + 40];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node41, &node40, diagramgroup41, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node42 = graphNodes[ihel * ndiagramgroups + 41];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node42, &node41, diagramgroup42, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node43 = graphNodes[ihel * ndiagramgroups + 42];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node43, &node42, diagramgroup43, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node44 = graphNodes[ihel * ndiagramgroups + 43];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node44, &node43, diagramgroup44, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node45 = graphNodes[ihel * ndiagramgroups + 44];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node45, &node44, diagramgroup45, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node46 = graphNodes[ihel * ndiagramgroups + 45];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node46, &node45, diagramgroup46, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node47 = graphNodes[ihel * ndiagramgroups + 46];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node47, &node46, diagramgroup47, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node48 = graphNodes[ihel * ndiagramgroups + 47];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node48, &node47, diagramgroup48, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node49 = graphNodes[ihel * ndiagramgroups + 48];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node49, &node48, diagramgroup49, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node50 = graphNodes[ihel * ndiagramgroups + 49];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node50, &node49, diagramgroup50, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node51 = graphNodes[ihel * ndiagramgroups + 50];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node51, &node50, diagramgroup51, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node52 = graphNodes[ihel * ndiagramgroups + 51];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node52, &node51, diagramgroup52, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node53 = graphNodes[ihel * ndiagramgroups + 52];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node53, &node52, diagramgroup53, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node54 = graphNodes[ihel * ndiagramgroups + 53];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node54, &node53, diagramgroup54, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node55 = graphNodes[ihel * ndiagramgroups + 54];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node55, &node54, diagramgroup55, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node56 = graphNodes[ihel * ndiagramgroups + 55];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node56, &node55, diagramgroup56, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node57 = graphNodes[ihel * ndiagramgroups + 56];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node57, &node56, diagramgroup57, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node58 = graphNodes[ihel * ndiagramgroups + 57];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node58, &node57, diagramgroup58, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node59 = graphNodes[ihel * ndiagramgroups + 58];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node59, &node58, diagramgroup59, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node60 = graphNodes[ihel * ndiagramgroups + 59];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node60, &node59, diagramgroup60, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node61 = graphNodes[ihel * ndiagramgroups + 60];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node61, &node60, diagramgroup61, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node62 = graphNodes[ihel * ndiagramgroups + 61];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node62, &node61, diagramgroup62, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node63 = graphNodes[ihel * ndiagramgroups + 62];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node63, &node62, diagramgroup63, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node64 = graphNodes[ihel * ndiagramgroups + 63];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node64, &node63, diagramgroup64, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node65 = graphNodes[ihel * ndiagramgroups + 64];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node65, &node64, diagramgroup65, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node66 = graphNodes[ihel * ndiagramgroups + 65];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node66, &node65, diagramgroup66, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node67 = graphNodes[ihel * ndiagramgroups + 66];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node67, &node66, diagramgroup67, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node68 = graphNodes[ihel * ndiagramgroups + 67];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node68, &node67, diagramgroup68, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node69 = graphNodes[ihel * ndiagramgroups + 68];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node69, &node68, diagramgroup69, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node70 = graphNodes[ihel * ndiagramgroups + 69];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node70, &node69, diagramgroup70, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node71 = graphNodes[ihel * ndiagramgroups + 70];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node71, &node70, diagramgroup71, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node72 = graphNodes[ihel * ndiagramgroups + 71];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node72, &node71, diagramgroup72, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node73 = graphNodes[ihel * ndiagramgroups + 72];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node73, &node72, diagramgroup73, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node74 = graphNodes[ihel * ndiagramgroups + 73];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node74, &node73, diagramgroup74, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node75 = graphNodes[ihel * ndiagramgroups + 74];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node75, &node74, diagramgroup75, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node76 = graphNodes[ihel * ndiagramgroups + 75];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node76, &node75, diagramgroup76, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node77 = graphNodes[ihel * ndiagramgroups + 76];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node77, &node76, diagramgroup77, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node78 = graphNodes[ihel * ndiagramgroups + 77];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node78, &node77, diagramgroup78, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node79 = graphNodes[ihel * ndiagramgroups + 78];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node79, &node78, diagramgroup79, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node80 = graphNodes[ihel * ndiagramgroups + 79];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node80, &node79, diagramgroup80, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node81 = graphNodes[ihel * ndiagramgroups + 80];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node81, &node80, diagramgroup81, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node82 = graphNodes[ihel * ndiagramgroups + 81];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node82, &node81, diagramgroup82, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node83 = graphNodes[ihel * ndiagramgroups + 82];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node83, &node82, diagramgroup83, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node84 = graphNodes[ihel * ndiagramgroups + 83];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node84, &node83, diagramgroup84, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node85 = graphNodes[ihel * ndiagramgroups + 84];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node85, &node84, diagramgroup85, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node86 = graphNodes[ihel * ndiagramgroups + 85];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node86, &node85, diagramgroup86, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node87 = graphNodes[ihel * ndiagramgroups + 86];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node87, &node86, diagramgroup87, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node88 = graphNodes[ihel * ndiagramgroups + 87];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node88, &node87, diagramgroup88, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node89 = graphNodes[ihel * ndiagramgroups + 88];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node89, &node88, diagramgroup89, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node90 = graphNodes[ihel * ndiagramgroups + 89];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node90, &node89, diagramgroup90, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node91 = graphNodes[ihel * ndiagramgroups + 90];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node91, &node90, diagramgroup91, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node92 = graphNodes[ihel * ndiagramgroups + 91];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node92, &node91, diagramgroup92, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node93 = graphNodes[ihel * ndiagramgroups + 92];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node93, &node92, diagramgroup93, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node94 = graphNodes[ihel * ndiagramgroups + 93];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node94, &node93, diagramgroup94, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node95 = graphNodes[ihel * ndiagramgroups + 94];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node95, &node94, diagramgroup95, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node96 = graphNodes[ihel * ndiagramgroups + 95];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node96, &node95, diagramgroup96, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node97 = graphNodes[ihel * ndiagramgroups + 96];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node97, &node96, diagramgroup97, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node98 = graphNodes[ihel * ndiagramgroups + 97];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node98, &node97, diagramgroup98, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node99 = graphNodes[ihel * ndiagramgroups + 98];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node99, &node98, diagramgroup99, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node100 = graphNodes[ihel * ndiagramgroups + 99];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node100, &node99, diagramgroup100, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node101 = graphNodes[ihel * ndiagramgroups + 100];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node101, &node100, diagramgroup101, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node102 = graphNodes[ihel * ndiagramgroups + 101];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node102, &node101, diagramgroup102, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node103 = graphNodes[ihel * ndiagramgroups + 102];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node103, &node102, diagramgroup103, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node104 = graphNodes[ihel * ndiagramgroups + 103];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node104, &node103, diagramgroup104, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node105 = graphNodes[ihel * ndiagramgroups + 104];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node105, &node104, diagramgroup105, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node106 = graphNodes[ihel * ndiagramgroups + 105];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node106, &node105, diagramgroup106, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node107 = graphNodes[ihel * ndiagramgroups + 106];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node107, &node106, diagramgroup107, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node108 = graphNodes[ihel * ndiagramgroups + 107];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node108, &node107, diagramgroup108, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node109 = graphNodes[ihel * ndiagramgroups + 108];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node109, &node108, diagramgroup109, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node110 = graphNodes[ihel * ndiagramgroups + 109];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node110, &node109, diagramgroup110, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node111 = graphNodes[ihel * ndiagramgroups + 110];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node111, &node110, diagramgroup111, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node112 = graphNodes[ihel * ndiagramgroups + 111];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node112, &node111, diagramgroup112, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node113 = graphNodes[ihel * ndiagramgroups + 112];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node113, &node112, diagramgroup113, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node114 = graphNodes[ihel * ndiagramgroups + 113];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node114, &node113, diagramgroup114, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node115 = graphNodes[ihel * ndiagramgroups + 114];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node115, &node114, diagramgroup115, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node116 = graphNodes[ihel * ndiagramgroups + 115];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node116, &node115, diagramgroup116, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node117 = graphNodes[ihel * ndiagramgroups + 116];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node117, &node116, diagramgroup117, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node118 = graphNodes[ihel * ndiagramgroups + 117];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node118, &node117, diagramgroup118, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node119 = graphNodes[ihel * ndiagramgroups + 118];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node119, &node118, diagramgroup119, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node120 = graphNodes[ihel * ndiagramgroups + 119];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node120, &node119, diagramgroup120, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node121 = graphNodes[ihel * ndiagramgroups + 120];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node121, &node120, diagramgroup121, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node122 = graphNodes[ihel * ndiagramgroups + 121];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node122, &node121, diagramgroup122, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node123 = graphNodes[ihel * ndiagramgroups + 122];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node123, &node122, diagramgroup123, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node124 = graphNodes[ihel * ndiagramgroups + 123];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node124, &node123, diagramgroup124, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node125 = graphNodes[ihel * ndiagramgroups + 124];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node125, &node124, diagramgroup125, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node126 = graphNodes[ihel * ndiagramgroups + 125];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node126, &node125, diagramgroup126, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node127 = graphNodes[ihel * ndiagramgroups + 126];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node127, &node126, diagramgroup127, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node128 = graphNodes[ihel * ndiagramgroups + 127];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node128, &node127, diagramgroup128, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node129 = graphNodes[ihel * ndiagramgroups + 128];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node129, &node128, diagramgroup129, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node130 = graphNodes[ihel * ndiagramgroups + 129];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node130, &node129, diagramgroup130, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node131 = graphNodes[ihel * ndiagramgroups + 130];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node131, &node130, diagramgroup131, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node132 = graphNodes[ihel * ndiagramgroups + 131];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node132, &node131, diagramgroup132, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node133 = graphNodes[ihel * ndiagramgroups + 132];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node133, &node132, diagramgroup133, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node134 = graphNodes[ihel * ndiagramgroups + 133];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node134, &node133, diagramgroup134, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node135 = graphNodes[ihel * ndiagramgroups + 134];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node135, &node134, diagramgroup135, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node136 = graphNodes[ihel * ndiagramgroups + 135];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node136, &node135, diagramgroup136, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node137 = graphNodes[ihel * ndiagramgroups + 136];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node137, &node136, diagramgroup137, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node138 = graphNodes[ihel * ndiagramgroups + 137];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node138, &node137, diagramgroup138, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node139 = graphNodes[ihel * ndiagramgroups + 138];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node139, &node138, diagramgroup139, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node140 = graphNodes[ihel * ndiagramgroups + 139];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node140, &node139, diagramgroup140, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node141 = graphNodes[ihel * ndiagramgroups + 140];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node141, &node140, diagramgroup141, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node142 = graphNodes[ihel * ndiagramgroups + 141];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node142, &node141, diagramgroup142, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node143 = graphNodes[ihel * ndiagramgroups + 142];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node143, &node142, diagramgroup143, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node144 = graphNodes[ihel * ndiagramgroups + 143];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node144, &node143, diagramgroup144, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node145 = graphNodes[ihel * ndiagramgroups + 144];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node145, &node144, diagramgroup145, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node146 = graphNodes[ihel * ndiagramgroups + 145];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node146, &node145, diagramgroup146, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node147 = graphNodes[ihel * ndiagramgroups + 146];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node147, &node146, diagramgroup147, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node148 = graphNodes[ihel * ndiagramgroups + 147];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node148, &node147, diagramgroup148, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node149 = graphNodes[ihel * ndiagramgroups + 148];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node149, &node148, diagramgroup149, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node150 = graphNodes[ihel * ndiagramgroups + 149];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node150, &node149, diagramgroup150, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node151 = graphNodes[ihel * ndiagramgroups + 150];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node151, &node150, diagramgroup151, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node152 = graphNodes[ihel * ndiagramgroups + 151];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node152, &node151, diagramgroup152, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node153 = graphNodes[ihel * ndiagramgroups + 152];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node153, &node152, diagramgroup153, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node154 = graphNodes[ihel * ndiagramgroups + 153];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node154, &node153, diagramgroup154, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node155 = graphNodes[ihel * ndiagramgroups + 154];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node155, &node154, diagramgroup155, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node156 = graphNodes[ihel * ndiagramgroups + 155];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node156, &node155, diagramgroup156, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node157 = graphNodes[ihel * ndiagramgroups + 156];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node157, &node156, diagramgroup157, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node158 = graphNodes[ihel * ndiagramgroups + 157];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node158, &node157, diagramgroup158, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node159 = graphNodes[ihel * ndiagramgroups + 158];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node159, &node158, diagramgroup159, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node160 = graphNodes[ihel * ndiagramgroups + 159];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node160, &node159, diagramgroup160, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node161 = graphNodes[ihel * ndiagramgroups + 160];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node161, &node160, diagramgroup161, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node162 = graphNodes[ihel * ndiagramgroups + 161];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node162, &node161, diagramgroup162, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node163 = graphNodes[ihel * ndiagramgroups + 162];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node163, &node162, diagramgroup163, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node164 = graphNodes[ihel * ndiagramgroups + 163];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node164, &node163, diagramgroup164, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node165 = graphNodes[ihel * ndiagramgroups + 164];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node165, &node164, diagramgroup165, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node166 = graphNodes[ihel * ndiagramgroups + 165];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node166, &node165, diagramgroup166, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node167 = graphNodes[ihel * ndiagramgroups + 166];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node167, &node166, diagramgroup167, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node168 = graphNodes[ihel * ndiagramgroups + 167];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node168, &node167, diagramgroup168, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node169 = graphNodes[ihel * ndiagramgroups + 168];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node169, &node168, diagramgroup169, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node170 = graphNodes[ihel * ndiagramgroups + 169];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node170, &node169, diagramgroup170, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node171 = graphNodes[ihel * ndiagramgroups + 170];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node171, &node170, diagramgroup171, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node172 = graphNodes[ihel * ndiagramgroups + 171];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node172, &node171, diagramgroup172, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node173 = graphNodes[ihel * ndiagramgroups + 172];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node173, &node172, diagramgroup173, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node174 = graphNodes[ihel * ndiagramgroups + 173];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node174, &node173, diagramgroup174, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node175 = graphNodes[ihel * ndiagramgroups + 174];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node175, &node174, diagramgroup175, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node176 = graphNodes[ihel * ndiagramgroups + 175];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node176, &node175, diagramgroup176, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node177 = graphNodes[ihel * ndiagramgroups + 176];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node177, &node176, diagramgroup177, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node178 = graphNodes[ihel * ndiagramgroups + 177];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node178, &node177, diagramgroup178, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node179 = graphNodes[ihel * ndiagramgroups + 178];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node179, &node178, diagramgroup179, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node180 = graphNodes[ihel * ndiagramgroups + 179];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node180, &node179, diagramgroup180, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node181 = graphNodes[ihel * ndiagramgroups + 180];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node181, &node180, diagramgroup181, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node182 = graphNodes[ihel * ndiagramgroups + 181];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node182, &node181, diagramgroup182, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node183 = graphNodes[ihel * ndiagramgroups + 182];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node183, &node182, diagramgroup183, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node184 = graphNodes[ihel * ndiagramgroups + 183];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node184, &node183, diagramgroup184, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node185 = graphNodes[ihel * ndiagramgroups + 184];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node185, &node184, diagramgroup185, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node186 = graphNodes[ihel * ndiagramgroups + 185];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node186, &node185, diagramgroup186, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node187 = graphNodes[ihel * ndiagramgroups + 186];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node187, &node186, diagramgroup187, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node188 = graphNodes[ihel * ndiagramgroups + 187];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node188, &node187, diagramgroup188, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node189 = graphNodes[ihel * ndiagramgroups + 188];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node189, &node188, diagramgroup189, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node190 = graphNodes[ihel * ndiagramgroups + 189];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node190, &node189, diagramgroup190, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node191 = graphNodes[ihel * ndiagramgroups + 190];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node191, &node190, diagramgroup191, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node192 = graphNodes[ihel * ndiagramgroups + 191];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node192, &node191, diagramgroup192, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node193 = graphNodes[ihel * ndiagramgroups + 192];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node193, &node192, diagramgroup193, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node194 = graphNodes[ihel * ndiagramgroups + 193];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node194, &node193, diagramgroup194, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node195 = graphNodes[ihel * ndiagramgroups + 194];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node195, &node194, diagramgroup195, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node196 = graphNodes[ihel * ndiagramgroups + 195];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node196, &node195, diagramgroup196, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node197 = graphNodes[ihel * ndiagramgroups + 196];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node197, &node196, diagramgroup197, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node198 = graphNodes[ihel * ndiagramgroups + 197];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node198, &node197, diagramgroup198, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node199 = graphNodes[ihel * ndiagramgroups + 198];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node199, &node198, diagramgroup199, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node200 = graphNodes[ihel * ndiagramgroups + 199];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node200, &node199, diagramgroup200, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node201 = graphNodes[ihel * ndiagramgroups + 200];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node201, &node200, diagramgroup201, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node202 = graphNodes[ihel * ndiagramgroups + 201];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node202, &node201, diagramgroup202, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node203 = graphNodes[ihel * ndiagramgroups + 202];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node203, &node202, diagramgroup203, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node204 = graphNodes[ihel * ndiagramgroups + 203];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node204, &node203, diagramgroup204, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node205 = graphNodes[ihel * ndiagramgroups + 204];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node205, &node204, diagramgroup205, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node206 = graphNodes[ihel * ndiagramgroups + 205];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node206, &node205, diagramgroup206, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node207 = graphNodes[ihel * ndiagramgroups + 206];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node207, &node206, diagramgroup207, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node208 = graphNodes[ihel * ndiagramgroups + 207];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node208, &node207, diagramgroup208, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node209 = graphNodes[ihel * ndiagramgroups + 208];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node209, &node208, diagramgroup209, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node210 = graphNodes[ihel * ndiagramgroups + 209];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node210, &node209, diagramgroup210, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node211 = graphNodes[ihel * ndiagramgroups + 210];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node211, &node210, diagramgroup211, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node212 = graphNodes[ihel * ndiagramgroups + 211];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node212, &node211, diagramgroup212, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node213 = graphNodes[ihel * ndiagramgroups + 212];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node213, &node212, diagramgroup213, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node214 = graphNodes[ihel * ndiagramgroups + 213];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node214, &node213, diagramgroup214, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node215 = graphNodes[ihel * ndiagramgroups + 214];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node215, &node214, diagramgroup215, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node216 = graphNodes[ihel * ndiagramgroups + 215];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node216, &node215, diagramgroup216, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node217 = graphNodes[ihel * ndiagramgroups + 216];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node217, &node216, diagramgroup217, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node218 = graphNodes[ihel * ndiagramgroups + 217];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node218, &node217, diagramgroup218, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node219 = graphNodes[ihel * ndiagramgroups + 218];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node219, &node218, diagramgroup219, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node220 = graphNodes[ihel * ndiagramgroups + 219];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node220, &node219, diagramgroup220, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node221 = graphNodes[ihel * ndiagramgroups + 220];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node221, &node220, diagramgroup221, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node222 = graphNodes[ihel * ndiagramgroups + 221];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node222, &node221, diagramgroup222, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node223 = graphNodes[ihel * ndiagramgroups + 222];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node223, &node222, diagramgroup223, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node224 = graphNodes[ihel * ndiagramgroups + 223];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node224, &node223, diagramgroup224, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node225 = graphNodes[ihel * ndiagramgroups + 224];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node225, &node224, diagramgroup225, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node226 = graphNodes[ihel * ndiagramgroups + 225];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node226, &node225, diagramgroup226, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node227 = graphNodes[ihel * ndiagramgroups + 226];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node227, &node226, diagramgroup227, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node228 = graphNodes[ihel * ndiagramgroups + 227];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node228, &node227, diagramgroup228, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node229 = graphNodes[ihel * ndiagramgroups + 228];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node229, &node228, diagramgroup229, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node230 = graphNodes[ihel * ndiagramgroups + 229];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node230, &node229, diagramgroup230, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node231 = graphNodes[ihel * ndiagramgroups + 230];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node231, &node230, diagramgroup231, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node232 = graphNodes[ihel * ndiagramgroups + 231];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node232, &node231, diagramgroup232, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node233 = graphNodes[ihel * ndiagramgroups + 232];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node233, &node232, diagramgroup233, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node234 = graphNodes[ihel * ndiagramgroups + 233];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node234, &node233, diagramgroup234, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node235 = graphNodes[ihel * ndiagramgroups + 234];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node235, &node234, diagramgroup235, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node236 = graphNodes[ihel * ndiagramgroups + 235];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node236, &node235, diagramgroup236, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node237 = graphNodes[ihel * ndiagramgroups + 236];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node237, &node236, diagramgroup237, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node238 = graphNodes[ihel * ndiagramgroups + 237];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node238, &node237, diagramgroup238, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node239 = graphNodes[ihel * ndiagramgroups + 238];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node239, &node238, diagramgroup239, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node240 = graphNodes[ihel * ndiagramgroups + 239];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node240, &node239, diagramgroup240, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node241 = graphNodes[ihel * ndiagramgroups + 240];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node241, &node240, diagramgroup241, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node242 = graphNodes[ihel * ndiagramgroups + 241];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node242, &node241, diagramgroup242, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node243 = graphNodes[ihel * ndiagramgroups + 242];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node243, &node242, diagramgroup243, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node244 = graphNodes[ihel * ndiagramgroups + 243];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node244, &node243, diagramgroup244, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node245 = graphNodes[ihel * ndiagramgroups + 244];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node245, &node244, diagramgroup245, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node246 = graphNodes[ihel * ndiagramgroups + 245];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node246, &node245, diagramgroup246, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node247 = graphNodes[ihel * ndiagramgroups + 246];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node247, &node246, diagramgroup247, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
      gpuGraphNode_t& node248 = graphNodes[ihel * ndiagramgroups + 247];
      gpuDiagrams( useGraphs, &graph, &graphExec, &node248, &node247, diagramgroup248, gpublocks, gputhreads, gpustream, wfs, jamps, couplings, channelIds, numerators, denominators, cIPC, cIPD );
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
    gpuMemcpyToSymbol( dcHel, tHel, ncomb * npar * sizeof( short ) );
    gpuGetSymbolAddress( (void**)( &cHelFlat ), dcHel );
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
#ifdef MGONGPUCPP_GPUIMPL
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
      color_sum_gpu( allMEs, allJamps, nullptr, nullptr, nullptr, gpublocks, gputhreads );
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
        calculate_jamps( ihel, allmomenta, allcouplings, jamp_sv_1or2, channelId, allNumerators, allDenominators, ievt00 ); //maxtry?
#else
        calculate_jamps( ihel, allmomenta, allcouplings, jamp_sv_1or2, ievt00 ); //maxtry?
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
            gpuBlasHandle_t* ghelBlasHandles,   // input: cuBLAS/hipBLAS handles (index is ighel: only the first nGoodHel <= ncomb are non-null)
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
    // (1b) Then, in multichannel mode, also compute the running sums over helicities of squared jamp2s within each helicity stream
    for( int ighel = 0; ighel < cNGoodHel; ighel++ )
    {
      fptype* hAllJamps = ghelAllJamps + ighel * nevt * ncolor * mgOnGpu::nx2;
      gpuLaunchKernelStream( update_jamp2s, gpublocks, gputhreads, ghelStreams[ighel], hAllJamps, colAllJamp2s );
    }
#endif
    // (2) Then, within each helicity stream, compute the ME for that helicity from the color sum of QCD partial amplitudes jamps
    if( !ghelBlasHandles )
      assert( ghelAllBlasTmp == nullptr ); // HASBLAS=hasNoBlas or CUDACPP_RUNTIME_BLASCOLORSUM not set
    else
      assert( ghelAllBlasTmp != nullptr ); // note: this should never happen for HASBLAS=hasNoBlas (a sanity check is in color_sum_gpu)
    for( int ighel = 0; ighel < cNGoodHel; ighel++ )
    {
      fptype* hAllMEs = ghelAllMEs + ighel * nevt;
      fptype* hAllJamps = ghelAllJamps + ighel * nevt * ncolor * mgOnGpu::nx2;
#if defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
      fptype2* hAllBlasTmp = ( ghelAllBlasTmp != nullptr ? ghelAllBlasTmp + ighel * nevt * ( 2 * ncolor * mgOnGpu::nx2 + 1 ) : nullptr );
      if( hAllBlasTmp )
        gpuMemset( hAllBlasTmp, 0, nevt * ( 2 * ncolor * mgOnGpu::nx2 + 1 ) * sizeof( fptype2 ) ); // reset the tmp buffer (bug fix: reset MEs=0)
#else
      fptype2* hAllBlasTmp = ( ghelAllBlasTmp != nullptr ? ghelAllBlasTmp + ighel * nevt * ncolor * mgOnGpu::nx2 : nullptr );
      if( hAllBlasTmp )
        gpuMemset( hAllBlasTmp, 0, nevt * ncolor * mgOnGpu::nx2 * sizeof( fptype2 ) ); // reset the tmp buffer (just in case...)
#endif
#ifndef MGONGPU_HAS_NO_BLAS
      gpuBlasHandle_t* pBlasHandle = ( ghelBlasHandles ? &( ghelBlasHandles[ighel] ) : nullptr );
#else /* clang-format off */
      assert( ghelBlasHandles == nullptr ); // sanity check
      gpuBlasHandle_t* pBlasHandle = nullptr; // this is a void* (hack to keep the same API in noBLAS builds)
#endif /* clang-format on */
      color_sum_gpu( hAllMEs, hAllJamps, hAllBlasTmp, ghelStreams[ighel], pBlasHandle, gpublocks, gputhreads );
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
        calculate_jamps( ihel, allmomenta, allcouplings, jamp_sv_1or2, channelId, allNumerators, allDenominators, ievt00 );
#else
        calculate_jamps( ihel, allmomenta, allcouplings, jamp_sv_1or2, ievt00 );
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
