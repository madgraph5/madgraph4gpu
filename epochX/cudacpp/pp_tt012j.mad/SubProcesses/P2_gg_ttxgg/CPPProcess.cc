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
// Process: g g > t t~ g g WEIGHTED<=4 @2

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

  constexpr int ndiagrams = CPPProcess::ndiagrams; // the number of Feynman diagrams

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

#ifdef MGONGPUCPP_GPUIMPL
  // Launch a Feynman diagram as a standalone kernel (sigmaKin_getGoodHel) or within a CUDA/HIP graph (sigmaKin)
  template<typename Func, typename... Args>
  void
  gpuDiagram( gpuGraph_t* pGraph,
              gpuGraphExec_t* pGraphExec,
              gpuGraphNode_t* pNode,
              gpuGraphNode_t* pNodeDep,
              Func diagram,
              int gpublocks,
              int gputhreads,
              gpuStream_t gpustream,
              Args... args )
  {
    // CASE 0: WITHOUT GRAPHS (sigmaKin_getGoodHel)
    if( gpustream == 0 )
    {
      gpuLaunchKernelStream( diagram, gpublocks, gputhreads, gpustream, args... );
    }
    // CASE 1: WITH GRAPHS (sigmaKin)
    else
    {
      // Define the parameters for the graph node for this Feynman diagram
      gpuKernelNodeParams params = {};
      void* kParams[] = { static_cast<void*>( &args )... };
      params.func = (void*)diagram;
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

      // *** DIAGRAMS 1 TO 123 ***
#ifdef MGONGPUCPP_GPUIMPL
      static gpuGraph_t graphs[ncomb] = {};
      static gpuGraphExec_t graphExecs[ncomb] = {};
      static gpuGraphNode_t graphNodes[ncomb * ndiagrams] = {};
      gpuGraph_t& graph = graphs[ihel];
      gpuGraphExec_t& graphExec = graphExecs[ihel];
      // Case 1 with graphs (gpustream!=0, sigmaKin): create the graph if not yet done
      if( gpustream != 0 )
      {
        if( !graph )
        {
          checkGpu( gpuGraphCreate( &graph, 0 ) );
          //std::cout << "(ihel=" << ihel << ") Created graph " << graph << std::endl;
        }
      }
      // Case 0 without graphs (gpustream==0, sigmaKin_getGoodHel): launch all diagram kernels
      // Case 1 with graphs (gpustream!=0, sigmaKin): create graph nodes if not yet done, else update them with new parameters
      gpuGraphNode_t& node1 = graphNodes[ihel * ndiagrams + 0];
      gpuDiagram( &graph, &graphExec, &node1, nullptr, diagram1, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators, momenta, ihel );
      gpuGraphNode_t& node2 = graphNodes[ihel * ndiagrams + 1];
      gpuDiagram( &graph, &graphExec, &node2, &node1, diagram2, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node3 = graphNodes[ihel * ndiagrams + 2];
      gpuDiagram( &graph, &graphExec, &node3, &node2, diagram3, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node4 = graphNodes[ihel * ndiagrams + 3];
      gpuDiagram( &graph, &graphExec, &node4, &node3, diagram4, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node5 = graphNodes[ihel * ndiagrams + 4];
      gpuDiagram( &graph, &graphExec, &node5, &node4, diagram5, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node6 = graphNodes[ihel * ndiagrams + 5];
      gpuDiagram( &graph, &graphExec, &node6, &node5, diagram6, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node7 = graphNodes[ihel * ndiagrams + 6];
      gpuDiagram( &graph, &graphExec, &node7, &node6, diagram7, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node8 = graphNodes[ihel * ndiagrams + 7];
      gpuDiagram( &graph, &graphExec, &node8, &node7, diagram8, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node9 = graphNodes[ihel * ndiagrams + 8];
      gpuDiagram( &graph, &graphExec, &node9, &node8, diagram9, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node10 = graphNodes[ihel * ndiagrams + 9];
      gpuDiagram( &graph, &graphExec, &node10, &node9, diagram10, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node11 = graphNodes[ihel * ndiagrams + 10];
      gpuDiagram( &graph, &graphExec, &node11, &node10, diagram11, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node12 = graphNodes[ihel * ndiagrams + 11];
      gpuDiagram( &graph, &graphExec, &node12, &node11, diagram12, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node13 = graphNodes[ihel * ndiagrams + 12];
      gpuDiagram( &graph, &graphExec, &node13, &node12, diagram13, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node14 = graphNodes[ihel * ndiagrams + 13];
      gpuDiagram( &graph, &graphExec, &node14, &node13, diagram14, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node15 = graphNodes[ihel * ndiagrams + 14];
      gpuDiagram( &graph, &graphExec, &node15, &node14, diagram15, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node16 = graphNodes[ihel * ndiagrams + 15];
      gpuDiagram( &graph, &graphExec, &node16, &node15, diagram16, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node17 = graphNodes[ihel * ndiagrams + 16];
      gpuDiagram( &graph, &graphExec, &node17, &node16, diagram17, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node18 = graphNodes[ihel * ndiagrams + 17];
      gpuDiagram( &graph, &graphExec, &node18, &node17, diagram18, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node19 = graphNodes[ihel * ndiagrams + 18];
      gpuDiagram( &graph, &graphExec, &node19, &node18, diagram19, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node20 = graphNodes[ihel * ndiagrams + 19];
      gpuDiagram( &graph, &graphExec, &node20, &node19, diagram20, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node21 = graphNodes[ihel * ndiagrams + 20];
      gpuDiagram( &graph, &graphExec, &node21, &node20, diagram21, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node22 = graphNodes[ihel * ndiagrams + 21];
      gpuDiagram( &graph, &graphExec, &node22, &node21, diagram22, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node23 = graphNodes[ihel * ndiagrams + 22];
      gpuDiagram( &graph, &graphExec, &node23, &node22, diagram23, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node24 = graphNodes[ihel * ndiagrams + 23];
      gpuDiagram( &graph, &graphExec, &node24, &node23, diagram24, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node25 = graphNodes[ihel * ndiagrams + 24];
      gpuDiagram( &graph, &graphExec, &node25, &node24, diagram25, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node26 = graphNodes[ihel * ndiagrams + 25];
      gpuDiagram( &graph, &graphExec, &node26, &node25, diagram26, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node27 = graphNodes[ihel * ndiagrams + 26];
      gpuDiagram( &graph, &graphExec, &node27, &node26, diagram27, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node28 = graphNodes[ihel * ndiagrams + 27];
      gpuDiagram( &graph, &graphExec, &node28, &node27, diagram28, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node29 = graphNodes[ihel * ndiagrams + 28];
      gpuDiagram( &graph, &graphExec, &node29, &node28, diagram29, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node30 = graphNodes[ihel * ndiagrams + 29];
      gpuDiagram( &graph, &graphExec, &node30, &node29, diagram30, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node31 = graphNodes[ihel * ndiagrams + 30];
      gpuDiagram( &graph, &graphExec, &node31, &node30, diagram31, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node32 = graphNodes[ihel * ndiagrams + 31];
      gpuDiagram( &graph, &graphExec, &node32, &node31, diagram32, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node33 = graphNodes[ihel * ndiagrams + 32];
      gpuDiagram( &graph, &graphExec, &node33, &node32, diagram33, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node34 = graphNodes[ihel * ndiagrams + 33];
      gpuDiagram( &graph, &graphExec, &node34, &node33, diagram34, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node35 = graphNodes[ihel * ndiagrams + 34];
      gpuDiagram( &graph, &graphExec, &node35, &node34, diagram35, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node36 = graphNodes[ihel * ndiagrams + 35];
      gpuDiagram( &graph, &graphExec, &node36, &node35, diagram36, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node37 = graphNodes[ihel * ndiagrams + 36];
      gpuDiagram( &graph, &graphExec, &node37, &node36, diagram37, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node38 = graphNodes[ihel * ndiagrams + 37];
      gpuDiagram( &graph, &graphExec, &node38, &node37, diagram38, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node39 = graphNodes[ihel * ndiagrams + 38];
      gpuDiagram( &graph, &graphExec, &node39, &node38, diagram39, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node40 = graphNodes[ihel * ndiagrams + 39];
      gpuDiagram( &graph, &graphExec, &node40, &node39, diagram40, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node41 = graphNodes[ihel * ndiagrams + 40];
      gpuDiagram( &graph, &graphExec, &node41, &node40, diagram41, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node42 = graphNodes[ihel * ndiagrams + 41];
      gpuDiagram( &graph, &graphExec, &node42, &node41, diagram42, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node43 = graphNodes[ihel * ndiagrams + 42];
      gpuDiagram( &graph, &graphExec, &node43, &node42, diagram43, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node44 = graphNodes[ihel * ndiagrams + 43];
      gpuDiagram( &graph, &graphExec, &node44, &node43, diagram44, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node45 = graphNodes[ihel * ndiagrams + 44];
      gpuDiagram( &graph, &graphExec, &node45, &node44, diagram45, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node46 = graphNodes[ihel * ndiagrams + 45];
      gpuDiagram( &graph, &graphExec, &node46, &node45, diagram46, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node47 = graphNodes[ihel * ndiagrams + 46];
      gpuDiagram( &graph, &graphExec, &node47, &node46, diagram47, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node48 = graphNodes[ihel * ndiagrams + 47];
      gpuDiagram( &graph, &graphExec, &node48, &node47, diagram48, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node49 = graphNodes[ihel * ndiagrams + 48];
      gpuDiagram( &graph, &graphExec, &node49, &node48, diagram49, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node50 = graphNodes[ihel * ndiagrams + 49];
      gpuDiagram( &graph, &graphExec, &node50, &node49, diagram50, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node51 = graphNodes[ihel * ndiagrams + 50];
      gpuDiagram( &graph, &graphExec, &node51, &node50, diagram51, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node52 = graphNodes[ihel * ndiagrams + 51];
      gpuDiagram( &graph, &graphExec, &node52, &node51, diagram52, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node53 = graphNodes[ihel * ndiagrams + 52];
      gpuDiagram( &graph, &graphExec, &node53, &node52, diagram53, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node54 = graphNodes[ihel * ndiagrams + 53];
      gpuDiagram( &graph, &graphExec, &node54, &node53, diagram54, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node55 = graphNodes[ihel * ndiagrams + 54];
      gpuDiagram( &graph, &graphExec, &node55, &node54, diagram55, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node56 = graphNodes[ihel * ndiagrams + 55];
      gpuDiagram( &graph, &graphExec, &node56, &node55, diagram56, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node57 = graphNodes[ihel * ndiagrams + 56];
      gpuDiagram( &graph, &graphExec, &node57, &node56, diagram57, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node58 = graphNodes[ihel * ndiagrams + 57];
      gpuDiagram( &graph, &graphExec, &node58, &node57, diagram58, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node59 = graphNodes[ihel * ndiagrams + 58];
      gpuDiagram( &graph, &graphExec, &node59, &node58, diagram59, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node60 = graphNodes[ihel * ndiagrams + 59];
      gpuDiagram( &graph, &graphExec, &node60, &node59, diagram60, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node61 = graphNodes[ihel * ndiagrams + 60];
      gpuDiagram( &graph, &graphExec, &node61, &node60, diagram61, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node62 = graphNodes[ihel * ndiagrams + 61];
      gpuDiagram( &graph, &graphExec, &node62, &node61, diagram62, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node63 = graphNodes[ihel * ndiagrams + 62];
      gpuDiagram( &graph, &graphExec, &node63, &node62, diagram63, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node64 = graphNodes[ihel * ndiagrams + 63];
      gpuDiagram( &graph, &graphExec, &node64, &node63, diagram64, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node65 = graphNodes[ihel * ndiagrams + 64];
      gpuDiagram( &graph, &graphExec, &node65, &node64, diagram65, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node66 = graphNodes[ihel * ndiagrams + 65];
      gpuDiagram( &graph, &graphExec, &node66, &node65, diagram66, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node67 = graphNodes[ihel * ndiagrams + 66];
      gpuDiagram( &graph, &graphExec, &node67, &node66, diagram67, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node68 = graphNodes[ihel * ndiagrams + 67];
      gpuDiagram( &graph, &graphExec, &node68, &node67, diagram68, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node69 = graphNodes[ihel * ndiagrams + 68];
      gpuDiagram( &graph, &graphExec, &node69, &node68, diagram69, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node70 = graphNodes[ihel * ndiagrams + 69];
      gpuDiagram( &graph, &graphExec, &node70, &node69, diagram70, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node71 = graphNodes[ihel * ndiagrams + 70];
      gpuDiagram( &graph, &graphExec, &node71, &node70, diagram71, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node72 = graphNodes[ihel * ndiagrams + 71];
      gpuDiagram( &graph, &graphExec, &node72, &node71, diagram72, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node73 = graphNodes[ihel * ndiagrams + 72];
      gpuDiagram( &graph, &graphExec, &node73, &node72, diagram73, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node74 = graphNodes[ihel * ndiagrams + 73];
      gpuDiagram( &graph, &graphExec, &node74, &node73, diagram74, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node75 = graphNodes[ihel * ndiagrams + 74];
      gpuDiagram( &graph, &graphExec, &node75, &node74, diagram75, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node76 = graphNodes[ihel * ndiagrams + 75];
      gpuDiagram( &graph, &graphExec, &node76, &node75, diagram76, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node77 = graphNodes[ihel * ndiagrams + 76];
      gpuDiagram( &graph, &graphExec, &node77, &node76, diagram77, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node78 = graphNodes[ihel * ndiagrams + 77];
      gpuDiagram( &graph, &graphExec, &node78, &node77, diagram78, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node79 = graphNodes[ihel * ndiagrams + 78];
      gpuDiagram( &graph, &graphExec, &node79, &node78, diagram79, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node80 = graphNodes[ihel * ndiagrams + 79];
      gpuDiagram( &graph, &graphExec, &node80, &node79, diagram80, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node81 = graphNodes[ihel * ndiagrams + 80];
      gpuDiagram( &graph, &graphExec, &node81, &node80, diagram81, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node82 = graphNodes[ihel * ndiagrams + 81];
      gpuDiagram( &graph, &graphExec, &node82, &node81, diagram82, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node83 = graphNodes[ihel * ndiagrams + 82];
      gpuDiagram( &graph, &graphExec, &node83, &node82, diagram83, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node84 = graphNodes[ihel * ndiagrams + 83];
      gpuDiagram( &graph, &graphExec, &node84, &node83, diagram84, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node85 = graphNodes[ihel * ndiagrams + 84];
      gpuDiagram( &graph, &graphExec, &node85, &node84, diagram85, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node86 = graphNodes[ihel * ndiagrams + 85];
      gpuDiagram( &graph, &graphExec, &node86, &node85, diagram86, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node87 = graphNodes[ihel * ndiagrams + 86];
      gpuDiagram( &graph, &graphExec, &node87, &node86, diagram87, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node88 = graphNodes[ihel * ndiagrams + 87];
      gpuDiagram( &graph, &graphExec, &node88, &node87, diagram88, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node89 = graphNodes[ihel * ndiagrams + 88];
      gpuDiagram( &graph, &graphExec, &node89, &node88, diagram89, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node90 = graphNodes[ihel * ndiagrams + 89];
      gpuDiagram( &graph, &graphExec, &node90, &node89, diagram90, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node91 = graphNodes[ihel * ndiagrams + 90];
      gpuDiagram( &graph, &graphExec, &node91, &node90, diagram91, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node92 = graphNodes[ihel * ndiagrams + 91];
      gpuDiagram( &graph, &graphExec, &node92, &node91, diagram92, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node93 = graphNodes[ihel * ndiagrams + 92];
      gpuDiagram( &graph, &graphExec, &node93, &node92, diagram93, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node94 = graphNodes[ihel * ndiagrams + 93];
      gpuDiagram( &graph, &graphExec, &node94, &node93, diagram94, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node95 = graphNodes[ihel * ndiagrams + 94];
      gpuDiagram( &graph, &graphExec, &node95, &node94, diagram95, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node96 = graphNodes[ihel * ndiagrams + 95];
      gpuDiagram( &graph, &graphExec, &node96, &node95, diagram96, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node97 = graphNodes[ihel * ndiagrams + 96];
      gpuDiagram( &graph, &graphExec, &node97, &node96, diagram97, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node98 = graphNodes[ihel * ndiagrams + 97];
      gpuDiagram( &graph, &graphExec, &node98, &node97, diagram98, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node99 = graphNodes[ihel * ndiagrams + 98];
      gpuDiagram( &graph, &graphExec, &node99, &node98, diagram99, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node100 = graphNodes[ihel * ndiagrams + 99];
      gpuDiagram( &graph, &graphExec, &node100, &node99, diagram100, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node101 = graphNodes[ihel * ndiagrams + 100];
      gpuDiagram( &graph, &graphExec, &node101, &node100, diagram101, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node102 = graphNodes[ihel * ndiagrams + 101];
      gpuDiagram( &graph, &graphExec, &node102, &node101, diagram102, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node103 = graphNodes[ihel * ndiagrams + 102];
      gpuDiagram( &graph, &graphExec, &node103, &node102, diagram103, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node104 = graphNodes[ihel * ndiagrams + 103];
      gpuDiagram( &graph, &graphExec, &node104, &node103, diagram104, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node105 = graphNodes[ihel * ndiagrams + 104];
      gpuDiagram( &graph, &graphExec, &node105, &node104, diagram105, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node106 = graphNodes[ihel * ndiagrams + 105];
      gpuDiagram( &graph, &graphExec, &node106, &node105, diagram106, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node107 = graphNodes[ihel * ndiagrams + 106];
      gpuDiagram( &graph, &graphExec, &node107, &node106, diagram107, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node108 = graphNodes[ihel * ndiagrams + 107];
      gpuDiagram( &graph, &graphExec, &node108, &node107, diagram108, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node109 = graphNodes[ihel * ndiagrams + 108];
      gpuDiagram( &graph, &graphExec, &node109, &node108, diagram109, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node110 = graphNodes[ihel * ndiagrams + 109];
      gpuDiagram( &graph, &graphExec, &node110, &node109, diagram110, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node111 = graphNodes[ihel * ndiagrams + 110];
      gpuDiagram( &graph, &graphExec, &node111, &node110, diagram111, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node112 = graphNodes[ihel * ndiagrams + 111];
      gpuDiagram( &graph, &graphExec, &node112, &node111, diagram112, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node113 = graphNodes[ihel * ndiagrams + 112];
      gpuDiagram( &graph, &graphExec, &node113, &node112, diagram113, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node114 = graphNodes[ihel * ndiagrams + 113];
      gpuDiagram( &graph, &graphExec, &node114, &node113, diagram114, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node115 = graphNodes[ihel * ndiagrams + 114];
      gpuDiagram( &graph, &graphExec, &node115, &node114, diagram115, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node116 = graphNodes[ihel * ndiagrams + 115];
      gpuDiagram( &graph, &graphExec, &node116, &node115, diagram116, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node117 = graphNodes[ihel * ndiagrams + 116];
      gpuDiagram( &graph, &graphExec, &node117, &node116, diagram117, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node118 = graphNodes[ihel * ndiagrams + 117];
      gpuDiagram( &graph, &graphExec, &node118, &node117, diagram118, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node119 = graphNodes[ihel * ndiagrams + 118];
      gpuDiagram( &graph, &graphExec, &node119, &node118, diagram119, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node120 = graphNodes[ihel * ndiagrams + 119];
      gpuDiagram( &graph, &graphExec, &node120, &node119, diagram120, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node121 = graphNodes[ihel * ndiagrams + 120];
      gpuDiagram( &graph, &graphExec, &node121, &node120, diagram121, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node122 = graphNodes[ihel * ndiagrams + 121];
      gpuDiagram( &graph, &graphExec, &node122, &node121, diagram122, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node123 = graphNodes[ihel * ndiagrams + 122];
      gpuDiagram( &graph, &graphExec, &node123, &node122, diagram123, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      // Case 1 with graphs (gpustream!=0, sigmaKin): create the graph executor if not yet done, then launch the graph executor
      if( gpustream != 0 )
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
      { -1, -1, -1, 1, -1, -1 },
      { -1, -1, -1, 1, -1, 1 },
      { -1, -1, -1, 1, 1, -1 },
      { -1, -1, -1, 1, 1, 1 },
      { -1, -1, -1, -1, -1, -1 },
      { -1, -1, -1, -1, -1, 1 },
      { -1, -1, -1, -1, 1, -1 },
      { -1, -1, -1, -1, 1, 1 },
      { -1, -1, 1, 1, -1, -1 },
      { -1, -1, 1, 1, -1, 1 },
      { -1, -1, 1, 1, 1, -1 },
      { -1, -1, 1, 1, 1, 1 },
      { -1, -1, 1, -1, -1, -1 },
      { -1, -1, 1, -1, -1, 1 },
      { -1, -1, 1, -1, 1, -1 },
      { -1, -1, 1, -1, 1, 1 },
      { -1, 1, -1, 1, -1, -1 },
      { -1, 1, -1, 1, -1, 1 },
      { -1, 1, -1, 1, 1, -1 },
      { -1, 1, -1, 1, 1, 1 },
      { -1, 1, -1, -1, -1, -1 },
      { -1, 1, -1, -1, -1, 1 },
      { -1, 1, -1, -1, 1, -1 },
      { -1, 1, -1, -1, 1, 1 },
      { -1, 1, 1, 1, -1, -1 },
      { -1, 1, 1, 1, -1, 1 },
      { -1, 1, 1, 1, 1, -1 },
      { -1, 1, 1, 1, 1, 1 },
      { -1, 1, 1, -1, -1, -1 },
      { -1, 1, 1, -1, -1, 1 },
      { -1, 1, 1, -1, 1, -1 },
      { -1, 1, 1, -1, 1, 1 },
      { 1, -1, -1, 1, -1, -1 },
      { 1, -1, -1, 1, -1, 1 },
      { 1, -1, -1, 1, 1, -1 },
      { 1, -1, -1, 1, 1, 1 },
      { 1, -1, -1, -1, -1, -1 },
      { 1, -1, -1, -1, -1, 1 },
      { 1, -1, -1, -1, 1, -1 },
      { 1, -1, -1, -1, 1, 1 },
      { 1, -1, 1, 1, -1, -1 },
      { 1, -1, 1, 1, -1, 1 },
      { 1, -1, 1, 1, 1, -1 },
      { 1, -1, 1, 1, 1, 1 },
      { 1, -1, 1, -1, -1, -1 },
      { 1, -1, 1, -1, -1, 1 },
      { 1, -1, 1, -1, 1, -1 },
      { 1, -1, 1, -1, 1, 1 },
      { 1, 1, -1, 1, -1, -1 },
      { 1, 1, -1, 1, -1, 1 },
      { 1, 1, -1, 1, 1, -1 },
      { 1, 1, -1, 1, 1, 1 },
      { 1, 1, -1, -1, -1, -1 },
      { 1, 1, -1, -1, -1, 1 },
      { 1, 1, -1, -1, 1, -1 },
      { 1, 1, -1, -1, 1, 1 },
      { 1, 1, 1, 1, -1, -1 },
      { 1, 1, 1, 1, -1, 1 },
      { 1, 1, 1, 1, 1, -1 },
      { 1, 1, 1, 1, 1, 1 },
      { 1, 1, 1, -1, -1, -1 },
      { 1, 1, 1, -1, -1, 1 },
      { 1, 1, 1, -1, 1, -1 },
      { 1, 1, 1, -1, 1, 1 } };
#ifdef MGONGPUCPP_GPUIMPL
    gpuMemcpyToSymbol( cHel, tHel, ncomb * npar * sizeof( short ) );
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
#ifdef MGONGPUCPP_GPUIMPL
    // Create the normalized color matrix in device memory
    createNormalizedColorMatrix();
#endif
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
#ifdef MGONGPUCPP_GPUIMPL
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
        color_sum_cpu( allMEs, jamp_sv, ievt00 );
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
    constexpr int helcolDenominators[1] = { 512 }; // assume nprocesses == 1 (#272 and #343)

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
        color_sum_cpu( allMEs, jamp_sv, ievt00 );
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
