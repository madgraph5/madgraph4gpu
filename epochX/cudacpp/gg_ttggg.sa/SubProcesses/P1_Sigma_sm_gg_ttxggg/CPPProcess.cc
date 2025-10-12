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

      // *** DIAGRAMS 1 TO 1240 ***
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
      gpuGraphNode_t& node124 = graphNodes[ihel * ndiagrams + 123];
      gpuDiagram( &graph, &graphExec, &node124, &node123, diagram124, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node125 = graphNodes[ihel * ndiagrams + 124];
      gpuDiagram( &graph, &graphExec, &node125, &node124, diagram125, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node126 = graphNodes[ihel * ndiagrams + 125];
      gpuDiagram( &graph, &graphExec, &node126, &node125, diagram126, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node127 = graphNodes[ihel * ndiagrams + 126];
      gpuDiagram( &graph, &graphExec, &node127, &node126, diagram127, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node128 = graphNodes[ihel * ndiagrams + 127];
      gpuDiagram( &graph, &graphExec, &node128, &node127, diagram128, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node129 = graphNodes[ihel * ndiagrams + 128];
      gpuDiagram( &graph, &graphExec, &node129, &node128, diagram129, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node130 = graphNodes[ihel * ndiagrams + 129];
      gpuDiagram( &graph, &graphExec, &node130, &node129, diagram130, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node131 = graphNodes[ihel * ndiagrams + 130];
      gpuDiagram( &graph, &graphExec, &node131, &node130, diagram131, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node132 = graphNodes[ihel * ndiagrams + 131];
      gpuDiagram( &graph, &graphExec, &node132, &node131, diagram132, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node133 = graphNodes[ihel * ndiagrams + 132];
      gpuDiagram( &graph, &graphExec, &node133, &node132, diagram133, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node134 = graphNodes[ihel * ndiagrams + 133];
      gpuDiagram( &graph, &graphExec, &node134, &node133, diagram134, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node135 = graphNodes[ihel * ndiagrams + 134];
      gpuDiagram( &graph, &graphExec, &node135, &node134, diagram135, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node136 = graphNodes[ihel * ndiagrams + 135];
      gpuDiagram( &graph, &graphExec, &node136, &node135, diagram136, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node137 = graphNodes[ihel * ndiagrams + 136];
      gpuDiagram( &graph, &graphExec, &node137, &node136, diagram137, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node138 = graphNodes[ihel * ndiagrams + 137];
      gpuDiagram( &graph, &graphExec, &node138, &node137, diagram138, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node139 = graphNodes[ihel * ndiagrams + 138];
      gpuDiagram( &graph, &graphExec, &node139, &node138, diagram139, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node140 = graphNodes[ihel * ndiagrams + 139];
      gpuDiagram( &graph, &graphExec, &node140, &node139, diagram140, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node141 = graphNodes[ihel * ndiagrams + 140];
      gpuDiagram( &graph, &graphExec, &node141, &node140, diagram141, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node142 = graphNodes[ihel * ndiagrams + 141];
      gpuDiagram( &graph, &graphExec, &node142, &node141, diagram142, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node143 = graphNodes[ihel * ndiagrams + 142];
      gpuDiagram( &graph, &graphExec, &node143, &node142, diagram143, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node144 = graphNodes[ihel * ndiagrams + 143];
      gpuDiagram( &graph, &graphExec, &node144, &node143, diagram144, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node145 = graphNodes[ihel * ndiagrams + 144];
      gpuDiagram( &graph, &graphExec, &node145, &node144, diagram145, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node146 = graphNodes[ihel * ndiagrams + 145];
      gpuDiagram( &graph, &graphExec, &node146, &node145, diagram146, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node147 = graphNodes[ihel * ndiagrams + 146];
      gpuDiagram( &graph, &graphExec, &node147, &node146, diagram147, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node148 = graphNodes[ihel * ndiagrams + 147];
      gpuDiagram( &graph, &graphExec, &node148, &node147, diagram148, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node149 = graphNodes[ihel * ndiagrams + 148];
      gpuDiagram( &graph, &graphExec, &node149, &node148, diagram149, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node150 = graphNodes[ihel * ndiagrams + 149];
      gpuDiagram( &graph, &graphExec, &node150, &node149, diagram150, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node151 = graphNodes[ihel * ndiagrams + 150];
      gpuDiagram( &graph, &graphExec, &node151, &node150, diagram151, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node152 = graphNodes[ihel * ndiagrams + 151];
      gpuDiagram( &graph, &graphExec, &node152, &node151, diagram152, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node153 = graphNodes[ihel * ndiagrams + 152];
      gpuDiagram( &graph, &graphExec, &node153, &node152, diagram153, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node154 = graphNodes[ihel * ndiagrams + 153];
      gpuDiagram( &graph, &graphExec, &node154, &node153, diagram154, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node155 = graphNodes[ihel * ndiagrams + 154];
      gpuDiagram( &graph, &graphExec, &node155, &node154, diagram155, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node156 = graphNodes[ihel * ndiagrams + 155];
      gpuDiagram( &graph, &graphExec, &node156, &node155, diagram156, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node157 = graphNodes[ihel * ndiagrams + 156];
      gpuDiagram( &graph, &graphExec, &node157, &node156, diagram157, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node158 = graphNodes[ihel * ndiagrams + 157];
      gpuDiagram( &graph, &graphExec, &node158, &node157, diagram158, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node159 = graphNodes[ihel * ndiagrams + 158];
      gpuDiagram( &graph, &graphExec, &node159, &node158, diagram159, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node160 = graphNodes[ihel * ndiagrams + 159];
      gpuDiagram( &graph, &graphExec, &node160, &node159, diagram160, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node161 = graphNodes[ihel * ndiagrams + 160];
      gpuDiagram( &graph, &graphExec, &node161, &node160, diagram161, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node162 = graphNodes[ihel * ndiagrams + 161];
      gpuDiagram( &graph, &graphExec, &node162, &node161, diagram162, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node163 = graphNodes[ihel * ndiagrams + 162];
      gpuDiagram( &graph, &graphExec, &node163, &node162, diagram163, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node164 = graphNodes[ihel * ndiagrams + 163];
      gpuDiagram( &graph, &graphExec, &node164, &node163, diagram164, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node165 = graphNodes[ihel * ndiagrams + 164];
      gpuDiagram( &graph, &graphExec, &node165, &node164, diagram165, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node166 = graphNodes[ihel * ndiagrams + 165];
      gpuDiagram( &graph, &graphExec, &node166, &node165, diagram166, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node167 = graphNodes[ihel * ndiagrams + 166];
      gpuDiagram( &graph, &graphExec, &node167, &node166, diagram167, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node168 = graphNodes[ihel * ndiagrams + 167];
      gpuDiagram( &graph, &graphExec, &node168, &node167, diagram168, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node169 = graphNodes[ihel * ndiagrams + 168];
      gpuDiagram( &graph, &graphExec, &node169, &node168, diagram169, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node170 = graphNodes[ihel * ndiagrams + 169];
      gpuDiagram( &graph, &graphExec, &node170, &node169, diagram170, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node171 = graphNodes[ihel * ndiagrams + 170];
      gpuDiagram( &graph, &graphExec, &node171, &node170, diagram171, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node172 = graphNodes[ihel * ndiagrams + 171];
      gpuDiagram( &graph, &graphExec, &node172, &node171, diagram172, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node173 = graphNodes[ihel * ndiagrams + 172];
      gpuDiagram( &graph, &graphExec, &node173, &node172, diagram173, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node174 = graphNodes[ihel * ndiagrams + 173];
      gpuDiagram( &graph, &graphExec, &node174, &node173, diagram174, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node175 = graphNodes[ihel * ndiagrams + 174];
      gpuDiagram( &graph, &graphExec, &node175, &node174, diagram175, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node176 = graphNodes[ihel * ndiagrams + 175];
      gpuDiagram( &graph, &graphExec, &node176, &node175, diagram176, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node177 = graphNodes[ihel * ndiagrams + 176];
      gpuDiagram( &graph, &graphExec, &node177, &node176, diagram177, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node178 = graphNodes[ihel * ndiagrams + 177];
      gpuDiagram( &graph, &graphExec, &node178, &node177, diagram178, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node179 = graphNodes[ihel * ndiagrams + 178];
      gpuDiagram( &graph, &graphExec, &node179, &node178, diagram179, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node180 = graphNodes[ihel * ndiagrams + 179];
      gpuDiagram( &graph, &graphExec, &node180, &node179, diagram180, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node181 = graphNodes[ihel * ndiagrams + 180];
      gpuDiagram( &graph, &graphExec, &node181, &node180, diagram181, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node182 = graphNodes[ihel * ndiagrams + 181];
      gpuDiagram( &graph, &graphExec, &node182, &node181, diagram182, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node183 = graphNodes[ihel * ndiagrams + 182];
      gpuDiagram( &graph, &graphExec, &node183, &node182, diagram183, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node184 = graphNodes[ihel * ndiagrams + 183];
      gpuDiagram( &graph, &graphExec, &node184, &node183, diagram184, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node185 = graphNodes[ihel * ndiagrams + 184];
      gpuDiagram( &graph, &graphExec, &node185, &node184, diagram185, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node186 = graphNodes[ihel * ndiagrams + 185];
      gpuDiagram( &graph, &graphExec, &node186, &node185, diagram186, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node187 = graphNodes[ihel * ndiagrams + 186];
      gpuDiagram( &graph, &graphExec, &node187, &node186, diagram187, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node188 = graphNodes[ihel * ndiagrams + 187];
      gpuDiagram( &graph, &graphExec, &node188, &node187, diagram188, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node189 = graphNodes[ihel * ndiagrams + 188];
      gpuDiagram( &graph, &graphExec, &node189, &node188, diagram189, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node190 = graphNodes[ihel * ndiagrams + 189];
      gpuDiagram( &graph, &graphExec, &node190, &node189, diagram190, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node191 = graphNodes[ihel * ndiagrams + 190];
      gpuDiagram( &graph, &graphExec, &node191, &node190, diagram191, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node192 = graphNodes[ihel * ndiagrams + 191];
      gpuDiagram( &graph, &graphExec, &node192, &node191, diagram192, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node193 = graphNodes[ihel * ndiagrams + 192];
      gpuDiagram( &graph, &graphExec, &node193, &node192, diagram193, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node194 = graphNodes[ihel * ndiagrams + 193];
      gpuDiagram( &graph, &graphExec, &node194, &node193, diagram194, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node195 = graphNodes[ihel * ndiagrams + 194];
      gpuDiagram( &graph, &graphExec, &node195, &node194, diagram195, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node196 = graphNodes[ihel * ndiagrams + 195];
      gpuDiagram( &graph, &graphExec, &node196, &node195, diagram196, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node197 = graphNodes[ihel * ndiagrams + 196];
      gpuDiagram( &graph, &graphExec, &node197, &node196, diagram197, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node198 = graphNodes[ihel * ndiagrams + 197];
      gpuDiagram( &graph, &graphExec, &node198, &node197, diagram198, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node199 = graphNodes[ihel * ndiagrams + 198];
      gpuDiagram( &graph, &graphExec, &node199, &node198, diagram199, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node200 = graphNodes[ihel * ndiagrams + 199];
      gpuDiagram( &graph, &graphExec, &node200, &node199, diagram200, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node201 = graphNodes[ihel * ndiagrams + 200];
      gpuDiagram( &graph, &graphExec, &node201, &node200, diagram201, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node202 = graphNodes[ihel * ndiagrams + 201];
      gpuDiagram( &graph, &graphExec, &node202, &node201, diagram202, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node203 = graphNodes[ihel * ndiagrams + 202];
      gpuDiagram( &graph, &graphExec, &node203, &node202, diagram203, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node204 = graphNodes[ihel * ndiagrams + 203];
      gpuDiagram( &graph, &graphExec, &node204, &node203, diagram204, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node205 = graphNodes[ihel * ndiagrams + 204];
      gpuDiagram( &graph, &graphExec, &node205, &node204, diagram205, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node206 = graphNodes[ihel * ndiagrams + 205];
      gpuDiagram( &graph, &graphExec, &node206, &node205, diagram206, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node207 = graphNodes[ihel * ndiagrams + 206];
      gpuDiagram( &graph, &graphExec, &node207, &node206, diagram207, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node208 = graphNodes[ihel * ndiagrams + 207];
      gpuDiagram( &graph, &graphExec, &node208, &node207, diagram208, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node209 = graphNodes[ihel * ndiagrams + 208];
      gpuDiagram( &graph, &graphExec, &node209, &node208, diagram209, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node210 = graphNodes[ihel * ndiagrams + 209];
      gpuDiagram( &graph, &graphExec, &node210, &node209, diagram210, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node211 = graphNodes[ihel * ndiagrams + 210];
      gpuDiagram( &graph, &graphExec, &node211, &node210, diagram211, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node212 = graphNodes[ihel * ndiagrams + 211];
      gpuDiagram( &graph, &graphExec, &node212, &node211, diagram212, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node213 = graphNodes[ihel * ndiagrams + 212];
      gpuDiagram( &graph, &graphExec, &node213, &node212, diagram213, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node214 = graphNodes[ihel * ndiagrams + 213];
      gpuDiagram( &graph, &graphExec, &node214, &node213, diagram214, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node215 = graphNodes[ihel * ndiagrams + 214];
      gpuDiagram( &graph, &graphExec, &node215, &node214, diagram215, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node216 = graphNodes[ihel * ndiagrams + 215];
      gpuDiagram( &graph, &graphExec, &node216, &node215, diagram216, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node217 = graphNodes[ihel * ndiagrams + 216];
      gpuDiagram( &graph, &graphExec, &node217, &node216, diagram217, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node218 = graphNodes[ihel * ndiagrams + 217];
      gpuDiagram( &graph, &graphExec, &node218, &node217, diagram218, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node219 = graphNodes[ihel * ndiagrams + 218];
      gpuDiagram( &graph, &graphExec, &node219, &node218, diagram219, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node220 = graphNodes[ihel * ndiagrams + 219];
      gpuDiagram( &graph, &graphExec, &node220, &node219, diagram220, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node221 = graphNodes[ihel * ndiagrams + 220];
      gpuDiagram( &graph, &graphExec, &node221, &node220, diagram221, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node222 = graphNodes[ihel * ndiagrams + 221];
      gpuDiagram( &graph, &graphExec, &node222, &node221, diagram222, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node223 = graphNodes[ihel * ndiagrams + 222];
      gpuDiagram( &graph, &graphExec, &node223, &node222, diagram223, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node224 = graphNodes[ihel * ndiagrams + 223];
      gpuDiagram( &graph, &graphExec, &node224, &node223, diagram224, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node225 = graphNodes[ihel * ndiagrams + 224];
      gpuDiagram( &graph, &graphExec, &node225, &node224, diagram225, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node226 = graphNodes[ihel * ndiagrams + 225];
      gpuDiagram( &graph, &graphExec, &node226, &node225, diagram226, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node227 = graphNodes[ihel * ndiagrams + 226];
      gpuDiagram( &graph, &graphExec, &node227, &node226, diagram227, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node228 = graphNodes[ihel * ndiagrams + 227];
      gpuDiagram( &graph, &graphExec, &node228, &node227, diagram228, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node229 = graphNodes[ihel * ndiagrams + 228];
      gpuDiagram( &graph, &graphExec, &node229, &node228, diagram229, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node230 = graphNodes[ihel * ndiagrams + 229];
      gpuDiagram( &graph, &graphExec, &node230, &node229, diagram230, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node231 = graphNodes[ihel * ndiagrams + 230];
      gpuDiagram( &graph, &graphExec, &node231, &node230, diagram231, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node232 = graphNodes[ihel * ndiagrams + 231];
      gpuDiagram( &graph, &graphExec, &node232, &node231, diagram232, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node233 = graphNodes[ihel * ndiagrams + 232];
      gpuDiagram( &graph, &graphExec, &node233, &node232, diagram233, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node234 = graphNodes[ihel * ndiagrams + 233];
      gpuDiagram( &graph, &graphExec, &node234, &node233, diagram234, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node235 = graphNodes[ihel * ndiagrams + 234];
      gpuDiagram( &graph, &graphExec, &node235, &node234, diagram235, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node236 = graphNodes[ihel * ndiagrams + 235];
      gpuDiagram( &graph, &graphExec, &node236, &node235, diagram236, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node237 = graphNodes[ihel * ndiagrams + 236];
      gpuDiagram( &graph, &graphExec, &node237, &node236, diagram237, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node238 = graphNodes[ihel * ndiagrams + 237];
      gpuDiagram( &graph, &graphExec, &node238, &node237, diagram238, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node239 = graphNodes[ihel * ndiagrams + 238];
      gpuDiagram( &graph, &graphExec, &node239, &node238, diagram239, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node240 = graphNodes[ihel * ndiagrams + 239];
      gpuDiagram( &graph, &graphExec, &node240, &node239, diagram240, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node241 = graphNodes[ihel * ndiagrams + 240];
      gpuDiagram( &graph, &graphExec, &node241, &node240, diagram241, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node242 = graphNodes[ihel * ndiagrams + 241];
      gpuDiagram( &graph, &graphExec, &node242, &node241, diagram242, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node243 = graphNodes[ihel * ndiagrams + 242];
      gpuDiagram( &graph, &graphExec, &node243, &node242, diagram243, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node244 = graphNodes[ihel * ndiagrams + 243];
      gpuDiagram( &graph, &graphExec, &node244, &node243, diagram244, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node245 = graphNodes[ihel * ndiagrams + 244];
      gpuDiagram( &graph, &graphExec, &node245, &node244, diagram245, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node246 = graphNodes[ihel * ndiagrams + 245];
      gpuDiagram( &graph, &graphExec, &node246, &node245, diagram246, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node247 = graphNodes[ihel * ndiagrams + 246];
      gpuDiagram( &graph, &graphExec, &node247, &node246, diagram247, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node248 = graphNodes[ihel * ndiagrams + 247];
      gpuDiagram( &graph, &graphExec, &node248, &node247, diagram248, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node249 = graphNodes[ihel * ndiagrams + 248];
      gpuDiagram( &graph, &graphExec, &node249, &node248, diagram249, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node250 = graphNodes[ihel * ndiagrams + 249];
      gpuDiagram( &graph, &graphExec, &node250, &node249, diagram250, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node251 = graphNodes[ihel * ndiagrams + 250];
      gpuDiagram( &graph, &graphExec, &node251, &node250, diagram251, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node252 = graphNodes[ihel * ndiagrams + 251];
      gpuDiagram( &graph, &graphExec, &node252, &node251, diagram252, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node253 = graphNodes[ihel * ndiagrams + 252];
      gpuDiagram( &graph, &graphExec, &node253, &node252, diagram253, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node254 = graphNodes[ihel * ndiagrams + 253];
      gpuDiagram( &graph, &graphExec, &node254, &node253, diagram254, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node255 = graphNodes[ihel * ndiagrams + 254];
      gpuDiagram( &graph, &graphExec, &node255, &node254, diagram255, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node256 = graphNodes[ihel * ndiagrams + 255];
      gpuDiagram( &graph, &graphExec, &node256, &node255, diagram256, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node257 = graphNodes[ihel * ndiagrams + 256];
      gpuDiagram( &graph, &graphExec, &node257, &node256, diagram257, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node258 = graphNodes[ihel * ndiagrams + 257];
      gpuDiagram( &graph, &graphExec, &node258, &node257, diagram258, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node259 = graphNodes[ihel * ndiagrams + 258];
      gpuDiagram( &graph, &graphExec, &node259, &node258, diagram259, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node260 = graphNodes[ihel * ndiagrams + 259];
      gpuDiagram( &graph, &graphExec, &node260, &node259, diagram260, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node261 = graphNodes[ihel * ndiagrams + 260];
      gpuDiagram( &graph, &graphExec, &node261, &node260, diagram261, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node262 = graphNodes[ihel * ndiagrams + 261];
      gpuDiagram( &graph, &graphExec, &node262, &node261, diagram262, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node263 = graphNodes[ihel * ndiagrams + 262];
      gpuDiagram( &graph, &graphExec, &node263, &node262, diagram263, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node264 = graphNodes[ihel * ndiagrams + 263];
      gpuDiagram( &graph, &graphExec, &node264, &node263, diagram264, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node265 = graphNodes[ihel * ndiagrams + 264];
      gpuDiagram( &graph, &graphExec, &node265, &node264, diagram265, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node266 = graphNodes[ihel * ndiagrams + 265];
      gpuDiagram( &graph, &graphExec, &node266, &node265, diagram266, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node267 = graphNodes[ihel * ndiagrams + 266];
      gpuDiagram( &graph, &graphExec, &node267, &node266, diagram267, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node268 = graphNodes[ihel * ndiagrams + 267];
      gpuDiagram( &graph, &graphExec, &node268, &node267, diagram268, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node269 = graphNodes[ihel * ndiagrams + 268];
      gpuDiagram( &graph, &graphExec, &node269, &node268, diagram269, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node270 = graphNodes[ihel * ndiagrams + 269];
      gpuDiagram( &graph, &graphExec, &node270, &node269, diagram270, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node271 = graphNodes[ihel * ndiagrams + 270];
      gpuDiagram( &graph, &graphExec, &node271, &node270, diagram271, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node272 = graphNodes[ihel * ndiagrams + 271];
      gpuDiagram( &graph, &graphExec, &node272, &node271, diagram272, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node273 = graphNodes[ihel * ndiagrams + 272];
      gpuDiagram( &graph, &graphExec, &node273, &node272, diagram273, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node274 = graphNodes[ihel * ndiagrams + 273];
      gpuDiagram( &graph, &graphExec, &node274, &node273, diagram274, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node275 = graphNodes[ihel * ndiagrams + 274];
      gpuDiagram( &graph, &graphExec, &node275, &node274, diagram275, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node276 = graphNodes[ihel * ndiagrams + 275];
      gpuDiagram( &graph, &graphExec, &node276, &node275, diagram276, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node277 = graphNodes[ihel * ndiagrams + 276];
      gpuDiagram( &graph, &graphExec, &node277, &node276, diagram277, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node278 = graphNodes[ihel * ndiagrams + 277];
      gpuDiagram( &graph, &graphExec, &node278, &node277, diagram278, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node279 = graphNodes[ihel * ndiagrams + 278];
      gpuDiagram( &graph, &graphExec, &node279, &node278, diagram279, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node280 = graphNodes[ihel * ndiagrams + 279];
      gpuDiagram( &graph, &graphExec, &node280, &node279, diagram280, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node281 = graphNodes[ihel * ndiagrams + 280];
      gpuDiagram( &graph, &graphExec, &node281, &node280, diagram281, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node282 = graphNodes[ihel * ndiagrams + 281];
      gpuDiagram( &graph, &graphExec, &node282, &node281, diagram282, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node283 = graphNodes[ihel * ndiagrams + 282];
      gpuDiagram( &graph, &graphExec, &node283, &node282, diagram283, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node284 = graphNodes[ihel * ndiagrams + 283];
      gpuDiagram( &graph, &graphExec, &node284, &node283, diagram284, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node285 = graphNodes[ihel * ndiagrams + 284];
      gpuDiagram( &graph, &graphExec, &node285, &node284, diagram285, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node286 = graphNodes[ihel * ndiagrams + 285];
      gpuDiagram( &graph, &graphExec, &node286, &node285, diagram286, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node287 = graphNodes[ihel * ndiagrams + 286];
      gpuDiagram( &graph, &graphExec, &node287, &node286, diagram287, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node288 = graphNodes[ihel * ndiagrams + 287];
      gpuDiagram( &graph, &graphExec, &node288, &node287, diagram288, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node289 = graphNodes[ihel * ndiagrams + 288];
      gpuDiagram( &graph, &graphExec, &node289, &node288, diagram289, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node290 = graphNodes[ihel * ndiagrams + 289];
      gpuDiagram( &graph, &graphExec, &node290, &node289, diagram290, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node291 = graphNodes[ihel * ndiagrams + 290];
      gpuDiagram( &graph, &graphExec, &node291, &node290, diagram291, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node292 = graphNodes[ihel * ndiagrams + 291];
      gpuDiagram( &graph, &graphExec, &node292, &node291, diagram292, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node293 = graphNodes[ihel * ndiagrams + 292];
      gpuDiagram( &graph, &graphExec, &node293, &node292, diagram293, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node294 = graphNodes[ihel * ndiagrams + 293];
      gpuDiagram( &graph, &graphExec, &node294, &node293, diagram294, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node295 = graphNodes[ihel * ndiagrams + 294];
      gpuDiagram( &graph, &graphExec, &node295, &node294, diagram295, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node296 = graphNodes[ihel * ndiagrams + 295];
      gpuDiagram( &graph, &graphExec, &node296, &node295, diagram296, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node297 = graphNodes[ihel * ndiagrams + 296];
      gpuDiagram( &graph, &graphExec, &node297, &node296, diagram297, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node298 = graphNodes[ihel * ndiagrams + 297];
      gpuDiagram( &graph, &graphExec, &node298, &node297, diagram298, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node299 = graphNodes[ihel * ndiagrams + 298];
      gpuDiagram( &graph, &graphExec, &node299, &node298, diagram299, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node300 = graphNodes[ihel * ndiagrams + 299];
      gpuDiagram( &graph, &graphExec, &node300, &node299, diagram300, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node301 = graphNodes[ihel * ndiagrams + 300];
      gpuDiagram( &graph, &graphExec, &node301, &node300, diagram301, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node302 = graphNodes[ihel * ndiagrams + 301];
      gpuDiagram( &graph, &graphExec, &node302, &node301, diagram302, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node303 = graphNodes[ihel * ndiagrams + 302];
      gpuDiagram( &graph, &graphExec, &node303, &node302, diagram303, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node304 = graphNodes[ihel * ndiagrams + 303];
      gpuDiagram( &graph, &graphExec, &node304, &node303, diagram304, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node305 = graphNodes[ihel * ndiagrams + 304];
      gpuDiagram( &graph, &graphExec, &node305, &node304, diagram305, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node306 = graphNodes[ihel * ndiagrams + 305];
      gpuDiagram( &graph, &graphExec, &node306, &node305, diagram306, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node307 = graphNodes[ihel * ndiagrams + 306];
      gpuDiagram( &graph, &graphExec, &node307, &node306, diagram307, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node308 = graphNodes[ihel * ndiagrams + 307];
      gpuDiagram( &graph, &graphExec, &node308, &node307, diagram308, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node309 = graphNodes[ihel * ndiagrams + 308];
      gpuDiagram( &graph, &graphExec, &node309, &node308, diagram309, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node310 = graphNodes[ihel * ndiagrams + 309];
      gpuDiagram( &graph, &graphExec, &node310, &node309, diagram310, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node311 = graphNodes[ihel * ndiagrams + 310];
      gpuDiagram( &graph, &graphExec, &node311, &node310, diagram311, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node312 = graphNodes[ihel * ndiagrams + 311];
      gpuDiagram( &graph, &graphExec, &node312, &node311, diagram312, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node313 = graphNodes[ihel * ndiagrams + 312];
      gpuDiagram( &graph, &graphExec, &node313, &node312, diagram313, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node314 = graphNodes[ihel * ndiagrams + 313];
      gpuDiagram( &graph, &graphExec, &node314, &node313, diagram314, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node315 = graphNodes[ihel * ndiagrams + 314];
      gpuDiagram( &graph, &graphExec, &node315, &node314, diagram315, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node316 = graphNodes[ihel * ndiagrams + 315];
      gpuDiagram( &graph, &graphExec, &node316, &node315, diagram316, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node317 = graphNodes[ihel * ndiagrams + 316];
      gpuDiagram( &graph, &graphExec, &node317, &node316, diagram317, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node318 = graphNodes[ihel * ndiagrams + 317];
      gpuDiagram( &graph, &graphExec, &node318, &node317, diagram318, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node319 = graphNodes[ihel * ndiagrams + 318];
      gpuDiagram( &graph, &graphExec, &node319, &node318, diagram319, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node320 = graphNodes[ihel * ndiagrams + 319];
      gpuDiagram( &graph, &graphExec, &node320, &node319, diagram320, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node321 = graphNodes[ihel * ndiagrams + 320];
      gpuDiagram( &graph, &graphExec, &node321, &node320, diagram321, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node322 = graphNodes[ihel * ndiagrams + 321];
      gpuDiagram( &graph, &graphExec, &node322, &node321, diagram322, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node323 = graphNodes[ihel * ndiagrams + 322];
      gpuDiagram( &graph, &graphExec, &node323, &node322, diagram323, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node324 = graphNodes[ihel * ndiagrams + 323];
      gpuDiagram( &graph, &graphExec, &node324, &node323, diagram324, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node325 = graphNodes[ihel * ndiagrams + 324];
      gpuDiagram( &graph, &graphExec, &node325, &node324, diagram325, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node326 = graphNodes[ihel * ndiagrams + 325];
      gpuDiagram( &graph, &graphExec, &node326, &node325, diagram326, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node327 = graphNodes[ihel * ndiagrams + 326];
      gpuDiagram( &graph, &graphExec, &node327, &node326, diagram327, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node328 = graphNodes[ihel * ndiagrams + 327];
      gpuDiagram( &graph, &graphExec, &node328, &node327, diagram328, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node329 = graphNodes[ihel * ndiagrams + 328];
      gpuDiagram( &graph, &graphExec, &node329, &node328, diagram329, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node330 = graphNodes[ihel * ndiagrams + 329];
      gpuDiagram( &graph, &graphExec, &node330, &node329, diagram330, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node331 = graphNodes[ihel * ndiagrams + 330];
      gpuDiagram( &graph, &graphExec, &node331, &node330, diagram331, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node332 = graphNodes[ihel * ndiagrams + 331];
      gpuDiagram( &graph, &graphExec, &node332, &node331, diagram332, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node333 = graphNodes[ihel * ndiagrams + 332];
      gpuDiagram( &graph, &graphExec, &node333, &node332, diagram333, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node334 = graphNodes[ihel * ndiagrams + 333];
      gpuDiagram( &graph, &graphExec, &node334, &node333, diagram334, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node335 = graphNodes[ihel * ndiagrams + 334];
      gpuDiagram( &graph, &graphExec, &node335, &node334, diagram335, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node336 = graphNodes[ihel * ndiagrams + 335];
      gpuDiagram( &graph, &graphExec, &node336, &node335, diagram336, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node337 = graphNodes[ihel * ndiagrams + 336];
      gpuDiagram( &graph, &graphExec, &node337, &node336, diagram337, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node338 = graphNodes[ihel * ndiagrams + 337];
      gpuDiagram( &graph, &graphExec, &node338, &node337, diagram338, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node339 = graphNodes[ihel * ndiagrams + 338];
      gpuDiagram( &graph, &graphExec, &node339, &node338, diagram339, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node340 = graphNodes[ihel * ndiagrams + 339];
      gpuDiagram( &graph, &graphExec, &node340, &node339, diagram340, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node341 = graphNodes[ihel * ndiagrams + 340];
      gpuDiagram( &graph, &graphExec, &node341, &node340, diagram341, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node342 = graphNodes[ihel * ndiagrams + 341];
      gpuDiagram( &graph, &graphExec, &node342, &node341, diagram342, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node343 = graphNodes[ihel * ndiagrams + 342];
      gpuDiagram( &graph, &graphExec, &node343, &node342, diagram343, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node344 = graphNodes[ihel * ndiagrams + 343];
      gpuDiagram( &graph, &graphExec, &node344, &node343, diagram344, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node345 = graphNodes[ihel * ndiagrams + 344];
      gpuDiagram( &graph, &graphExec, &node345, &node344, diagram345, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node346 = graphNodes[ihel * ndiagrams + 345];
      gpuDiagram( &graph, &graphExec, &node346, &node345, diagram346, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node347 = graphNodes[ihel * ndiagrams + 346];
      gpuDiagram( &graph, &graphExec, &node347, &node346, diagram347, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node348 = graphNodes[ihel * ndiagrams + 347];
      gpuDiagram( &graph, &graphExec, &node348, &node347, diagram348, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node349 = graphNodes[ihel * ndiagrams + 348];
      gpuDiagram( &graph, &graphExec, &node349, &node348, diagram349, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node350 = graphNodes[ihel * ndiagrams + 349];
      gpuDiagram( &graph, &graphExec, &node350, &node349, diagram350, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node351 = graphNodes[ihel * ndiagrams + 350];
      gpuDiagram( &graph, &graphExec, &node351, &node350, diagram351, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node352 = graphNodes[ihel * ndiagrams + 351];
      gpuDiagram( &graph, &graphExec, &node352, &node351, diagram352, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node353 = graphNodes[ihel * ndiagrams + 352];
      gpuDiagram( &graph, &graphExec, &node353, &node352, diagram353, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node354 = graphNodes[ihel * ndiagrams + 353];
      gpuDiagram( &graph, &graphExec, &node354, &node353, diagram354, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node355 = graphNodes[ihel * ndiagrams + 354];
      gpuDiagram( &graph, &graphExec, &node355, &node354, diagram355, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node356 = graphNodes[ihel * ndiagrams + 355];
      gpuDiagram( &graph, &graphExec, &node356, &node355, diagram356, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node357 = graphNodes[ihel * ndiagrams + 356];
      gpuDiagram( &graph, &graphExec, &node357, &node356, diagram357, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node358 = graphNodes[ihel * ndiagrams + 357];
      gpuDiagram( &graph, &graphExec, &node358, &node357, diagram358, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node359 = graphNodes[ihel * ndiagrams + 358];
      gpuDiagram( &graph, &graphExec, &node359, &node358, diagram359, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node360 = graphNodes[ihel * ndiagrams + 359];
      gpuDiagram( &graph, &graphExec, &node360, &node359, diagram360, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node361 = graphNodes[ihel * ndiagrams + 360];
      gpuDiagram( &graph, &graphExec, &node361, &node360, diagram361, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node362 = graphNodes[ihel * ndiagrams + 361];
      gpuDiagram( &graph, &graphExec, &node362, &node361, diagram362, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node363 = graphNodes[ihel * ndiagrams + 362];
      gpuDiagram( &graph, &graphExec, &node363, &node362, diagram363, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node364 = graphNodes[ihel * ndiagrams + 363];
      gpuDiagram( &graph, &graphExec, &node364, &node363, diagram364, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node365 = graphNodes[ihel * ndiagrams + 364];
      gpuDiagram( &graph, &graphExec, &node365, &node364, diagram365, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node366 = graphNodes[ihel * ndiagrams + 365];
      gpuDiagram( &graph, &graphExec, &node366, &node365, diagram366, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node367 = graphNodes[ihel * ndiagrams + 366];
      gpuDiagram( &graph, &graphExec, &node367, &node366, diagram367, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node368 = graphNodes[ihel * ndiagrams + 367];
      gpuDiagram( &graph, &graphExec, &node368, &node367, diagram368, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node369 = graphNodes[ihel * ndiagrams + 368];
      gpuDiagram( &graph, &graphExec, &node369, &node368, diagram369, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node370 = graphNodes[ihel * ndiagrams + 369];
      gpuDiagram( &graph, &graphExec, &node370, &node369, diagram370, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node371 = graphNodes[ihel * ndiagrams + 370];
      gpuDiagram( &graph, &graphExec, &node371, &node370, diagram371, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node372 = graphNodes[ihel * ndiagrams + 371];
      gpuDiagram( &graph, &graphExec, &node372, &node371, diagram372, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node373 = graphNodes[ihel * ndiagrams + 372];
      gpuDiagram( &graph, &graphExec, &node373, &node372, diagram373, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node374 = graphNodes[ihel * ndiagrams + 373];
      gpuDiagram( &graph, &graphExec, &node374, &node373, diagram374, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node375 = graphNodes[ihel * ndiagrams + 374];
      gpuDiagram( &graph, &graphExec, &node375, &node374, diagram375, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node376 = graphNodes[ihel * ndiagrams + 375];
      gpuDiagram( &graph, &graphExec, &node376, &node375, diagram376, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node377 = graphNodes[ihel * ndiagrams + 376];
      gpuDiagram( &graph, &graphExec, &node377, &node376, diagram377, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node378 = graphNodes[ihel * ndiagrams + 377];
      gpuDiagram( &graph, &graphExec, &node378, &node377, diagram378, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node379 = graphNodes[ihel * ndiagrams + 378];
      gpuDiagram( &graph, &graphExec, &node379, &node378, diagram379, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node380 = graphNodes[ihel * ndiagrams + 379];
      gpuDiagram( &graph, &graphExec, &node380, &node379, diagram380, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node381 = graphNodes[ihel * ndiagrams + 380];
      gpuDiagram( &graph, &graphExec, &node381, &node380, diagram381, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node382 = graphNodes[ihel * ndiagrams + 381];
      gpuDiagram( &graph, &graphExec, &node382, &node381, diagram382, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node383 = graphNodes[ihel * ndiagrams + 382];
      gpuDiagram( &graph, &graphExec, &node383, &node382, diagram383, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node384 = graphNodes[ihel * ndiagrams + 383];
      gpuDiagram( &graph, &graphExec, &node384, &node383, diagram384, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node385 = graphNodes[ihel * ndiagrams + 384];
      gpuDiagram( &graph, &graphExec, &node385, &node384, diagram385, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node386 = graphNodes[ihel * ndiagrams + 385];
      gpuDiagram( &graph, &graphExec, &node386, &node385, diagram386, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node387 = graphNodes[ihel * ndiagrams + 386];
      gpuDiagram( &graph, &graphExec, &node387, &node386, diagram387, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node388 = graphNodes[ihel * ndiagrams + 387];
      gpuDiagram( &graph, &graphExec, &node388, &node387, diagram388, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node389 = graphNodes[ihel * ndiagrams + 388];
      gpuDiagram( &graph, &graphExec, &node389, &node388, diagram389, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node390 = graphNodes[ihel * ndiagrams + 389];
      gpuDiagram( &graph, &graphExec, &node390, &node389, diagram390, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node391 = graphNodes[ihel * ndiagrams + 390];
      gpuDiagram( &graph, &graphExec, &node391, &node390, diagram391, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node392 = graphNodes[ihel * ndiagrams + 391];
      gpuDiagram( &graph, &graphExec, &node392, &node391, diagram392, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node393 = graphNodes[ihel * ndiagrams + 392];
      gpuDiagram( &graph, &graphExec, &node393, &node392, diagram393, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node394 = graphNodes[ihel * ndiagrams + 393];
      gpuDiagram( &graph, &graphExec, &node394, &node393, diagram394, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node395 = graphNodes[ihel * ndiagrams + 394];
      gpuDiagram( &graph, &graphExec, &node395, &node394, diagram395, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node396 = graphNodes[ihel * ndiagrams + 395];
      gpuDiagram( &graph, &graphExec, &node396, &node395, diagram396, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node397 = graphNodes[ihel * ndiagrams + 396];
      gpuDiagram( &graph, &graphExec, &node397, &node396, diagram397, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node398 = graphNodes[ihel * ndiagrams + 397];
      gpuDiagram( &graph, &graphExec, &node398, &node397, diagram398, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node399 = graphNodes[ihel * ndiagrams + 398];
      gpuDiagram( &graph, &graphExec, &node399, &node398, diagram399, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node400 = graphNodes[ihel * ndiagrams + 399];
      gpuDiagram( &graph, &graphExec, &node400, &node399, diagram400, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node401 = graphNodes[ihel * ndiagrams + 400];
      gpuDiagram( &graph, &graphExec, &node401, &node400, diagram401, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node402 = graphNodes[ihel * ndiagrams + 401];
      gpuDiagram( &graph, &graphExec, &node402, &node401, diagram402, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node403 = graphNodes[ihel * ndiagrams + 402];
      gpuDiagram( &graph, &graphExec, &node403, &node402, diagram403, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node404 = graphNodes[ihel * ndiagrams + 403];
      gpuDiagram( &graph, &graphExec, &node404, &node403, diagram404, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node405 = graphNodes[ihel * ndiagrams + 404];
      gpuDiagram( &graph, &graphExec, &node405, &node404, diagram405, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node406 = graphNodes[ihel * ndiagrams + 405];
      gpuDiagram( &graph, &graphExec, &node406, &node405, diagram406, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node407 = graphNodes[ihel * ndiagrams + 406];
      gpuDiagram( &graph, &graphExec, &node407, &node406, diagram407, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node408 = graphNodes[ihel * ndiagrams + 407];
      gpuDiagram( &graph, &graphExec, &node408, &node407, diagram408, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node409 = graphNodes[ihel * ndiagrams + 408];
      gpuDiagram( &graph, &graphExec, &node409, &node408, diagram409, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node410 = graphNodes[ihel * ndiagrams + 409];
      gpuDiagram( &graph, &graphExec, &node410, &node409, diagram410, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node411 = graphNodes[ihel * ndiagrams + 410];
      gpuDiagram( &graph, &graphExec, &node411, &node410, diagram411, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node412 = graphNodes[ihel * ndiagrams + 411];
      gpuDiagram( &graph, &graphExec, &node412, &node411, diagram412, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node413 = graphNodes[ihel * ndiagrams + 412];
      gpuDiagram( &graph, &graphExec, &node413, &node412, diagram413, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node414 = graphNodes[ihel * ndiagrams + 413];
      gpuDiagram( &graph, &graphExec, &node414, &node413, diagram414, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node415 = graphNodes[ihel * ndiagrams + 414];
      gpuDiagram( &graph, &graphExec, &node415, &node414, diagram415, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node416 = graphNodes[ihel * ndiagrams + 415];
      gpuDiagram( &graph, &graphExec, &node416, &node415, diagram416, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node417 = graphNodes[ihel * ndiagrams + 416];
      gpuDiagram( &graph, &graphExec, &node417, &node416, diagram417, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node418 = graphNodes[ihel * ndiagrams + 417];
      gpuDiagram( &graph, &graphExec, &node418, &node417, diagram418, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node419 = graphNodes[ihel * ndiagrams + 418];
      gpuDiagram( &graph, &graphExec, &node419, &node418, diagram419, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node420 = graphNodes[ihel * ndiagrams + 419];
      gpuDiagram( &graph, &graphExec, &node420, &node419, diagram420, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node421 = graphNodes[ihel * ndiagrams + 420];
      gpuDiagram( &graph, &graphExec, &node421, &node420, diagram421, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node422 = graphNodes[ihel * ndiagrams + 421];
      gpuDiagram( &graph, &graphExec, &node422, &node421, diagram422, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node423 = graphNodes[ihel * ndiagrams + 422];
      gpuDiagram( &graph, &graphExec, &node423, &node422, diagram423, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node424 = graphNodes[ihel * ndiagrams + 423];
      gpuDiagram( &graph, &graphExec, &node424, &node423, diagram424, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node425 = graphNodes[ihel * ndiagrams + 424];
      gpuDiagram( &graph, &graphExec, &node425, &node424, diagram425, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node426 = graphNodes[ihel * ndiagrams + 425];
      gpuDiagram( &graph, &graphExec, &node426, &node425, diagram426, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node427 = graphNodes[ihel * ndiagrams + 426];
      gpuDiagram( &graph, &graphExec, &node427, &node426, diagram427, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node428 = graphNodes[ihel * ndiagrams + 427];
      gpuDiagram( &graph, &graphExec, &node428, &node427, diagram428, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node429 = graphNodes[ihel * ndiagrams + 428];
      gpuDiagram( &graph, &graphExec, &node429, &node428, diagram429, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node430 = graphNodes[ihel * ndiagrams + 429];
      gpuDiagram( &graph, &graphExec, &node430, &node429, diagram430, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node431 = graphNodes[ihel * ndiagrams + 430];
      gpuDiagram( &graph, &graphExec, &node431, &node430, diagram431, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node432 = graphNodes[ihel * ndiagrams + 431];
      gpuDiagram( &graph, &graphExec, &node432, &node431, diagram432, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node433 = graphNodes[ihel * ndiagrams + 432];
      gpuDiagram( &graph, &graphExec, &node433, &node432, diagram433, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node434 = graphNodes[ihel * ndiagrams + 433];
      gpuDiagram( &graph, &graphExec, &node434, &node433, diagram434, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node435 = graphNodes[ihel * ndiagrams + 434];
      gpuDiagram( &graph, &graphExec, &node435, &node434, diagram435, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node436 = graphNodes[ihel * ndiagrams + 435];
      gpuDiagram( &graph, &graphExec, &node436, &node435, diagram436, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node437 = graphNodes[ihel * ndiagrams + 436];
      gpuDiagram( &graph, &graphExec, &node437, &node436, diagram437, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node438 = graphNodes[ihel * ndiagrams + 437];
      gpuDiagram( &graph, &graphExec, &node438, &node437, diagram438, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node439 = graphNodes[ihel * ndiagrams + 438];
      gpuDiagram( &graph, &graphExec, &node439, &node438, diagram439, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node440 = graphNodes[ihel * ndiagrams + 439];
      gpuDiagram( &graph, &graphExec, &node440, &node439, diagram440, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node441 = graphNodes[ihel * ndiagrams + 440];
      gpuDiagram( &graph, &graphExec, &node441, &node440, diagram441, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node442 = graphNodes[ihel * ndiagrams + 441];
      gpuDiagram( &graph, &graphExec, &node442, &node441, diagram442, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node443 = graphNodes[ihel * ndiagrams + 442];
      gpuDiagram( &graph, &graphExec, &node443, &node442, diagram443, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node444 = graphNodes[ihel * ndiagrams + 443];
      gpuDiagram( &graph, &graphExec, &node444, &node443, diagram444, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node445 = graphNodes[ihel * ndiagrams + 444];
      gpuDiagram( &graph, &graphExec, &node445, &node444, diagram445, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node446 = graphNodes[ihel * ndiagrams + 445];
      gpuDiagram( &graph, &graphExec, &node446, &node445, diagram446, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node447 = graphNodes[ihel * ndiagrams + 446];
      gpuDiagram( &graph, &graphExec, &node447, &node446, diagram447, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node448 = graphNodes[ihel * ndiagrams + 447];
      gpuDiagram( &graph, &graphExec, &node448, &node447, diagram448, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node449 = graphNodes[ihel * ndiagrams + 448];
      gpuDiagram( &graph, &graphExec, &node449, &node448, diagram449, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node450 = graphNodes[ihel * ndiagrams + 449];
      gpuDiagram( &graph, &graphExec, &node450, &node449, diagram450, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node451 = graphNodes[ihel * ndiagrams + 450];
      gpuDiagram( &graph, &graphExec, &node451, &node450, diagram451, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node452 = graphNodes[ihel * ndiagrams + 451];
      gpuDiagram( &graph, &graphExec, &node452, &node451, diagram452, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node453 = graphNodes[ihel * ndiagrams + 452];
      gpuDiagram( &graph, &graphExec, &node453, &node452, diagram453, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node454 = graphNodes[ihel * ndiagrams + 453];
      gpuDiagram( &graph, &graphExec, &node454, &node453, diagram454, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node455 = graphNodes[ihel * ndiagrams + 454];
      gpuDiagram( &graph, &graphExec, &node455, &node454, diagram455, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node456 = graphNodes[ihel * ndiagrams + 455];
      gpuDiagram( &graph, &graphExec, &node456, &node455, diagram456, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node457 = graphNodes[ihel * ndiagrams + 456];
      gpuDiagram( &graph, &graphExec, &node457, &node456, diagram457, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node458 = graphNodes[ihel * ndiagrams + 457];
      gpuDiagram( &graph, &graphExec, &node458, &node457, diagram458, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node459 = graphNodes[ihel * ndiagrams + 458];
      gpuDiagram( &graph, &graphExec, &node459, &node458, diagram459, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node460 = graphNodes[ihel * ndiagrams + 459];
      gpuDiagram( &graph, &graphExec, &node460, &node459, diagram460, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node461 = graphNodes[ihel * ndiagrams + 460];
      gpuDiagram( &graph, &graphExec, &node461, &node460, diagram461, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node462 = graphNodes[ihel * ndiagrams + 461];
      gpuDiagram( &graph, &graphExec, &node462, &node461, diagram462, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node463 = graphNodes[ihel * ndiagrams + 462];
      gpuDiagram( &graph, &graphExec, &node463, &node462, diagram463, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node464 = graphNodes[ihel * ndiagrams + 463];
      gpuDiagram( &graph, &graphExec, &node464, &node463, diagram464, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node465 = graphNodes[ihel * ndiagrams + 464];
      gpuDiagram( &graph, &graphExec, &node465, &node464, diagram465, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node466 = graphNodes[ihel * ndiagrams + 465];
      gpuDiagram( &graph, &graphExec, &node466, &node465, diagram466, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node467 = graphNodes[ihel * ndiagrams + 466];
      gpuDiagram( &graph, &graphExec, &node467, &node466, diagram467, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node468 = graphNodes[ihel * ndiagrams + 467];
      gpuDiagram( &graph, &graphExec, &node468, &node467, diagram468, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node469 = graphNodes[ihel * ndiagrams + 468];
      gpuDiagram( &graph, &graphExec, &node469, &node468, diagram469, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node470 = graphNodes[ihel * ndiagrams + 469];
      gpuDiagram( &graph, &graphExec, &node470, &node469, diagram470, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node471 = graphNodes[ihel * ndiagrams + 470];
      gpuDiagram( &graph, &graphExec, &node471, &node470, diagram471, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node472 = graphNodes[ihel * ndiagrams + 471];
      gpuDiagram( &graph, &graphExec, &node472, &node471, diagram472, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node473 = graphNodes[ihel * ndiagrams + 472];
      gpuDiagram( &graph, &graphExec, &node473, &node472, diagram473, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node474 = graphNodes[ihel * ndiagrams + 473];
      gpuDiagram( &graph, &graphExec, &node474, &node473, diagram474, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node475 = graphNodes[ihel * ndiagrams + 474];
      gpuDiagram( &graph, &graphExec, &node475, &node474, diagram475, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node476 = graphNodes[ihel * ndiagrams + 475];
      gpuDiagram( &graph, &graphExec, &node476, &node475, diagram476, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node477 = graphNodes[ihel * ndiagrams + 476];
      gpuDiagram( &graph, &graphExec, &node477, &node476, diagram477, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node478 = graphNodes[ihel * ndiagrams + 477];
      gpuDiagram( &graph, &graphExec, &node478, &node477, diagram478, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node479 = graphNodes[ihel * ndiagrams + 478];
      gpuDiagram( &graph, &graphExec, &node479, &node478, diagram479, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node480 = graphNodes[ihel * ndiagrams + 479];
      gpuDiagram( &graph, &graphExec, &node480, &node479, diagram480, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node481 = graphNodes[ihel * ndiagrams + 480];
      gpuDiagram( &graph, &graphExec, &node481, &node480, diagram481, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node482 = graphNodes[ihel * ndiagrams + 481];
      gpuDiagram( &graph, &graphExec, &node482, &node481, diagram482, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node483 = graphNodes[ihel * ndiagrams + 482];
      gpuDiagram( &graph, &graphExec, &node483, &node482, diagram483, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node484 = graphNodes[ihel * ndiagrams + 483];
      gpuDiagram( &graph, &graphExec, &node484, &node483, diagram484, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node485 = graphNodes[ihel * ndiagrams + 484];
      gpuDiagram( &graph, &graphExec, &node485, &node484, diagram485, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node486 = graphNodes[ihel * ndiagrams + 485];
      gpuDiagram( &graph, &graphExec, &node486, &node485, diagram486, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node487 = graphNodes[ihel * ndiagrams + 486];
      gpuDiagram( &graph, &graphExec, &node487, &node486, diagram487, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node488 = graphNodes[ihel * ndiagrams + 487];
      gpuDiagram( &graph, &graphExec, &node488, &node487, diagram488, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node489 = graphNodes[ihel * ndiagrams + 488];
      gpuDiagram( &graph, &graphExec, &node489, &node488, diagram489, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node490 = graphNodes[ihel * ndiagrams + 489];
      gpuDiagram( &graph, &graphExec, &node490, &node489, diagram490, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node491 = graphNodes[ihel * ndiagrams + 490];
      gpuDiagram( &graph, &graphExec, &node491, &node490, diagram491, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node492 = graphNodes[ihel * ndiagrams + 491];
      gpuDiagram( &graph, &graphExec, &node492, &node491, diagram492, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node493 = graphNodes[ihel * ndiagrams + 492];
      gpuDiagram( &graph, &graphExec, &node493, &node492, diagram493, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node494 = graphNodes[ihel * ndiagrams + 493];
      gpuDiagram( &graph, &graphExec, &node494, &node493, diagram494, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node495 = graphNodes[ihel * ndiagrams + 494];
      gpuDiagram( &graph, &graphExec, &node495, &node494, diagram495, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node496 = graphNodes[ihel * ndiagrams + 495];
      gpuDiagram( &graph, &graphExec, &node496, &node495, diagram496, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node497 = graphNodes[ihel * ndiagrams + 496];
      gpuDiagram( &graph, &graphExec, &node497, &node496, diagram497, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node498 = graphNodes[ihel * ndiagrams + 497];
      gpuDiagram( &graph, &graphExec, &node498, &node497, diagram498, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node499 = graphNodes[ihel * ndiagrams + 498];
      gpuDiagram( &graph, &graphExec, &node499, &node498, diagram499, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node500 = graphNodes[ihel * ndiagrams + 499];
      gpuDiagram( &graph, &graphExec, &node500, &node499, diagram500, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node501 = graphNodes[ihel * ndiagrams + 500];
      gpuDiagram( &graph, &graphExec, &node501, &node500, diagram501, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node502 = graphNodes[ihel * ndiagrams + 501];
      gpuDiagram( &graph, &graphExec, &node502, &node501, diagram502, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node503 = graphNodes[ihel * ndiagrams + 502];
      gpuDiagram( &graph, &graphExec, &node503, &node502, diagram503, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node504 = graphNodes[ihel * ndiagrams + 503];
      gpuDiagram( &graph, &graphExec, &node504, &node503, diagram504, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node505 = graphNodes[ihel * ndiagrams + 504];
      gpuDiagram( &graph, &graphExec, &node505, &node504, diagram505, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node506 = graphNodes[ihel * ndiagrams + 505];
      gpuDiagram( &graph, &graphExec, &node506, &node505, diagram506, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node507 = graphNodes[ihel * ndiagrams + 506];
      gpuDiagram( &graph, &graphExec, &node507, &node506, diagram507, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node508 = graphNodes[ihel * ndiagrams + 507];
      gpuDiagram( &graph, &graphExec, &node508, &node507, diagram508, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node509 = graphNodes[ihel * ndiagrams + 508];
      gpuDiagram( &graph, &graphExec, &node509, &node508, diagram509, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node510 = graphNodes[ihel * ndiagrams + 509];
      gpuDiagram( &graph, &graphExec, &node510, &node509, diagram510, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node511 = graphNodes[ihel * ndiagrams + 510];
      gpuDiagram( &graph, &graphExec, &node511, &node510, diagram511, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node512 = graphNodes[ihel * ndiagrams + 511];
      gpuDiagram( &graph, &graphExec, &node512, &node511, diagram512, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node513 = graphNodes[ihel * ndiagrams + 512];
      gpuDiagram( &graph, &graphExec, &node513, &node512, diagram513, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node514 = graphNodes[ihel * ndiagrams + 513];
      gpuDiagram( &graph, &graphExec, &node514, &node513, diagram514, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node515 = graphNodes[ihel * ndiagrams + 514];
      gpuDiagram( &graph, &graphExec, &node515, &node514, diagram515, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node516 = graphNodes[ihel * ndiagrams + 515];
      gpuDiagram( &graph, &graphExec, &node516, &node515, diagram516, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node517 = graphNodes[ihel * ndiagrams + 516];
      gpuDiagram( &graph, &graphExec, &node517, &node516, diagram517, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node518 = graphNodes[ihel * ndiagrams + 517];
      gpuDiagram( &graph, &graphExec, &node518, &node517, diagram518, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node519 = graphNodes[ihel * ndiagrams + 518];
      gpuDiagram( &graph, &graphExec, &node519, &node518, diagram519, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node520 = graphNodes[ihel * ndiagrams + 519];
      gpuDiagram( &graph, &graphExec, &node520, &node519, diagram520, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node521 = graphNodes[ihel * ndiagrams + 520];
      gpuDiagram( &graph, &graphExec, &node521, &node520, diagram521, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node522 = graphNodes[ihel * ndiagrams + 521];
      gpuDiagram( &graph, &graphExec, &node522, &node521, diagram522, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node523 = graphNodes[ihel * ndiagrams + 522];
      gpuDiagram( &graph, &graphExec, &node523, &node522, diagram523, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node524 = graphNodes[ihel * ndiagrams + 523];
      gpuDiagram( &graph, &graphExec, &node524, &node523, diagram524, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node525 = graphNodes[ihel * ndiagrams + 524];
      gpuDiagram( &graph, &graphExec, &node525, &node524, diagram525, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node526 = graphNodes[ihel * ndiagrams + 525];
      gpuDiagram( &graph, &graphExec, &node526, &node525, diagram526, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node527 = graphNodes[ihel * ndiagrams + 526];
      gpuDiagram( &graph, &graphExec, &node527, &node526, diagram527, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node528 = graphNodes[ihel * ndiagrams + 527];
      gpuDiagram( &graph, &graphExec, &node528, &node527, diagram528, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node529 = graphNodes[ihel * ndiagrams + 528];
      gpuDiagram( &graph, &graphExec, &node529, &node528, diagram529, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node530 = graphNodes[ihel * ndiagrams + 529];
      gpuDiagram( &graph, &graphExec, &node530, &node529, diagram530, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node531 = graphNodes[ihel * ndiagrams + 530];
      gpuDiagram( &graph, &graphExec, &node531, &node530, diagram531, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node532 = graphNodes[ihel * ndiagrams + 531];
      gpuDiagram( &graph, &graphExec, &node532, &node531, diagram532, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node533 = graphNodes[ihel * ndiagrams + 532];
      gpuDiagram( &graph, &graphExec, &node533, &node532, diagram533, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node534 = graphNodes[ihel * ndiagrams + 533];
      gpuDiagram( &graph, &graphExec, &node534, &node533, diagram534, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node535 = graphNodes[ihel * ndiagrams + 534];
      gpuDiagram( &graph, &graphExec, &node535, &node534, diagram535, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node536 = graphNodes[ihel * ndiagrams + 535];
      gpuDiagram( &graph, &graphExec, &node536, &node535, diagram536, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node537 = graphNodes[ihel * ndiagrams + 536];
      gpuDiagram( &graph, &graphExec, &node537, &node536, diagram537, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node538 = graphNodes[ihel * ndiagrams + 537];
      gpuDiagram( &graph, &graphExec, &node538, &node537, diagram538, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node539 = graphNodes[ihel * ndiagrams + 538];
      gpuDiagram( &graph, &graphExec, &node539, &node538, diagram539, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node540 = graphNodes[ihel * ndiagrams + 539];
      gpuDiagram( &graph, &graphExec, &node540, &node539, diagram540, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node541 = graphNodes[ihel * ndiagrams + 540];
      gpuDiagram( &graph, &graphExec, &node541, &node540, diagram541, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node542 = graphNodes[ihel * ndiagrams + 541];
      gpuDiagram( &graph, &graphExec, &node542, &node541, diagram542, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node543 = graphNodes[ihel * ndiagrams + 542];
      gpuDiagram( &graph, &graphExec, &node543, &node542, diagram543, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node544 = graphNodes[ihel * ndiagrams + 543];
      gpuDiagram( &graph, &graphExec, &node544, &node543, diagram544, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node545 = graphNodes[ihel * ndiagrams + 544];
      gpuDiagram( &graph, &graphExec, &node545, &node544, diagram545, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node546 = graphNodes[ihel * ndiagrams + 545];
      gpuDiagram( &graph, &graphExec, &node546, &node545, diagram546, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node547 = graphNodes[ihel * ndiagrams + 546];
      gpuDiagram( &graph, &graphExec, &node547, &node546, diagram547, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node548 = graphNodes[ihel * ndiagrams + 547];
      gpuDiagram( &graph, &graphExec, &node548, &node547, diagram548, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node549 = graphNodes[ihel * ndiagrams + 548];
      gpuDiagram( &graph, &graphExec, &node549, &node548, diagram549, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node550 = graphNodes[ihel * ndiagrams + 549];
      gpuDiagram( &graph, &graphExec, &node550, &node549, diagram550, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node551 = graphNodes[ihel * ndiagrams + 550];
      gpuDiagram( &graph, &graphExec, &node551, &node550, diagram551, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node552 = graphNodes[ihel * ndiagrams + 551];
      gpuDiagram( &graph, &graphExec, &node552, &node551, diagram552, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node553 = graphNodes[ihel * ndiagrams + 552];
      gpuDiagram( &graph, &graphExec, &node553, &node552, diagram553, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node554 = graphNodes[ihel * ndiagrams + 553];
      gpuDiagram( &graph, &graphExec, &node554, &node553, diagram554, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node555 = graphNodes[ihel * ndiagrams + 554];
      gpuDiagram( &graph, &graphExec, &node555, &node554, diagram555, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node556 = graphNodes[ihel * ndiagrams + 555];
      gpuDiagram( &graph, &graphExec, &node556, &node555, diagram556, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node557 = graphNodes[ihel * ndiagrams + 556];
      gpuDiagram( &graph, &graphExec, &node557, &node556, diagram557, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node558 = graphNodes[ihel * ndiagrams + 557];
      gpuDiagram( &graph, &graphExec, &node558, &node557, diagram558, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node559 = graphNodes[ihel * ndiagrams + 558];
      gpuDiagram( &graph, &graphExec, &node559, &node558, diagram559, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node560 = graphNodes[ihel * ndiagrams + 559];
      gpuDiagram( &graph, &graphExec, &node560, &node559, diagram560, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node561 = graphNodes[ihel * ndiagrams + 560];
      gpuDiagram( &graph, &graphExec, &node561, &node560, diagram561, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node562 = graphNodes[ihel * ndiagrams + 561];
      gpuDiagram( &graph, &graphExec, &node562, &node561, diagram562, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node563 = graphNodes[ihel * ndiagrams + 562];
      gpuDiagram( &graph, &graphExec, &node563, &node562, diagram563, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node564 = graphNodes[ihel * ndiagrams + 563];
      gpuDiagram( &graph, &graphExec, &node564, &node563, diagram564, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node565 = graphNodes[ihel * ndiagrams + 564];
      gpuDiagram( &graph, &graphExec, &node565, &node564, diagram565, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node566 = graphNodes[ihel * ndiagrams + 565];
      gpuDiagram( &graph, &graphExec, &node566, &node565, diagram566, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node567 = graphNodes[ihel * ndiagrams + 566];
      gpuDiagram( &graph, &graphExec, &node567, &node566, diagram567, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node568 = graphNodes[ihel * ndiagrams + 567];
      gpuDiagram( &graph, &graphExec, &node568, &node567, diagram568, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node569 = graphNodes[ihel * ndiagrams + 568];
      gpuDiagram( &graph, &graphExec, &node569, &node568, diagram569, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node570 = graphNodes[ihel * ndiagrams + 569];
      gpuDiagram( &graph, &graphExec, &node570, &node569, diagram570, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node571 = graphNodes[ihel * ndiagrams + 570];
      gpuDiagram( &graph, &graphExec, &node571, &node570, diagram571, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node572 = graphNodes[ihel * ndiagrams + 571];
      gpuDiagram( &graph, &graphExec, &node572, &node571, diagram572, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node573 = graphNodes[ihel * ndiagrams + 572];
      gpuDiagram( &graph, &graphExec, &node573, &node572, diagram573, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node574 = graphNodes[ihel * ndiagrams + 573];
      gpuDiagram( &graph, &graphExec, &node574, &node573, diagram574, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node575 = graphNodes[ihel * ndiagrams + 574];
      gpuDiagram( &graph, &graphExec, &node575, &node574, diagram575, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node576 = graphNodes[ihel * ndiagrams + 575];
      gpuDiagram( &graph, &graphExec, &node576, &node575, diagram576, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node577 = graphNodes[ihel * ndiagrams + 576];
      gpuDiagram( &graph, &graphExec, &node577, &node576, diagram577, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node578 = graphNodes[ihel * ndiagrams + 577];
      gpuDiagram( &graph, &graphExec, &node578, &node577, diagram578, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node579 = graphNodes[ihel * ndiagrams + 578];
      gpuDiagram( &graph, &graphExec, &node579, &node578, diagram579, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node580 = graphNodes[ihel * ndiagrams + 579];
      gpuDiagram( &graph, &graphExec, &node580, &node579, diagram580, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node581 = graphNodes[ihel * ndiagrams + 580];
      gpuDiagram( &graph, &graphExec, &node581, &node580, diagram581, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node582 = graphNodes[ihel * ndiagrams + 581];
      gpuDiagram( &graph, &graphExec, &node582, &node581, diagram582, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node583 = graphNodes[ihel * ndiagrams + 582];
      gpuDiagram( &graph, &graphExec, &node583, &node582, diagram583, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node584 = graphNodes[ihel * ndiagrams + 583];
      gpuDiagram( &graph, &graphExec, &node584, &node583, diagram584, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node585 = graphNodes[ihel * ndiagrams + 584];
      gpuDiagram( &graph, &graphExec, &node585, &node584, diagram585, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node586 = graphNodes[ihel * ndiagrams + 585];
      gpuDiagram( &graph, &graphExec, &node586, &node585, diagram586, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node587 = graphNodes[ihel * ndiagrams + 586];
      gpuDiagram( &graph, &graphExec, &node587, &node586, diagram587, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node588 = graphNodes[ihel * ndiagrams + 587];
      gpuDiagram( &graph, &graphExec, &node588, &node587, diagram588, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node589 = graphNodes[ihel * ndiagrams + 588];
      gpuDiagram( &graph, &graphExec, &node589, &node588, diagram589, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node590 = graphNodes[ihel * ndiagrams + 589];
      gpuDiagram( &graph, &graphExec, &node590, &node589, diagram590, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node591 = graphNodes[ihel * ndiagrams + 590];
      gpuDiagram( &graph, &graphExec, &node591, &node590, diagram591, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node592 = graphNodes[ihel * ndiagrams + 591];
      gpuDiagram( &graph, &graphExec, &node592, &node591, diagram592, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node593 = graphNodes[ihel * ndiagrams + 592];
      gpuDiagram( &graph, &graphExec, &node593, &node592, diagram593, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node594 = graphNodes[ihel * ndiagrams + 593];
      gpuDiagram( &graph, &graphExec, &node594, &node593, diagram594, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node595 = graphNodes[ihel * ndiagrams + 594];
      gpuDiagram( &graph, &graphExec, &node595, &node594, diagram595, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node596 = graphNodes[ihel * ndiagrams + 595];
      gpuDiagram( &graph, &graphExec, &node596, &node595, diagram596, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node597 = graphNodes[ihel * ndiagrams + 596];
      gpuDiagram( &graph, &graphExec, &node597, &node596, diagram597, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node598 = graphNodes[ihel * ndiagrams + 597];
      gpuDiagram( &graph, &graphExec, &node598, &node597, diagram598, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node599 = graphNodes[ihel * ndiagrams + 598];
      gpuDiagram( &graph, &graphExec, &node599, &node598, diagram599, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node600 = graphNodes[ihel * ndiagrams + 599];
      gpuDiagram( &graph, &graphExec, &node600, &node599, diagram600, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node601 = graphNodes[ihel * ndiagrams + 600];
      gpuDiagram( &graph, &graphExec, &node601, &node600, diagram601, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node602 = graphNodes[ihel * ndiagrams + 601];
      gpuDiagram( &graph, &graphExec, &node602, &node601, diagram602, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node603 = graphNodes[ihel * ndiagrams + 602];
      gpuDiagram( &graph, &graphExec, &node603, &node602, diagram603, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node604 = graphNodes[ihel * ndiagrams + 603];
      gpuDiagram( &graph, &graphExec, &node604, &node603, diagram604, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node605 = graphNodes[ihel * ndiagrams + 604];
      gpuDiagram( &graph, &graphExec, &node605, &node604, diagram605, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node606 = graphNodes[ihel * ndiagrams + 605];
      gpuDiagram( &graph, &graphExec, &node606, &node605, diagram606, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node607 = graphNodes[ihel * ndiagrams + 606];
      gpuDiagram( &graph, &graphExec, &node607, &node606, diagram607, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node608 = graphNodes[ihel * ndiagrams + 607];
      gpuDiagram( &graph, &graphExec, &node608, &node607, diagram608, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node609 = graphNodes[ihel * ndiagrams + 608];
      gpuDiagram( &graph, &graphExec, &node609, &node608, diagram609, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node610 = graphNodes[ihel * ndiagrams + 609];
      gpuDiagram( &graph, &graphExec, &node610, &node609, diagram610, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node611 = graphNodes[ihel * ndiagrams + 610];
      gpuDiagram( &graph, &graphExec, &node611, &node610, diagram611, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node612 = graphNodes[ihel * ndiagrams + 611];
      gpuDiagram( &graph, &graphExec, &node612, &node611, diagram612, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node613 = graphNodes[ihel * ndiagrams + 612];
      gpuDiagram( &graph, &graphExec, &node613, &node612, diagram613, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node614 = graphNodes[ihel * ndiagrams + 613];
      gpuDiagram( &graph, &graphExec, &node614, &node613, diagram614, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node615 = graphNodes[ihel * ndiagrams + 614];
      gpuDiagram( &graph, &graphExec, &node615, &node614, diagram615, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node616 = graphNodes[ihel * ndiagrams + 615];
      gpuDiagram( &graph, &graphExec, &node616, &node615, diagram616, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node617 = graphNodes[ihel * ndiagrams + 616];
      gpuDiagram( &graph, &graphExec, &node617, &node616, diagram617, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node618 = graphNodes[ihel * ndiagrams + 617];
      gpuDiagram( &graph, &graphExec, &node618, &node617, diagram618, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node619 = graphNodes[ihel * ndiagrams + 618];
      gpuDiagram( &graph, &graphExec, &node619, &node618, diagram619, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node620 = graphNodes[ihel * ndiagrams + 619];
      gpuDiagram( &graph, &graphExec, &node620, &node619, diagram620, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node621 = graphNodes[ihel * ndiagrams + 620];
      gpuDiagram( &graph, &graphExec, &node621, &node620, diagram621, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node622 = graphNodes[ihel * ndiagrams + 621];
      gpuDiagram( &graph, &graphExec, &node622, &node621, diagram622, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node623 = graphNodes[ihel * ndiagrams + 622];
      gpuDiagram( &graph, &graphExec, &node623, &node622, diagram623, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node624 = graphNodes[ihel * ndiagrams + 623];
      gpuDiagram( &graph, &graphExec, &node624, &node623, diagram624, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node625 = graphNodes[ihel * ndiagrams + 624];
      gpuDiagram( &graph, &graphExec, &node625, &node624, diagram625, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node626 = graphNodes[ihel * ndiagrams + 625];
      gpuDiagram( &graph, &graphExec, &node626, &node625, diagram626, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node627 = graphNodes[ihel * ndiagrams + 626];
      gpuDiagram( &graph, &graphExec, &node627, &node626, diagram627, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node628 = graphNodes[ihel * ndiagrams + 627];
      gpuDiagram( &graph, &graphExec, &node628, &node627, diagram628, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node629 = graphNodes[ihel * ndiagrams + 628];
      gpuDiagram( &graph, &graphExec, &node629, &node628, diagram629, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node630 = graphNodes[ihel * ndiagrams + 629];
      gpuDiagram( &graph, &graphExec, &node630, &node629, diagram630, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node631 = graphNodes[ihel * ndiagrams + 630];
      gpuDiagram( &graph, &graphExec, &node631, &node630, diagram631, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node632 = graphNodes[ihel * ndiagrams + 631];
      gpuDiagram( &graph, &graphExec, &node632, &node631, diagram632, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node633 = graphNodes[ihel * ndiagrams + 632];
      gpuDiagram( &graph, &graphExec, &node633, &node632, diagram633, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node634 = graphNodes[ihel * ndiagrams + 633];
      gpuDiagram( &graph, &graphExec, &node634, &node633, diagram634, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node635 = graphNodes[ihel * ndiagrams + 634];
      gpuDiagram( &graph, &graphExec, &node635, &node634, diagram635, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node636 = graphNodes[ihel * ndiagrams + 635];
      gpuDiagram( &graph, &graphExec, &node636, &node635, diagram636, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node637 = graphNodes[ihel * ndiagrams + 636];
      gpuDiagram( &graph, &graphExec, &node637, &node636, diagram637, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node638 = graphNodes[ihel * ndiagrams + 637];
      gpuDiagram( &graph, &graphExec, &node638, &node637, diagram638, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node639 = graphNodes[ihel * ndiagrams + 638];
      gpuDiagram( &graph, &graphExec, &node639, &node638, diagram639, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node640 = graphNodes[ihel * ndiagrams + 639];
      gpuDiagram( &graph, &graphExec, &node640, &node639, diagram640, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node641 = graphNodes[ihel * ndiagrams + 640];
      gpuDiagram( &graph, &graphExec, &node641, &node640, diagram641, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node642 = graphNodes[ihel * ndiagrams + 641];
      gpuDiagram( &graph, &graphExec, &node642, &node641, diagram642, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node643 = graphNodes[ihel * ndiagrams + 642];
      gpuDiagram( &graph, &graphExec, &node643, &node642, diagram643, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node644 = graphNodes[ihel * ndiagrams + 643];
      gpuDiagram( &graph, &graphExec, &node644, &node643, diagram644, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node645 = graphNodes[ihel * ndiagrams + 644];
      gpuDiagram( &graph, &graphExec, &node645, &node644, diagram645, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node646 = graphNodes[ihel * ndiagrams + 645];
      gpuDiagram( &graph, &graphExec, &node646, &node645, diagram646, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node647 = graphNodes[ihel * ndiagrams + 646];
      gpuDiagram( &graph, &graphExec, &node647, &node646, diagram647, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node648 = graphNodes[ihel * ndiagrams + 647];
      gpuDiagram( &graph, &graphExec, &node648, &node647, diagram648, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node649 = graphNodes[ihel * ndiagrams + 648];
      gpuDiagram( &graph, &graphExec, &node649, &node648, diagram649, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node650 = graphNodes[ihel * ndiagrams + 649];
      gpuDiagram( &graph, &graphExec, &node650, &node649, diagram650, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node651 = graphNodes[ihel * ndiagrams + 650];
      gpuDiagram( &graph, &graphExec, &node651, &node650, diagram651, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node652 = graphNodes[ihel * ndiagrams + 651];
      gpuDiagram( &graph, &graphExec, &node652, &node651, diagram652, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node653 = graphNodes[ihel * ndiagrams + 652];
      gpuDiagram( &graph, &graphExec, &node653, &node652, diagram653, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node654 = graphNodes[ihel * ndiagrams + 653];
      gpuDiagram( &graph, &graphExec, &node654, &node653, diagram654, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node655 = graphNodes[ihel * ndiagrams + 654];
      gpuDiagram( &graph, &graphExec, &node655, &node654, diagram655, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node656 = graphNodes[ihel * ndiagrams + 655];
      gpuDiagram( &graph, &graphExec, &node656, &node655, diagram656, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node657 = graphNodes[ihel * ndiagrams + 656];
      gpuDiagram( &graph, &graphExec, &node657, &node656, diagram657, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node658 = graphNodes[ihel * ndiagrams + 657];
      gpuDiagram( &graph, &graphExec, &node658, &node657, diagram658, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node659 = graphNodes[ihel * ndiagrams + 658];
      gpuDiagram( &graph, &graphExec, &node659, &node658, diagram659, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node660 = graphNodes[ihel * ndiagrams + 659];
      gpuDiagram( &graph, &graphExec, &node660, &node659, diagram660, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node661 = graphNodes[ihel * ndiagrams + 660];
      gpuDiagram( &graph, &graphExec, &node661, &node660, diagram661, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node662 = graphNodes[ihel * ndiagrams + 661];
      gpuDiagram( &graph, &graphExec, &node662, &node661, diagram662, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node663 = graphNodes[ihel * ndiagrams + 662];
      gpuDiagram( &graph, &graphExec, &node663, &node662, diagram663, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node664 = graphNodes[ihel * ndiagrams + 663];
      gpuDiagram( &graph, &graphExec, &node664, &node663, diagram664, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node665 = graphNodes[ihel * ndiagrams + 664];
      gpuDiagram( &graph, &graphExec, &node665, &node664, diagram665, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node666 = graphNodes[ihel * ndiagrams + 665];
      gpuDiagram( &graph, &graphExec, &node666, &node665, diagram666, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node667 = graphNodes[ihel * ndiagrams + 666];
      gpuDiagram( &graph, &graphExec, &node667, &node666, diagram667, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node668 = graphNodes[ihel * ndiagrams + 667];
      gpuDiagram( &graph, &graphExec, &node668, &node667, diagram668, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node669 = graphNodes[ihel * ndiagrams + 668];
      gpuDiagram( &graph, &graphExec, &node669, &node668, diagram669, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node670 = graphNodes[ihel * ndiagrams + 669];
      gpuDiagram( &graph, &graphExec, &node670, &node669, diagram670, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node671 = graphNodes[ihel * ndiagrams + 670];
      gpuDiagram( &graph, &graphExec, &node671, &node670, diagram671, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node672 = graphNodes[ihel * ndiagrams + 671];
      gpuDiagram( &graph, &graphExec, &node672, &node671, diagram672, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node673 = graphNodes[ihel * ndiagrams + 672];
      gpuDiagram( &graph, &graphExec, &node673, &node672, diagram673, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node674 = graphNodes[ihel * ndiagrams + 673];
      gpuDiagram( &graph, &graphExec, &node674, &node673, diagram674, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node675 = graphNodes[ihel * ndiagrams + 674];
      gpuDiagram( &graph, &graphExec, &node675, &node674, diagram675, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node676 = graphNodes[ihel * ndiagrams + 675];
      gpuDiagram( &graph, &graphExec, &node676, &node675, diagram676, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node677 = graphNodes[ihel * ndiagrams + 676];
      gpuDiagram( &graph, &graphExec, &node677, &node676, diagram677, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node678 = graphNodes[ihel * ndiagrams + 677];
      gpuDiagram( &graph, &graphExec, &node678, &node677, diagram678, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node679 = graphNodes[ihel * ndiagrams + 678];
      gpuDiagram( &graph, &graphExec, &node679, &node678, diagram679, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node680 = graphNodes[ihel * ndiagrams + 679];
      gpuDiagram( &graph, &graphExec, &node680, &node679, diagram680, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node681 = graphNodes[ihel * ndiagrams + 680];
      gpuDiagram( &graph, &graphExec, &node681, &node680, diagram681, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node682 = graphNodes[ihel * ndiagrams + 681];
      gpuDiagram( &graph, &graphExec, &node682, &node681, diagram682, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node683 = graphNodes[ihel * ndiagrams + 682];
      gpuDiagram( &graph, &graphExec, &node683, &node682, diagram683, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node684 = graphNodes[ihel * ndiagrams + 683];
      gpuDiagram( &graph, &graphExec, &node684, &node683, diagram684, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node685 = graphNodes[ihel * ndiagrams + 684];
      gpuDiagram( &graph, &graphExec, &node685, &node684, diagram685, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node686 = graphNodes[ihel * ndiagrams + 685];
      gpuDiagram( &graph, &graphExec, &node686, &node685, diagram686, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node687 = graphNodes[ihel * ndiagrams + 686];
      gpuDiagram( &graph, &graphExec, &node687, &node686, diagram687, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node688 = graphNodes[ihel * ndiagrams + 687];
      gpuDiagram( &graph, &graphExec, &node688, &node687, diagram688, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node689 = graphNodes[ihel * ndiagrams + 688];
      gpuDiagram( &graph, &graphExec, &node689, &node688, diagram689, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node690 = graphNodes[ihel * ndiagrams + 689];
      gpuDiagram( &graph, &graphExec, &node690, &node689, diagram690, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node691 = graphNodes[ihel * ndiagrams + 690];
      gpuDiagram( &graph, &graphExec, &node691, &node690, diagram691, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node692 = graphNodes[ihel * ndiagrams + 691];
      gpuDiagram( &graph, &graphExec, &node692, &node691, diagram692, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node693 = graphNodes[ihel * ndiagrams + 692];
      gpuDiagram( &graph, &graphExec, &node693, &node692, diagram693, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node694 = graphNodes[ihel * ndiagrams + 693];
      gpuDiagram( &graph, &graphExec, &node694, &node693, diagram694, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node695 = graphNodes[ihel * ndiagrams + 694];
      gpuDiagram( &graph, &graphExec, &node695, &node694, diagram695, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node696 = graphNodes[ihel * ndiagrams + 695];
      gpuDiagram( &graph, &graphExec, &node696, &node695, diagram696, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node697 = graphNodes[ihel * ndiagrams + 696];
      gpuDiagram( &graph, &graphExec, &node697, &node696, diagram697, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node698 = graphNodes[ihel * ndiagrams + 697];
      gpuDiagram( &graph, &graphExec, &node698, &node697, diagram698, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node699 = graphNodes[ihel * ndiagrams + 698];
      gpuDiagram( &graph, &graphExec, &node699, &node698, diagram699, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node700 = graphNodes[ihel * ndiagrams + 699];
      gpuDiagram( &graph, &graphExec, &node700, &node699, diagram700, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node701 = graphNodes[ihel * ndiagrams + 700];
      gpuDiagram( &graph, &graphExec, &node701, &node700, diagram701, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node702 = graphNodes[ihel * ndiagrams + 701];
      gpuDiagram( &graph, &graphExec, &node702, &node701, diagram702, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node703 = graphNodes[ihel * ndiagrams + 702];
      gpuDiagram( &graph, &graphExec, &node703, &node702, diagram703, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node704 = graphNodes[ihel * ndiagrams + 703];
      gpuDiagram( &graph, &graphExec, &node704, &node703, diagram704, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node705 = graphNodes[ihel * ndiagrams + 704];
      gpuDiagram( &graph, &graphExec, &node705, &node704, diagram705, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node706 = graphNodes[ihel * ndiagrams + 705];
      gpuDiagram( &graph, &graphExec, &node706, &node705, diagram706, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node707 = graphNodes[ihel * ndiagrams + 706];
      gpuDiagram( &graph, &graphExec, &node707, &node706, diagram707, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node708 = graphNodes[ihel * ndiagrams + 707];
      gpuDiagram( &graph, &graphExec, &node708, &node707, diagram708, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node709 = graphNodes[ihel * ndiagrams + 708];
      gpuDiagram( &graph, &graphExec, &node709, &node708, diagram709, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node710 = graphNodes[ihel * ndiagrams + 709];
      gpuDiagram( &graph, &graphExec, &node710, &node709, diagram710, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node711 = graphNodes[ihel * ndiagrams + 710];
      gpuDiagram( &graph, &graphExec, &node711, &node710, diagram711, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node712 = graphNodes[ihel * ndiagrams + 711];
      gpuDiagram( &graph, &graphExec, &node712, &node711, diagram712, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node713 = graphNodes[ihel * ndiagrams + 712];
      gpuDiagram( &graph, &graphExec, &node713, &node712, diagram713, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node714 = graphNodes[ihel * ndiagrams + 713];
      gpuDiagram( &graph, &graphExec, &node714, &node713, diagram714, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node715 = graphNodes[ihel * ndiagrams + 714];
      gpuDiagram( &graph, &graphExec, &node715, &node714, diagram715, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node716 = graphNodes[ihel * ndiagrams + 715];
      gpuDiagram( &graph, &graphExec, &node716, &node715, diagram716, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node717 = graphNodes[ihel * ndiagrams + 716];
      gpuDiagram( &graph, &graphExec, &node717, &node716, diagram717, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node718 = graphNodes[ihel * ndiagrams + 717];
      gpuDiagram( &graph, &graphExec, &node718, &node717, diagram718, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node719 = graphNodes[ihel * ndiagrams + 718];
      gpuDiagram( &graph, &graphExec, &node719, &node718, diagram719, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node720 = graphNodes[ihel * ndiagrams + 719];
      gpuDiagram( &graph, &graphExec, &node720, &node719, diagram720, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node721 = graphNodes[ihel * ndiagrams + 720];
      gpuDiagram( &graph, &graphExec, &node721, &node720, diagram721, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node722 = graphNodes[ihel * ndiagrams + 721];
      gpuDiagram( &graph, &graphExec, &node722, &node721, diagram722, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node723 = graphNodes[ihel * ndiagrams + 722];
      gpuDiagram( &graph, &graphExec, &node723, &node722, diagram723, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node724 = graphNodes[ihel * ndiagrams + 723];
      gpuDiagram( &graph, &graphExec, &node724, &node723, diagram724, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node725 = graphNodes[ihel * ndiagrams + 724];
      gpuDiagram( &graph, &graphExec, &node725, &node724, diagram725, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node726 = graphNodes[ihel * ndiagrams + 725];
      gpuDiagram( &graph, &graphExec, &node726, &node725, diagram726, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node727 = graphNodes[ihel * ndiagrams + 726];
      gpuDiagram( &graph, &graphExec, &node727, &node726, diagram727, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node728 = graphNodes[ihel * ndiagrams + 727];
      gpuDiagram( &graph, &graphExec, &node728, &node727, diagram728, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node729 = graphNodes[ihel * ndiagrams + 728];
      gpuDiagram( &graph, &graphExec, &node729, &node728, diagram729, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node730 = graphNodes[ihel * ndiagrams + 729];
      gpuDiagram( &graph, &graphExec, &node730, &node729, diagram730, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node731 = graphNodes[ihel * ndiagrams + 730];
      gpuDiagram( &graph, &graphExec, &node731, &node730, diagram731, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node732 = graphNodes[ihel * ndiagrams + 731];
      gpuDiagram( &graph, &graphExec, &node732, &node731, diagram732, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node733 = graphNodes[ihel * ndiagrams + 732];
      gpuDiagram( &graph, &graphExec, &node733, &node732, diagram733, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node734 = graphNodes[ihel * ndiagrams + 733];
      gpuDiagram( &graph, &graphExec, &node734, &node733, diagram734, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node735 = graphNodes[ihel * ndiagrams + 734];
      gpuDiagram( &graph, &graphExec, &node735, &node734, diagram735, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node736 = graphNodes[ihel * ndiagrams + 735];
      gpuDiagram( &graph, &graphExec, &node736, &node735, diagram736, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node737 = graphNodes[ihel * ndiagrams + 736];
      gpuDiagram( &graph, &graphExec, &node737, &node736, diagram737, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node738 = graphNodes[ihel * ndiagrams + 737];
      gpuDiagram( &graph, &graphExec, &node738, &node737, diagram738, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node739 = graphNodes[ihel * ndiagrams + 738];
      gpuDiagram( &graph, &graphExec, &node739, &node738, diagram739, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node740 = graphNodes[ihel * ndiagrams + 739];
      gpuDiagram( &graph, &graphExec, &node740, &node739, diagram740, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node741 = graphNodes[ihel * ndiagrams + 740];
      gpuDiagram( &graph, &graphExec, &node741, &node740, diagram741, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node742 = graphNodes[ihel * ndiagrams + 741];
      gpuDiagram( &graph, &graphExec, &node742, &node741, diagram742, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node743 = graphNodes[ihel * ndiagrams + 742];
      gpuDiagram( &graph, &graphExec, &node743, &node742, diagram743, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node744 = graphNodes[ihel * ndiagrams + 743];
      gpuDiagram( &graph, &graphExec, &node744, &node743, diagram744, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node745 = graphNodes[ihel * ndiagrams + 744];
      gpuDiagram( &graph, &graphExec, &node745, &node744, diagram745, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node746 = graphNodes[ihel * ndiagrams + 745];
      gpuDiagram( &graph, &graphExec, &node746, &node745, diagram746, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node747 = graphNodes[ihel * ndiagrams + 746];
      gpuDiagram( &graph, &graphExec, &node747, &node746, diagram747, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node748 = graphNodes[ihel * ndiagrams + 747];
      gpuDiagram( &graph, &graphExec, &node748, &node747, diagram748, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node749 = graphNodes[ihel * ndiagrams + 748];
      gpuDiagram( &graph, &graphExec, &node749, &node748, diagram749, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node750 = graphNodes[ihel * ndiagrams + 749];
      gpuDiagram( &graph, &graphExec, &node750, &node749, diagram750, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node751 = graphNodes[ihel * ndiagrams + 750];
      gpuDiagram( &graph, &graphExec, &node751, &node750, diagram751, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node752 = graphNodes[ihel * ndiagrams + 751];
      gpuDiagram( &graph, &graphExec, &node752, &node751, diagram752, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node753 = graphNodes[ihel * ndiagrams + 752];
      gpuDiagram( &graph, &graphExec, &node753, &node752, diagram753, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node754 = graphNodes[ihel * ndiagrams + 753];
      gpuDiagram( &graph, &graphExec, &node754, &node753, diagram754, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node755 = graphNodes[ihel * ndiagrams + 754];
      gpuDiagram( &graph, &graphExec, &node755, &node754, diagram755, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node756 = graphNodes[ihel * ndiagrams + 755];
      gpuDiagram( &graph, &graphExec, &node756, &node755, diagram756, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node757 = graphNodes[ihel * ndiagrams + 756];
      gpuDiagram( &graph, &graphExec, &node757, &node756, diagram757, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node758 = graphNodes[ihel * ndiagrams + 757];
      gpuDiagram( &graph, &graphExec, &node758, &node757, diagram758, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node759 = graphNodes[ihel * ndiagrams + 758];
      gpuDiagram( &graph, &graphExec, &node759, &node758, diagram759, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node760 = graphNodes[ihel * ndiagrams + 759];
      gpuDiagram( &graph, &graphExec, &node760, &node759, diagram760, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node761 = graphNodes[ihel * ndiagrams + 760];
      gpuDiagram( &graph, &graphExec, &node761, &node760, diagram761, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node762 = graphNodes[ihel * ndiagrams + 761];
      gpuDiagram( &graph, &graphExec, &node762, &node761, diagram762, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node763 = graphNodes[ihel * ndiagrams + 762];
      gpuDiagram( &graph, &graphExec, &node763, &node762, diagram763, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node764 = graphNodes[ihel * ndiagrams + 763];
      gpuDiagram( &graph, &graphExec, &node764, &node763, diagram764, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node765 = graphNodes[ihel * ndiagrams + 764];
      gpuDiagram( &graph, &graphExec, &node765, &node764, diagram765, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node766 = graphNodes[ihel * ndiagrams + 765];
      gpuDiagram( &graph, &graphExec, &node766, &node765, diagram766, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node767 = graphNodes[ihel * ndiagrams + 766];
      gpuDiagram( &graph, &graphExec, &node767, &node766, diagram767, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node768 = graphNodes[ihel * ndiagrams + 767];
      gpuDiagram( &graph, &graphExec, &node768, &node767, diagram768, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node769 = graphNodes[ihel * ndiagrams + 768];
      gpuDiagram( &graph, &graphExec, &node769, &node768, diagram769, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node770 = graphNodes[ihel * ndiagrams + 769];
      gpuDiagram( &graph, &graphExec, &node770, &node769, diagram770, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node771 = graphNodes[ihel * ndiagrams + 770];
      gpuDiagram( &graph, &graphExec, &node771, &node770, diagram771, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node772 = graphNodes[ihel * ndiagrams + 771];
      gpuDiagram( &graph, &graphExec, &node772, &node771, diagram772, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node773 = graphNodes[ihel * ndiagrams + 772];
      gpuDiagram( &graph, &graphExec, &node773, &node772, diagram773, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node774 = graphNodes[ihel * ndiagrams + 773];
      gpuDiagram( &graph, &graphExec, &node774, &node773, diagram774, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node775 = graphNodes[ihel * ndiagrams + 774];
      gpuDiagram( &graph, &graphExec, &node775, &node774, diagram775, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node776 = graphNodes[ihel * ndiagrams + 775];
      gpuDiagram( &graph, &graphExec, &node776, &node775, diagram776, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node777 = graphNodes[ihel * ndiagrams + 776];
      gpuDiagram( &graph, &graphExec, &node777, &node776, diagram777, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node778 = graphNodes[ihel * ndiagrams + 777];
      gpuDiagram( &graph, &graphExec, &node778, &node777, diagram778, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node779 = graphNodes[ihel * ndiagrams + 778];
      gpuDiagram( &graph, &graphExec, &node779, &node778, diagram779, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node780 = graphNodes[ihel * ndiagrams + 779];
      gpuDiagram( &graph, &graphExec, &node780, &node779, diagram780, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node781 = graphNodes[ihel * ndiagrams + 780];
      gpuDiagram( &graph, &graphExec, &node781, &node780, diagram781, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node782 = graphNodes[ihel * ndiagrams + 781];
      gpuDiagram( &graph, &graphExec, &node782, &node781, diagram782, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node783 = graphNodes[ihel * ndiagrams + 782];
      gpuDiagram( &graph, &graphExec, &node783, &node782, diagram783, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node784 = graphNodes[ihel * ndiagrams + 783];
      gpuDiagram( &graph, &graphExec, &node784, &node783, diagram784, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node785 = graphNodes[ihel * ndiagrams + 784];
      gpuDiagram( &graph, &graphExec, &node785, &node784, diagram785, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node786 = graphNodes[ihel * ndiagrams + 785];
      gpuDiagram( &graph, &graphExec, &node786, &node785, diagram786, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node787 = graphNodes[ihel * ndiagrams + 786];
      gpuDiagram( &graph, &graphExec, &node787, &node786, diagram787, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node788 = graphNodes[ihel * ndiagrams + 787];
      gpuDiagram( &graph, &graphExec, &node788, &node787, diagram788, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node789 = graphNodes[ihel * ndiagrams + 788];
      gpuDiagram( &graph, &graphExec, &node789, &node788, diagram789, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node790 = graphNodes[ihel * ndiagrams + 789];
      gpuDiagram( &graph, &graphExec, &node790, &node789, diagram790, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node791 = graphNodes[ihel * ndiagrams + 790];
      gpuDiagram( &graph, &graphExec, &node791, &node790, diagram791, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node792 = graphNodes[ihel * ndiagrams + 791];
      gpuDiagram( &graph, &graphExec, &node792, &node791, diagram792, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node793 = graphNodes[ihel * ndiagrams + 792];
      gpuDiagram( &graph, &graphExec, &node793, &node792, diagram793, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node794 = graphNodes[ihel * ndiagrams + 793];
      gpuDiagram( &graph, &graphExec, &node794, &node793, diagram794, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node795 = graphNodes[ihel * ndiagrams + 794];
      gpuDiagram( &graph, &graphExec, &node795, &node794, diagram795, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node796 = graphNodes[ihel * ndiagrams + 795];
      gpuDiagram( &graph, &graphExec, &node796, &node795, diagram796, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node797 = graphNodes[ihel * ndiagrams + 796];
      gpuDiagram( &graph, &graphExec, &node797, &node796, diagram797, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node798 = graphNodes[ihel * ndiagrams + 797];
      gpuDiagram( &graph, &graphExec, &node798, &node797, diagram798, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node799 = graphNodes[ihel * ndiagrams + 798];
      gpuDiagram( &graph, &graphExec, &node799, &node798, diagram799, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node800 = graphNodes[ihel * ndiagrams + 799];
      gpuDiagram( &graph, &graphExec, &node800, &node799, diagram800, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node801 = graphNodes[ihel * ndiagrams + 800];
      gpuDiagram( &graph, &graphExec, &node801, &node800, diagram801, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node802 = graphNodes[ihel * ndiagrams + 801];
      gpuDiagram( &graph, &graphExec, &node802, &node801, diagram802, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node803 = graphNodes[ihel * ndiagrams + 802];
      gpuDiagram( &graph, &graphExec, &node803, &node802, diagram803, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node804 = graphNodes[ihel * ndiagrams + 803];
      gpuDiagram( &graph, &graphExec, &node804, &node803, diagram804, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node805 = graphNodes[ihel * ndiagrams + 804];
      gpuDiagram( &graph, &graphExec, &node805, &node804, diagram805, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node806 = graphNodes[ihel * ndiagrams + 805];
      gpuDiagram( &graph, &graphExec, &node806, &node805, diagram806, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node807 = graphNodes[ihel * ndiagrams + 806];
      gpuDiagram( &graph, &graphExec, &node807, &node806, diagram807, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node808 = graphNodes[ihel * ndiagrams + 807];
      gpuDiagram( &graph, &graphExec, &node808, &node807, diagram808, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node809 = graphNodes[ihel * ndiagrams + 808];
      gpuDiagram( &graph, &graphExec, &node809, &node808, diagram809, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node810 = graphNodes[ihel * ndiagrams + 809];
      gpuDiagram( &graph, &graphExec, &node810, &node809, diagram810, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node811 = graphNodes[ihel * ndiagrams + 810];
      gpuDiagram( &graph, &graphExec, &node811, &node810, diagram811, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node812 = graphNodes[ihel * ndiagrams + 811];
      gpuDiagram( &graph, &graphExec, &node812, &node811, diagram812, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node813 = graphNodes[ihel * ndiagrams + 812];
      gpuDiagram( &graph, &graphExec, &node813, &node812, diagram813, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node814 = graphNodes[ihel * ndiagrams + 813];
      gpuDiagram( &graph, &graphExec, &node814, &node813, diagram814, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node815 = graphNodes[ihel * ndiagrams + 814];
      gpuDiagram( &graph, &graphExec, &node815, &node814, diagram815, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node816 = graphNodes[ihel * ndiagrams + 815];
      gpuDiagram( &graph, &graphExec, &node816, &node815, diagram816, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node817 = graphNodes[ihel * ndiagrams + 816];
      gpuDiagram( &graph, &graphExec, &node817, &node816, diagram817, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node818 = graphNodes[ihel * ndiagrams + 817];
      gpuDiagram( &graph, &graphExec, &node818, &node817, diagram818, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node819 = graphNodes[ihel * ndiagrams + 818];
      gpuDiagram( &graph, &graphExec, &node819, &node818, diagram819, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node820 = graphNodes[ihel * ndiagrams + 819];
      gpuDiagram( &graph, &graphExec, &node820, &node819, diagram820, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node821 = graphNodes[ihel * ndiagrams + 820];
      gpuDiagram( &graph, &graphExec, &node821, &node820, diagram821, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node822 = graphNodes[ihel * ndiagrams + 821];
      gpuDiagram( &graph, &graphExec, &node822, &node821, diagram822, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node823 = graphNodes[ihel * ndiagrams + 822];
      gpuDiagram( &graph, &graphExec, &node823, &node822, diagram823, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node824 = graphNodes[ihel * ndiagrams + 823];
      gpuDiagram( &graph, &graphExec, &node824, &node823, diagram824, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node825 = graphNodes[ihel * ndiagrams + 824];
      gpuDiagram( &graph, &graphExec, &node825, &node824, diagram825, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node826 = graphNodes[ihel * ndiagrams + 825];
      gpuDiagram( &graph, &graphExec, &node826, &node825, diagram826, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node827 = graphNodes[ihel * ndiagrams + 826];
      gpuDiagram( &graph, &graphExec, &node827, &node826, diagram827, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node828 = graphNodes[ihel * ndiagrams + 827];
      gpuDiagram( &graph, &graphExec, &node828, &node827, diagram828, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node829 = graphNodes[ihel * ndiagrams + 828];
      gpuDiagram( &graph, &graphExec, &node829, &node828, diagram829, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node830 = graphNodes[ihel * ndiagrams + 829];
      gpuDiagram( &graph, &graphExec, &node830, &node829, diagram830, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node831 = graphNodes[ihel * ndiagrams + 830];
      gpuDiagram( &graph, &graphExec, &node831, &node830, diagram831, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node832 = graphNodes[ihel * ndiagrams + 831];
      gpuDiagram( &graph, &graphExec, &node832, &node831, diagram832, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node833 = graphNodes[ihel * ndiagrams + 832];
      gpuDiagram( &graph, &graphExec, &node833, &node832, diagram833, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node834 = graphNodes[ihel * ndiagrams + 833];
      gpuDiagram( &graph, &graphExec, &node834, &node833, diagram834, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node835 = graphNodes[ihel * ndiagrams + 834];
      gpuDiagram( &graph, &graphExec, &node835, &node834, diagram835, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node836 = graphNodes[ihel * ndiagrams + 835];
      gpuDiagram( &graph, &graphExec, &node836, &node835, diagram836, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node837 = graphNodes[ihel * ndiagrams + 836];
      gpuDiagram( &graph, &graphExec, &node837, &node836, diagram837, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node838 = graphNodes[ihel * ndiagrams + 837];
      gpuDiagram( &graph, &graphExec, &node838, &node837, diagram838, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node839 = graphNodes[ihel * ndiagrams + 838];
      gpuDiagram( &graph, &graphExec, &node839, &node838, diagram839, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node840 = graphNodes[ihel * ndiagrams + 839];
      gpuDiagram( &graph, &graphExec, &node840, &node839, diagram840, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node841 = graphNodes[ihel * ndiagrams + 840];
      gpuDiagram( &graph, &graphExec, &node841, &node840, diagram841, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node842 = graphNodes[ihel * ndiagrams + 841];
      gpuDiagram( &graph, &graphExec, &node842, &node841, diagram842, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node843 = graphNodes[ihel * ndiagrams + 842];
      gpuDiagram( &graph, &graphExec, &node843, &node842, diagram843, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node844 = graphNodes[ihel * ndiagrams + 843];
      gpuDiagram( &graph, &graphExec, &node844, &node843, diagram844, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node845 = graphNodes[ihel * ndiagrams + 844];
      gpuDiagram( &graph, &graphExec, &node845, &node844, diagram845, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node846 = graphNodes[ihel * ndiagrams + 845];
      gpuDiagram( &graph, &graphExec, &node846, &node845, diagram846, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node847 = graphNodes[ihel * ndiagrams + 846];
      gpuDiagram( &graph, &graphExec, &node847, &node846, diagram847, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node848 = graphNodes[ihel * ndiagrams + 847];
      gpuDiagram( &graph, &graphExec, &node848, &node847, diagram848, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node849 = graphNodes[ihel * ndiagrams + 848];
      gpuDiagram( &graph, &graphExec, &node849, &node848, diagram849, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node850 = graphNodes[ihel * ndiagrams + 849];
      gpuDiagram( &graph, &graphExec, &node850, &node849, diagram850, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node851 = graphNodes[ihel * ndiagrams + 850];
      gpuDiagram( &graph, &graphExec, &node851, &node850, diagram851, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node852 = graphNodes[ihel * ndiagrams + 851];
      gpuDiagram( &graph, &graphExec, &node852, &node851, diagram852, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node853 = graphNodes[ihel * ndiagrams + 852];
      gpuDiagram( &graph, &graphExec, &node853, &node852, diagram853, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node854 = graphNodes[ihel * ndiagrams + 853];
      gpuDiagram( &graph, &graphExec, &node854, &node853, diagram854, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node855 = graphNodes[ihel * ndiagrams + 854];
      gpuDiagram( &graph, &graphExec, &node855, &node854, diagram855, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node856 = graphNodes[ihel * ndiagrams + 855];
      gpuDiagram( &graph, &graphExec, &node856, &node855, diagram856, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node857 = graphNodes[ihel * ndiagrams + 856];
      gpuDiagram( &graph, &graphExec, &node857, &node856, diagram857, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node858 = graphNodes[ihel * ndiagrams + 857];
      gpuDiagram( &graph, &graphExec, &node858, &node857, diagram858, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node859 = graphNodes[ihel * ndiagrams + 858];
      gpuDiagram( &graph, &graphExec, &node859, &node858, diagram859, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node860 = graphNodes[ihel * ndiagrams + 859];
      gpuDiagram( &graph, &graphExec, &node860, &node859, diagram860, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node861 = graphNodes[ihel * ndiagrams + 860];
      gpuDiagram( &graph, &graphExec, &node861, &node860, diagram861, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node862 = graphNodes[ihel * ndiagrams + 861];
      gpuDiagram( &graph, &graphExec, &node862, &node861, diagram862, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node863 = graphNodes[ihel * ndiagrams + 862];
      gpuDiagram( &graph, &graphExec, &node863, &node862, diagram863, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node864 = graphNodes[ihel * ndiagrams + 863];
      gpuDiagram( &graph, &graphExec, &node864, &node863, diagram864, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node865 = graphNodes[ihel * ndiagrams + 864];
      gpuDiagram( &graph, &graphExec, &node865, &node864, diagram865, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node866 = graphNodes[ihel * ndiagrams + 865];
      gpuDiagram( &graph, &graphExec, &node866, &node865, diagram866, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node867 = graphNodes[ihel * ndiagrams + 866];
      gpuDiagram( &graph, &graphExec, &node867, &node866, diagram867, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node868 = graphNodes[ihel * ndiagrams + 867];
      gpuDiagram( &graph, &graphExec, &node868, &node867, diagram868, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node869 = graphNodes[ihel * ndiagrams + 868];
      gpuDiagram( &graph, &graphExec, &node869, &node868, diagram869, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node870 = graphNodes[ihel * ndiagrams + 869];
      gpuDiagram( &graph, &graphExec, &node870, &node869, diagram870, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node871 = graphNodes[ihel * ndiagrams + 870];
      gpuDiagram( &graph, &graphExec, &node871, &node870, diagram871, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node872 = graphNodes[ihel * ndiagrams + 871];
      gpuDiagram( &graph, &graphExec, &node872, &node871, diagram872, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node873 = graphNodes[ihel * ndiagrams + 872];
      gpuDiagram( &graph, &graphExec, &node873, &node872, diagram873, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node874 = graphNodes[ihel * ndiagrams + 873];
      gpuDiagram( &graph, &graphExec, &node874, &node873, diagram874, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node875 = graphNodes[ihel * ndiagrams + 874];
      gpuDiagram( &graph, &graphExec, &node875, &node874, diagram875, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node876 = graphNodes[ihel * ndiagrams + 875];
      gpuDiagram( &graph, &graphExec, &node876, &node875, diagram876, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node877 = graphNodes[ihel * ndiagrams + 876];
      gpuDiagram( &graph, &graphExec, &node877, &node876, diagram877, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node878 = graphNodes[ihel * ndiagrams + 877];
      gpuDiagram( &graph, &graphExec, &node878, &node877, diagram878, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node879 = graphNodes[ihel * ndiagrams + 878];
      gpuDiagram( &graph, &graphExec, &node879, &node878, diagram879, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node880 = graphNodes[ihel * ndiagrams + 879];
      gpuDiagram( &graph, &graphExec, &node880, &node879, diagram880, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node881 = graphNodes[ihel * ndiagrams + 880];
      gpuDiagram( &graph, &graphExec, &node881, &node880, diagram881, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node882 = graphNodes[ihel * ndiagrams + 881];
      gpuDiagram( &graph, &graphExec, &node882, &node881, diagram882, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node883 = graphNodes[ihel * ndiagrams + 882];
      gpuDiagram( &graph, &graphExec, &node883, &node882, diagram883, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node884 = graphNodes[ihel * ndiagrams + 883];
      gpuDiagram( &graph, &graphExec, &node884, &node883, diagram884, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node885 = graphNodes[ihel * ndiagrams + 884];
      gpuDiagram( &graph, &graphExec, &node885, &node884, diagram885, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node886 = graphNodes[ihel * ndiagrams + 885];
      gpuDiagram( &graph, &graphExec, &node886, &node885, diagram886, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node887 = graphNodes[ihel * ndiagrams + 886];
      gpuDiagram( &graph, &graphExec, &node887, &node886, diagram887, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node888 = graphNodes[ihel * ndiagrams + 887];
      gpuDiagram( &graph, &graphExec, &node888, &node887, diagram888, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node889 = graphNodes[ihel * ndiagrams + 888];
      gpuDiagram( &graph, &graphExec, &node889, &node888, diagram889, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node890 = graphNodes[ihel * ndiagrams + 889];
      gpuDiagram( &graph, &graphExec, &node890, &node889, diagram890, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node891 = graphNodes[ihel * ndiagrams + 890];
      gpuDiagram( &graph, &graphExec, &node891, &node890, diagram891, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node892 = graphNodes[ihel * ndiagrams + 891];
      gpuDiagram( &graph, &graphExec, &node892, &node891, diagram892, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node893 = graphNodes[ihel * ndiagrams + 892];
      gpuDiagram( &graph, &graphExec, &node893, &node892, diagram893, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node894 = graphNodes[ihel * ndiagrams + 893];
      gpuDiagram( &graph, &graphExec, &node894, &node893, diagram894, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node895 = graphNodes[ihel * ndiagrams + 894];
      gpuDiagram( &graph, &graphExec, &node895, &node894, diagram895, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node896 = graphNodes[ihel * ndiagrams + 895];
      gpuDiagram( &graph, &graphExec, &node896, &node895, diagram896, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node897 = graphNodes[ihel * ndiagrams + 896];
      gpuDiagram( &graph, &graphExec, &node897, &node896, diagram897, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node898 = graphNodes[ihel * ndiagrams + 897];
      gpuDiagram( &graph, &graphExec, &node898, &node897, diagram898, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node899 = graphNodes[ihel * ndiagrams + 898];
      gpuDiagram( &graph, &graphExec, &node899, &node898, diagram899, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node900 = graphNodes[ihel * ndiagrams + 899];
      gpuDiagram( &graph, &graphExec, &node900, &node899, diagram900, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node901 = graphNodes[ihel * ndiagrams + 900];
      gpuDiagram( &graph, &graphExec, &node901, &node900, diagram901, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node902 = graphNodes[ihel * ndiagrams + 901];
      gpuDiagram( &graph, &graphExec, &node902, &node901, diagram902, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node903 = graphNodes[ihel * ndiagrams + 902];
      gpuDiagram( &graph, &graphExec, &node903, &node902, diagram903, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node904 = graphNodes[ihel * ndiagrams + 903];
      gpuDiagram( &graph, &graphExec, &node904, &node903, diagram904, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node905 = graphNodes[ihel * ndiagrams + 904];
      gpuDiagram( &graph, &graphExec, &node905, &node904, diagram905, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node906 = graphNodes[ihel * ndiagrams + 905];
      gpuDiagram( &graph, &graphExec, &node906, &node905, diagram906, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node907 = graphNodes[ihel * ndiagrams + 906];
      gpuDiagram( &graph, &graphExec, &node907, &node906, diagram907, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node908 = graphNodes[ihel * ndiagrams + 907];
      gpuDiagram( &graph, &graphExec, &node908, &node907, diagram908, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node909 = graphNodes[ihel * ndiagrams + 908];
      gpuDiagram( &graph, &graphExec, &node909, &node908, diagram909, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node910 = graphNodes[ihel * ndiagrams + 909];
      gpuDiagram( &graph, &graphExec, &node910, &node909, diagram910, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node911 = graphNodes[ihel * ndiagrams + 910];
      gpuDiagram( &graph, &graphExec, &node911, &node910, diagram911, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node912 = graphNodes[ihel * ndiagrams + 911];
      gpuDiagram( &graph, &graphExec, &node912, &node911, diagram912, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node913 = graphNodes[ihel * ndiagrams + 912];
      gpuDiagram( &graph, &graphExec, &node913, &node912, diagram913, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node914 = graphNodes[ihel * ndiagrams + 913];
      gpuDiagram( &graph, &graphExec, &node914, &node913, diagram914, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node915 = graphNodes[ihel * ndiagrams + 914];
      gpuDiagram( &graph, &graphExec, &node915, &node914, diagram915, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node916 = graphNodes[ihel * ndiagrams + 915];
      gpuDiagram( &graph, &graphExec, &node916, &node915, diagram916, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node917 = graphNodes[ihel * ndiagrams + 916];
      gpuDiagram( &graph, &graphExec, &node917, &node916, diagram917, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node918 = graphNodes[ihel * ndiagrams + 917];
      gpuDiagram( &graph, &graphExec, &node918, &node917, diagram918, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node919 = graphNodes[ihel * ndiagrams + 918];
      gpuDiagram( &graph, &graphExec, &node919, &node918, diagram919, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node920 = graphNodes[ihel * ndiagrams + 919];
      gpuDiagram( &graph, &graphExec, &node920, &node919, diagram920, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node921 = graphNodes[ihel * ndiagrams + 920];
      gpuDiagram( &graph, &graphExec, &node921, &node920, diagram921, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node922 = graphNodes[ihel * ndiagrams + 921];
      gpuDiagram( &graph, &graphExec, &node922, &node921, diagram922, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node923 = graphNodes[ihel * ndiagrams + 922];
      gpuDiagram( &graph, &graphExec, &node923, &node922, diagram923, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node924 = graphNodes[ihel * ndiagrams + 923];
      gpuDiagram( &graph, &graphExec, &node924, &node923, diagram924, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node925 = graphNodes[ihel * ndiagrams + 924];
      gpuDiagram( &graph, &graphExec, &node925, &node924, diagram925, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node926 = graphNodes[ihel * ndiagrams + 925];
      gpuDiagram( &graph, &graphExec, &node926, &node925, diagram926, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node927 = graphNodes[ihel * ndiagrams + 926];
      gpuDiagram( &graph, &graphExec, &node927, &node926, diagram927, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node928 = graphNodes[ihel * ndiagrams + 927];
      gpuDiagram( &graph, &graphExec, &node928, &node927, diagram928, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node929 = graphNodes[ihel * ndiagrams + 928];
      gpuDiagram( &graph, &graphExec, &node929, &node928, diagram929, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node930 = graphNodes[ihel * ndiagrams + 929];
      gpuDiagram( &graph, &graphExec, &node930, &node929, diagram930, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node931 = graphNodes[ihel * ndiagrams + 930];
      gpuDiagram( &graph, &graphExec, &node931, &node930, diagram931, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node932 = graphNodes[ihel * ndiagrams + 931];
      gpuDiagram( &graph, &graphExec, &node932, &node931, diagram932, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node933 = graphNodes[ihel * ndiagrams + 932];
      gpuDiagram( &graph, &graphExec, &node933, &node932, diagram933, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node934 = graphNodes[ihel * ndiagrams + 933];
      gpuDiagram( &graph, &graphExec, &node934, &node933, diagram934, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node935 = graphNodes[ihel * ndiagrams + 934];
      gpuDiagram( &graph, &graphExec, &node935, &node934, diagram935, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node936 = graphNodes[ihel * ndiagrams + 935];
      gpuDiagram( &graph, &graphExec, &node936, &node935, diagram936, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node937 = graphNodes[ihel * ndiagrams + 936];
      gpuDiagram( &graph, &graphExec, &node937, &node936, diagram937, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node938 = graphNodes[ihel * ndiagrams + 937];
      gpuDiagram( &graph, &graphExec, &node938, &node937, diagram938, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node939 = graphNodes[ihel * ndiagrams + 938];
      gpuDiagram( &graph, &graphExec, &node939, &node938, diagram939, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node940 = graphNodes[ihel * ndiagrams + 939];
      gpuDiagram( &graph, &graphExec, &node940, &node939, diagram940, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node941 = graphNodes[ihel * ndiagrams + 940];
      gpuDiagram( &graph, &graphExec, &node941, &node940, diagram941, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node942 = graphNodes[ihel * ndiagrams + 941];
      gpuDiagram( &graph, &graphExec, &node942, &node941, diagram942, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node943 = graphNodes[ihel * ndiagrams + 942];
      gpuDiagram( &graph, &graphExec, &node943, &node942, diagram943, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node944 = graphNodes[ihel * ndiagrams + 943];
      gpuDiagram( &graph, &graphExec, &node944, &node943, diagram944, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node945 = graphNodes[ihel * ndiagrams + 944];
      gpuDiagram( &graph, &graphExec, &node945, &node944, diagram945, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node946 = graphNodes[ihel * ndiagrams + 945];
      gpuDiagram( &graph, &graphExec, &node946, &node945, diagram946, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node947 = graphNodes[ihel * ndiagrams + 946];
      gpuDiagram( &graph, &graphExec, &node947, &node946, diagram947, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node948 = graphNodes[ihel * ndiagrams + 947];
      gpuDiagram( &graph, &graphExec, &node948, &node947, diagram948, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node949 = graphNodes[ihel * ndiagrams + 948];
      gpuDiagram( &graph, &graphExec, &node949, &node948, diagram949, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node950 = graphNodes[ihel * ndiagrams + 949];
      gpuDiagram( &graph, &graphExec, &node950, &node949, diagram950, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node951 = graphNodes[ihel * ndiagrams + 950];
      gpuDiagram( &graph, &graphExec, &node951, &node950, diagram951, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node952 = graphNodes[ihel * ndiagrams + 951];
      gpuDiagram( &graph, &graphExec, &node952, &node951, diagram952, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node953 = graphNodes[ihel * ndiagrams + 952];
      gpuDiagram( &graph, &graphExec, &node953, &node952, diagram953, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node954 = graphNodes[ihel * ndiagrams + 953];
      gpuDiagram( &graph, &graphExec, &node954, &node953, diagram954, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node955 = graphNodes[ihel * ndiagrams + 954];
      gpuDiagram( &graph, &graphExec, &node955, &node954, diagram955, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node956 = graphNodes[ihel * ndiagrams + 955];
      gpuDiagram( &graph, &graphExec, &node956, &node955, diagram956, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node957 = graphNodes[ihel * ndiagrams + 956];
      gpuDiagram( &graph, &graphExec, &node957, &node956, diagram957, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node958 = graphNodes[ihel * ndiagrams + 957];
      gpuDiagram( &graph, &graphExec, &node958, &node957, diagram958, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node959 = graphNodes[ihel * ndiagrams + 958];
      gpuDiagram( &graph, &graphExec, &node959, &node958, diagram959, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node960 = graphNodes[ihel * ndiagrams + 959];
      gpuDiagram( &graph, &graphExec, &node960, &node959, diagram960, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node961 = graphNodes[ihel * ndiagrams + 960];
      gpuDiagram( &graph, &graphExec, &node961, &node960, diagram961, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node962 = graphNodes[ihel * ndiagrams + 961];
      gpuDiagram( &graph, &graphExec, &node962, &node961, diagram962, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node963 = graphNodes[ihel * ndiagrams + 962];
      gpuDiagram( &graph, &graphExec, &node963, &node962, diagram963, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node964 = graphNodes[ihel * ndiagrams + 963];
      gpuDiagram( &graph, &graphExec, &node964, &node963, diagram964, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node965 = graphNodes[ihel * ndiagrams + 964];
      gpuDiagram( &graph, &graphExec, &node965, &node964, diagram965, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node966 = graphNodes[ihel * ndiagrams + 965];
      gpuDiagram( &graph, &graphExec, &node966, &node965, diagram966, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node967 = graphNodes[ihel * ndiagrams + 966];
      gpuDiagram( &graph, &graphExec, &node967, &node966, diagram967, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node968 = graphNodes[ihel * ndiagrams + 967];
      gpuDiagram( &graph, &graphExec, &node968, &node967, diagram968, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node969 = graphNodes[ihel * ndiagrams + 968];
      gpuDiagram( &graph, &graphExec, &node969, &node968, diagram969, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node970 = graphNodes[ihel * ndiagrams + 969];
      gpuDiagram( &graph, &graphExec, &node970, &node969, diagram970, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node971 = graphNodes[ihel * ndiagrams + 970];
      gpuDiagram( &graph, &graphExec, &node971, &node970, diagram971, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node972 = graphNodes[ihel * ndiagrams + 971];
      gpuDiagram( &graph, &graphExec, &node972, &node971, diagram972, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node973 = graphNodes[ihel * ndiagrams + 972];
      gpuDiagram( &graph, &graphExec, &node973, &node972, diagram973, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node974 = graphNodes[ihel * ndiagrams + 973];
      gpuDiagram( &graph, &graphExec, &node974, &node973, diagram974, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node975 = graphNodes[ihel * ndiagrams + 974];
      gpuDiagram( &graph, &graphExec, &node975, &node974, diagram975, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node976 = graphNodes[ihel * ndiagrams + 975];
      gpuDiagram( &graph, &graphExec, &node976, &node975, diagram976, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node977 = graphNodes[ihel * ndiagrams + 976];
      gpuDiagram( &graph, &graphExec, &node977, &node976, diagram977, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node978 = graphNodes[ihel * ndiagrams + 977];
      gpuDiagram( &graph, &graphExec, &node978, &node977, diagram978, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node979 = graphNodes[ihel * ndiagrams + 978];
      gpuDiagram( &graph, &graphExec, &node979, &node978, diagram979, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node980 = graphNodes[ihel * ndiagrams + 979];
      gpuDiagram( &graph, &graphExec, &node980, &node979, diagram980, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node981 = graphNodes[ihel * ndiagrams + 980];
      gpuDiagram( &graph, &graphExec, &node981, &node980, diagram981, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node982 = graphNodes[ihel * ndiagrams + 981];
      gpuDiagram( &graph, &graphExec, &node982, &node981, diagram982, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node983 = graphNodes[ihel * ndiagrams + 982];
      gpuDiagram( &graph, &graphExec, &node983, &node982, diagram983, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node984 = graphNodes[ihel * ndiagrams + 983];
      gpuDiagram( &graph, &graphExec, &node984, &node983, diagram984, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node985 = graphNodes[ihel * ndiagrams + 984];
      gpuDiagram( &graph, &graphExec, &node985, &node984, diagram985, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node986 = graphNodes[ihel * ndiagrams + 985];
      gpuDiagram( &graph, &graphExec, &node986, &node985, diagram986, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node987 = graphNodes[ihel * ndiagrams + 986];
      gpuDiagram( &graph, &graphExec, &node987, &node986, diagram987, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node988 = graphNodes[ihel * ndiagrams + 987];
      gpuDiagram( &graph, &graphExec, &node988, &node987, diagram988, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node989 = graphNodes[ihel * ndiagrams + 988];
      gpuDiagram( &graph, &graphExec, &node989, &node988, diagram989, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node990 = graphNodes[ihel * ndiagrams + 989];
      gpuDiagram( &graph, &graphExec, &node990, &node989, diagram990, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node991 = graphNodes[ihel * ndiagrams + 990];
      gpuDiagram( &graph, &graphExec, &node991, &node990, diagram991, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node992 = graphNodes[ihel * ndiagrams + 991];
      gpuDiagram( &graph, &graphExec, &node992, &node991, diagram992, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node993 = graphNodes[ihel * ndiagrams + 992];
      gpuDiagram( &graph, &graphExec, &node993, &node992, diagram993, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node994 = graphNodes[ihel * ndiagrams + 993];
      gpuDiagram( &graph, &graphExec, &node994, &node993, diagram994, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node995 = graphNodes[ihel * ndiagrams + 994];
      gpuDiagram( &graph, &graphExec, &node995, &node994, diagram995, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node996 = graphNodes[ihel * ndiagrams + 995];
      gpuDiagram( &graph, &graphExec, &node996, &node995, diagram996, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node997 = graphNodes[ihel * ndiagrams + 996];
      gpuDiagram( &graph, &graphExec, &node997, &node996, diagram997, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node998 = graphNodes[ihel * ndiagrams + 997];
      gpuDiagram( &graph, &graphExec, &node998, &node997, diagram998, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node999 = graphNodes[ihel * ndiagrams + 998];
      gpuDiagram( &graph, &graphExec, &node999, &node998, diagram999, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1000 = graphNodes[ihel * ndiagrams + 999];
      gpuDiagram( &graph, &graphExec, &node1000, &node999, diagram1000, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1001 = graphNodes[ihel * ndiagrams + 1000];
      gpuDiagram( &graph, &graphExec, &node1001, &node1000, diagram1001, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1002 = graphNodes[ihel * ndiagrams + 1001];
      gpuDiagram( &graph, &graphExec, &node1002, &node1001, diagram1002, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1003 = graphNodes[ihel * ndiagrams + 1002];
      gpuDiagram( &graph, &graphExec, &node1003, &node1002, diagram1003, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1004 = graphNodes[ihel * ndiagrams + 1003];
      gpuDiagram( &graph, &graphExec, &node1004, &node1003, diagram1004, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1005 = graphNodes[ihel * ndiagrams + 1004];
      gpuDiagram( &graph, &graphExec, &node1005, &node1004, diagram1005, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1006 = graphNodes[ihel * ndiagrams + 1005];
      gpuDiagram( &graph, &graphExec, &node1006, &node1005, diagram1006, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1007 = graphNodes[ihel * ndiagrams + 1006];
      gpuDiagram( &graph, &graphExec, &node1007, &node1006, diagram1007, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1008 = graphNodes[ihel * ndiagrams + 1007];
      gpuDiagram( &graph, &graphExec, &node1008, &node1007, diagram1008, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1009 = graphNodes[ihel * ndiagrams + 1008];
      gpuDiagram( &graph, &graphExec, &node1009, &node1008, diagram1009, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1010 = graphNodes[ihel * ndiagrams + 1009];
      gpuDiagram( &graph, &graphExec, &node1010, &node1009, diagram1010, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1011 = graphNodes[ihel * ndiagrams + 1010];
      gpuDiagram( &graph, &graphExec, &node1011, &node1010, diagram1011, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1012 = graphNodes[ihel * ndiagrams + 1011];
      gpuDiagram( &graph, &graphExec, &node1012, &node1011, diagram1012, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1013 = graphNodes[ihel * ndiagrams + 1012];
      gpuDiagram( &graph, &graphExec, &node1013, &node1012, diagram1013, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1014 = graphNodes[ihel * ndiagrams + 1013];
      gpuDiagram( &graph, &graphExec, &node1014, &node1013, diagram1014, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1015 = graphNodes[ihel * ndiagrams + 1014];
      gpuDiagram( &graph, &graphExec, &node1015, &node1014, diagram1015, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1016 = graphNodes[ihel * ndiagrams + 1015];
      gpuDiagram( &graph, &graphExec, &node1016, &node1015, diagram1016, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1017 = graphNodes[ihel * ndiagrams + 1016];
      gpuDiagram( &graph, &graphExec, &node1017, &node1016, diagram1017, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1018 = graphNodes[ihel * ndiagrams + 1017];
      gpuDiagram( &graph, &graphExec, &node1018, &node1017, diagram1018, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1019 = graphNodes[ihel * ndiagrams + 1018];
      gpuDiagram( &graph, &graphExec, &node1019, &node1018, diagram1019, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1020 = graphNodes[ihel * ndiagrams + 1019];
      gpuDiagram( &graph, &graphExec, &node1020, &node1019, diagram1020, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1021 = graphNodes[ihel * ndiagrams + 1020];
      gpuDiagram( &graph, &graphExec, &node1021, &node1020, diagram1021, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1022 = graphNodes[ihel * ndiagrams + 1021];
      gpuDiagram( &graph, &graphExec, &node1022, &node1021, diagram1022, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1023 = graphNodes[ihel * ndiagrams + 1022];
      gpuDiagram( &graph, &graphExec, &node1023, &node1022, diagram1023, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1024 = graphNodes[ihel * ndiagrams + 1023];
      gpuDiagram( &graph, &graphExec, &node1024, &node1023, diagram1024, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1025 = graphNodes[ihel * ndiagrams + 1024];
      gpuDiagram( &graph, &graphExec, &node1025, &node1024, diagram1025, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1026 = graphNodes[ihel * ndiagrams + 1025];
      gpuDiagram( &graph, &graphExec, &node1026, &node1025, diagram1026, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1027 = graphNodes[ihel * ndiagrams + 1026];
      gpuDiagram( &graph, &graphExec, &node1027, &node1026, diagram1027, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1028 = graphNodes[ihel * ndiagrams + 1027];
      gpuDiagram( &graph, &graphExec, &node1028, &node1027, diagram1028, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1029 = graphNodes[ihel * ndiagrams + 1028];
      gpuDiagram( &graph, &graphExec, &node1029, &node1028, diagram1029, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1030 = graphNodes[ihel * ndiagrams + 1029];
      gpuDiagram( &graph, &graphExec, &node1030, &node1029, diagram1030, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1031 = graphNodes[ihel * ndiagrams + 1030];
      gpuDiagram( &graph, &graphExec, &node1031, &node1030, diagram1031, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1032 = graphNodes[ihel * ndiagrams + 1031];
      gpuDiagram( &graph, &graphExec, &node1032, &node1031, diagram1032, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1033 = graphNodes[ihel * ndiagrams + 1032];
      gpuDiagram( &graph, &graphExec, &node1033, &node1032, diagram1033, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1034 = graphNodes[ihel * ndiagrams + 1033];
      gpuDiagram( &graph, &graphExec, &node1034, &node1033, diagram1034, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1035 = graphNodes[ihel * ndiagrams + 1034];
      gpuDiagram( &graph, &graphExec, &node1035, &node1034, diagram1035, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1036 = graphNodes[ihel * ndiagrams + 1035];
      gpuDiagram( &graph, &graphExec, &node1036, &node1035, diagram1036, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1037 = graphNodes[ihel * ndiagrams + 1036];
      gpuDiagram( &graph, &graphExec, &node1037, &node1036, diagram1037, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1038 = graphNodes[ihel * ndiagrams + 1037];
      gpuDiagram( &graph, &graphExec, &node1038, &node1037, diagram1038, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1039 = graphNodes[ihel * ndiagrams + 1038];
      gpuDiagram( &graph, &graphExec, &node1039, &node1038, diagram1039, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1040 = graphNodes[ihel * ndiagrams + 1039];
      gpuDiagram( &graph, &graphExec, &node1040, &node1039, diagram1040, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1041 = graphNodes[ihel * ndiagrams + 1040];
      gpuDiagram( &graph, &graphExec, &node1041, &node1040, diagram1041, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1042 = graphNodes[ihel * ndiagrams + 1041];
      gpuDiagram( &graph, &graphExec, &node1042, &node1041, diagram1042, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1043 = graphNodes[ihel * ndiagrams + 1042];
      gpuDiagram( &graph, &graphExec, &node1043, &node1042, diagram1043, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1044 = graphNodes[ihel * ndiagrams + 1043];
      gpuDiagram( &graph, &graphExec, &node1044, &node1043, diagram1044, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1045 = graphNodes[ihel * ndiagrams + 1044];
      gpuDiagram( &graph, &graphExec, &node1045, &node1044, diagram1045, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1046 = graphNodes[ihel * ndiagrams + 1045];
      gpuDiagram( &graph, &graphExec, &node1046, &node1045, diagram1046, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1047 = graphNodes[ihel * ndiagrams + 1046];
      gpuDiagram( &graph, &graphExec, &node1047, &node1046, diagram1047, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1048 = graphNodes[ihel * ndiagrams + 1047];
      gpuDiagram( &graph, &graphExec, &node1048, &node1047, diagram1048, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1049 = graphNodes[ihel * ndiagrams + 1048];
      gpuDiagram( &graph, &graphExec, &node1049, &node1048, diagram1049, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1050 = graphNodes[ihel * ndiagrams + 1049];
      gpuDiagram( &graph, &graphExec, &node1050, &node1049, diagram1050, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1051 = graphNodes[ihel * ndiagrams + 1050];
      gpuDiagram( &graph, &graphExec, &node1051, &node1050, diagram1051, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1052 = graphNodes[ihel * ndiagrams + 1051];
      gpuDiagram( &graph, &graphExec, &node1052, &node1051, diagram1052, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1053 = graphNodes[ihel * ndiagrams + 1052];
      gpuDiagram( &graph, &graphExec, &node1053, &node1052, diagram1053, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1054 = graphNodes[ihel * ndiagrams + 1053];
      gpuDiagram( &graph, &graphExec, &node1054, &node1053, diagram1054, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1055 = graphNodes[ihel * ndiagrams + 1054];
      gpuDiagram( &graph, &graphExec, &node1055, &node1054, diagram1055, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1056 = graphNodes[ihel * ndiagrams + 1055];
      gpuDiagram( &graph, &graphExec, &node1056, &node1055, diagram1056, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1057 = graphNodes[ihel * ndiagrams + 1056];
      gpuDiagram( &graph, &graphExec, &node1057, &node1056, diagram1057, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1058 = graphNodes[ihel * ndiagrams + 1057];
      gpuDiagram( &graph, &graphExec, &node1058, &node1057, diagram1058, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1059 = graphNodes[ihel * ndiagrams + 1058];
      gpuDiagram( &graph, &graphExec, &node1059, &node1058, diagram1059, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1060 = graphNodes[ihel * ndiagrams + 1059];
      gpuDiagram( &graph, &graphExec, &node1060, &node1059, diagram1060, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1061 = graphNodes[ihel * ndiagrams + 1060];
      gpuDiagram( &graph, &graphExec, &node1061, &node1060, diagram1061, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1062 = graphNodes[ihel * ndiagrams + 1061];
      gpuDiagram( &graph, &graphExec, &node1062, &node1061, diagram1062, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1063 = graphNodes[ihel * ndiagrams + 1062];
      gpuDiagram( &graph, &graphExec, &node1063, &node1062, diagram1063, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1064 = graphNodes[ihel * ndiagrams + 1063];
      gpuDiagram( &graph, &graphExec, &node1064, &node1063, diagram1064, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1065 = graphNodes[ihel * ndiagrams + 1064];
      gpuDiagram( &graph, &graphExec, &node1065, &node1064, diagram1065, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1066 = graphNodes[ihel * ndiagrams + 1065];
      gpuDiagram( &graph, &graphExec, &node1066, &node1065, diagram1066, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1067 = graphNodes[ihel * ndiagrams + 1066];
      gpuDiagram( &graph, &graphExec, &node1067, &node1066, diagram1067, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1068 = graphNodes[ihel * ndiagrams + 1067];
      gpuDiagram( &graph, &graphExec, &node1068, &node1067, diagram1068, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1069 = graphNodes[ihel * ndiagrams + 1068];
      gpuDiagram( &graph, &graphExec, &node1069, &node1068, diagram1069, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1070 = graphNodes[ihel * ndiagrams + 1069];
      gpuDiagram( &graph, &graphExec, &node1070, &node1069, diagram1070, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1071 = graphNodes[ihel * ndiagrams + 1070];
      gpuDiagram( &graph, &graphExec, &node1071, &node1070, diagram1071, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1072 = graphNodes[ihel * ndiagrams + 1071];
      gpuDiagram( &graph, &graphExec, &node1072, &node1071, diagram1072, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1073 = graphNodes[ihel * ndiagrams + 1072];
      gpuDiagram( &graph, &graphExec, &node1073, &node1072, diagram1073, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1074 = graphNodes[ihel * ndiagrams + 1073];
      gpuDiagram( &graph, &graphExec, &node1074, &node1073, diagram1074, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1075 = graphNodes[ihel * ndiagrams + 1074];
      gpuDiagram( &graph, &graphExec, &node1075, &node1074, diagram1075, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1076 = graphNodes[ihel * ndiagrams + 1075];
      gpuDiagram( &graph, &graphExec, &node1076, &node1075, diagram1076, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1077 = graphNodes[ihel * ndiagrams + 1076];
      gpuDiagram( &graph, &graphExec, &node1077, &node1076, diagram1077, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1078 = graphNodes[ihel * ndiagrams + 1077];
      gpuDiagram( &graph, &graphExec, &node1078, &node1077, diagram1078, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1079 = graphNodes[ihel * ndiagrams + 1078];
      gpuDiagram( &graph, &graphExec, &node1079, &node1078, diagram1079, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1080 = graphNodes[ihel * ndiagrams + 1079];
      gpuDiagram( &graph, &graphExec, &node1080, &node1079, diagram1080, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1081 = graphNodes[ihel * ndiagrams + 1080];
      gpuDiagram( &graph, &graphExec, &node1081, &node1080, diagram1081, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1082 = graphNodes[ihel * ndiagrams + 1081];
      gpuDiagram( &graph, &graphExec, &node1082, &node1081, diagram1082, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1083 = graphNodes[ihel * ndiagrams + 1082];
      gpuDiagram( &graph, &graphExec, &node1083, &node1082, diagram1083, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1084 = graphNodes[ihel * ndiagrams + 1083];
      gpuDiagram( &graph, &graphExec, &node1084, &node1083, diagram1084, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1085 = graphNodes[ihel * ndiagrams + 1084];
      gpuDiagram( &graph, &graphExec, &node1085, &node1084, diagram1085, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1086 = graphNodes[ihel * ndiagrams + 1085];
      gpuDiagram( &graph, &graphExec, &node1086, &node1085, diagram1086, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1087 = graphNodes[ihel * ndiagrams + 1086];
      gpuDiagram( &graph, &graphExec, &node1087, &node1086, diagram1087, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1088 = graphNodes[ihel * ndiagrams + 1087];
      gpuDiagram( &graph, &graphExec, &node1088, &node1087, diagram1088, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1089 = graphNodes[ihel * ndiagrams + 1088];
      gpuDiagram( &graph, &graphExec, &node1089, &node1088, diagram1089, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1090 = graphNodes[ihel * ndiagrams + 1089];
      gpuDiagram( &graph, &graphExec, &node1090, &node1089, diagram1090, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1091 = graphNodes[ihel * ndiagrams + 1090];
      gpuDiagram( &graph, &graphExec, &node1091, &node1090, diagram1091, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1092 = graphNodes[ihel * ndiagrams + 1091];
      gpuDiagram( &graph, &graphExec, &node1092, &node1091, diagram1092, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1093 = graphNodes[ihel * ndiagrams + 1092];
      gpuDiagram( &graph, &graphExec, &node1093, &node1092, diagram1093, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1094 = graphNodes[ihel * ndiagrams + 1093];
      gpuDiagram( &graph, &graphExec, &node1094, &node1093, diagram1094, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1095 = graphNodes[ihel * ndiagrams + 1094];
      gpuDiagram( &graph, &graphExec, &node1095, &node1094, diagram1095, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1096 = graphNodes[ihel * ndiagrams + 1095];
      gpuDiagram( &graph, &graphExec, &node1096, &node1095, diagram1096, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1097 = graphNodes[ihel * ndiagrams + 1096];
      gpuDiagram( &graph, &graphExec, &node1097, &node1096, diagram1097, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1098 = graphNodes[ihel * ndiagrams + 1097];
      gpuDiagram( &graph, &graphExec, &node1098, &node1097, diagram1098, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1099 = graphNodes[ihel * ndiagrams + 1098];
      gpuDiagram( &graph, &graphExec, &node1099, &node1098, diagram1099, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1100 = graphNodes[ihel * ndiagrams + 1099];
      gpuDiagram( &graph, &graphExec, &node1100, &node1099, diagram1100, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1101 = graphNodes[ihel * ndiagrams + 1100];
      gpuDiagram( &graph, &graphExec, &node1101, &node1100, diagram1101, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1102 = graphNodes[ihel * ndiagrams + 1101];
      gpuDiagram( &graph, &graphExec, &node1102, &node1101, diagram1102, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1103 = graphNodes[ihel * ndiagrams + 1102];
      gpuDiagram( &graph, &graphExec, &node1103, &node1102, diagram1103, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1104 = graphNodes[ihel * ndiagrams + 1103];
      gpuDiagram( &graph, &graphExec, &node1104, &node1103, diagram1104, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1105 = graphNodes[ihel * ndiagrams + 1104];
      gpuDiagram( &graph, &graphExec, &node1105, &node1104, diagram1105, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1106 = graphNodes[ihel * ndiagrams + 1105];
      gpuDiagram( &graph, &graphExec, &node1106, &node1105, diagram1106, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1107 = graphNodes[ihel * ndiagrams + 1106];
      gpuDiagram( &graph, &graphExec, &node1107, &node1106, diagram1107, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1108 = graphNodes[ihel * ndiagrams + 1107];
      gpuDiagram( &graph, &graphExec, &node1108, &node1107, diagram1108, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1109 = graphNodes[ihel * ndiagrams + 1108];
      gpuDiagram( &graph, &graphExec, &node1109, &node1108, diagram1109, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1110 = graphNodes[ihel * ndiagrams + 1109];
      gpuDiagram( &graph, &graphExec, &node1110, &node1109, diagram1110, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1111 = graphNodes[ihel * ndiagrams + 1110];
      gpuDiagram( &graph, &graphExec, &node1111, &node1110, diagram1111, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1112 = graphNodes[ihel * ndiagrams + 1111];
      gpuDiagram( &graph, &graphExec, &node1112, &node1111, diagram1112, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1113 = graphNodes[ihel * ndiagrams + 1112];
      gpuDiagram( &graph, &graphExec, &node1113, &node1112, diagram1113, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1114 = graphNodes[ihel * ndiagrams + 1113];
      gpuDiagram( &graph, &graphExec, &node1114, &node1113, diagram1114, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1115 = graphNodes[ihel * ndiagrams + 1114];
      gpuDiagram( &graph, &graphExec, &node1115, &node1114, diagram1115, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1116 = graphNodes[ihel * ndiagrams + 1115];
      gpuDiagram( &graph, &graphExec, &node1116, &node1115, diagram1116, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1117 = graphNodes[ihel * ndiagrams + 1116];
      gpuDiagram( &graph, &graphExec, &node1117, &node1116, diagram1117, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1118 = graphNodes[ihel * ndiagrams + 1117];
      gpuDiagram( &graph, &graphExec, &node1118, &node1117, diagram1118, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1119 = graphNodes[ihel * ndiagrams + 1118];
      gpuDiagram( &graph, &graphExec, &node1119, &node1118, diagram1119, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1120 = graphNodes[ihel * ndiagrams + 1119];
      gpuDiagram( &graph, &graphExec, &node1120, &node1119, diagram1120, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1121 = graphNodes[ihel * ndiagrams + 1120];
      gpuDiagram( &graph, &graphExec, &node1121, &node1120, diagram1121, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1122 = graphNodes[ihel * ndiagrams + 1121];
      gpuDiagram( &graph, &graphExec, &node1122, &node1121, diagram1122, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1123 = graphNodes[ihel * ndiagrams + 1122];
      gpuDiagram( &graph, &graphExec, &node1123, &node1122, diagram1123, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1124 = graphNodes[ihel * ndiagrams + 1123];
      gpuDiagram( &graph, &graphExec, &node1124, &node1123, diagram1124, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1125 = graphNodes[ihel * ndiagrams + 1124];
      gpuDiagram( &graph, &graphExec, &node1125, &node1124, diagram1125, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1126 = graphNodes[ihel * ndiagrams + 1125];
      gpuDiagram( &graph, &graphExec, &node1126, &node1125, diagram1126, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1127 = graphNodes[ihel * ndiagrams + 1126];
      gpuDiagram( &graph, &graphExec, &node1127, &node1126, diagram1127, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1128 = graphNodes[ihel * ndiagrams + 1127];
      gpuDiagram( &graph, &graphExec, &node1128, &node1127, diagram1128, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1129 = graphNodes[ihel * ndiagrams + 1128];
      gpuDiagram( &graph, &graphExec, &node1129, &node1128, diagram1129, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1130 = graphNodes[ihel * ndiagrams + 1129];
      gpuDiagram( &graph, &graphExec, &node1130, &node1129, diagram1130, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1131 = graphNodes[ihel * ndiagrams + 1130];
      gpuDiagram( &graph, &graphExec, &node1131, &node1130, diagram1131, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1132 = graphNodes[ihel * ndiagrams + 1131];
      gpuDiagram( &graph, &graphExec, &node1132, &node1131, diagram1132, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1133 = graphNodes[ihel * ndiagrams + 1132];
      gpuDiagram( &graph, &graphExec, &node1133, &node1132, diagram1133, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1134 = graphNodes[ihel * ndiagrams + 1133];
      gpuDiagram( &graph, &graphExec, &node1134, &node1133, diagram1134, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1135 = graphNodes[ihel * ndiagrams + 1134];
      gpuDiagram( &graph, &graphExec, &node1135, &node1134, diagram1135, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1136 = graphNodes[ihel * ndiagrams + 1135];
      gpuDiagram( &graph, &graphExec, &node1136, &node1135, diagram1136, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1137 = graphNodes[ihel * ndiagrams + 1136];
      gpuDiagram( &graph, &graphExec, &node1137, &node1136, diagram1137, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1138 = graphNodes[ihel * ndiagrams + 1137];
      gpuDiagram( &graph, &graphExec, &node1138, &node1137, diagram1138, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1139 = graphNodes[ihel * ndiagrams + 1138];
      gpuDiagram( &graph, &graphExec, &node1139, &node1138, diagram1139, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1140 = graphNodes[ihel * ndiagrams + 1139];
      gpuDiagram( &graph, &graphExec, &node1140, &node1139, diagram1140, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1141 = graphNodes[ihel * ndiagrams + 1140];
      gpuDiagram( &graph, &graphExec, &node1141, &node1140, diagram1141, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1142 = graphNodes[ihel * ndiagrams + 1141];
      gpuDiagram( &graph, &graphExec, &node1142, &node1141, diagram1142, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1143 = graphNodes[ihel * ndiagrams + 1142];
      gpuDiagram( &graph, &graphExec, &node1143, &node1142, diagram1143, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1144 = graphNodes[ihel * ndiagrams + 1143];
      gpuDiagram( &graph, &graphExec, &node1144, &node1143, diagram1144, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1145 = graphNodes[ihel * ndiagrams + 1144];
      gpuDiagram( &graph, &graphExec, &node1145, &node1144, diagram1145, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1146 = graphNodes[ihel * ndiagrams + 1145];
      gpuDiagram( &graph, &graphExec, &node1146, &node1145, diagram1146, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1147 = graphNodes[ihel * ndiagrams + 1146];
      gpuDiagram( &graph, &graphExec, &node1147, &node1146, diagram1147, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1148 = graphNodes[ihel * ndiagrams + 1147];
      gpuDiagram( &graph, &graphExec, &node1148, &node1147, diagram1148, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1149 = graphNodes[ihel * ndiagrams + 1148];
      gpuDiagram( &graph, &graphExec, &node1149, &node1148, diagram1149, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1150 = graphNodes[ihel * ndiagrams + 1149];
      gpuDiagram( &graph, &graphExec, &node1150, &node1149, diagram1150, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1151 = graphNodes[ihel * ndiagrams + 1150];
      gpuDiagram( &graph, &graphExec, &node1151, &node1150, diagram1151, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1152 = graphNodes[ihel * ndiagrams + 1151];
      gpuDiagram( &graph, &graphExec, &node1152, &node1151, diagram1152, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1153 = graphNodes[ihel * ndiagrams + 1152];
      gpuDiagram( &graph, &graphExec, &node1153, &node1152, diagram1153, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1154 = graphNodes[ihel * ndiagrams + 1153];
      gpuDiagram( &graph, &graphExec, &node1154, &node1153, diagram1154, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1155 = graphNodes[ihel * ndiagrams + 1154];
      gpuDiagram( &graph, &graphExec, &node1155, &node1154, diagram1155, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1156 = graphNodes[ihel * ndiagrams + 1155];
      gpuDiagram( &graph, &graphExec, &node1156, &node1155, diagram1156, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1157 = graphNodes[ihel * ndiagrams + 1156];
      gpuDiagram( &graph, &graphExec, &node1157, &node1156, diagram1157, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1158 = graphNodes[ihel * ndiagrams + 1157];
      gpuDiagram( &graph, &graphExec, &node1158, &node1157, diagram1158, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1159 = graphNodes[ihel * ndiagrams + 1158];
      gpuDiagram( &graph, &graphExec, &node1159, &node1158, diagram1159, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1160 = graphNodes[ihel * ndiagrams + 1159];
      gpuDiagram( &graph, &graphExec, &node1160, &node1159, diagram1160, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1161 = graphNodes[ihel * ndiagrams + 1160];
      gpuDiagram( &graph, &graphExec, &node1161, &node1160, diagram1161, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1162 = graphNodes[ihel * ndiagrams + 1161];
      gpuDiagram( &graph, &graphExec, &node1162, &node1161, diagram1162, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1163 = graphNodes[ihel * ndiagrams + 1162];
      gpuDiagram( &graph, &graphExec, &node1163, &node1162, diagram1163, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1164 = graphNodes[ihel * ndiagrams + 1163];
      gpuDiagram( &graph, &graphExec, &node1164, &node1163, diagram1164, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1165 = graphNodes[ihel * ndiagrams + 1164];
      gpuDiagram( &graph, &graphExec, &node1165, &node1164, diagram1165, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1166 = graphNodes[ihel * ndiagrams + 1165];
      gpuDiagram( &graph, &graphExec, &node1166, &node1165, diagram1166, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1167 = graphNodes[ihel * ndiagrams + 1166];
      gpuDiagram( &graph, &graphExec, &node1167, &node1166, diagram1167, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1168 = graphNodes[ihel * ndiagrams + 1167];
      gpuDiagram( &graph, &graphExec, &node1168, &node1167, diagram1168, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1169 = graphNodes[ihel * ndiagrams + 1168];
      gpuDiagram( &graph, &graphExec, &node1169, &node1168, diagram1169, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1170 = graphNodes[ihel * ndiagrams + 1169];
      gpuDiagram( &graph, &graphExec, &node1170, &node1169, diagram1170, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1171 = graphNodes[ihel * ndiagrams + 1170];
      gpuDiagram( &graph, &graphExec, &node1171, &node1170, diagram1171, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1172 = graphNodes[ihel * ndiagrams + 1171];
      gpuDiagram( &graph, &graphExec, &node1172, &node1171, diagram1172, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1173 = graphNodes[ihel * ndiagrams + 1172];
      gpuDiagram( &graph, &graphExec, &node1173, &node1172, diagram1173, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1174 = graphNodes[ihel * ndiagrams + 1173];
      gpuDiagram( &graph, &graphExec, &node1174, &node1173, diagram1174, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1175 = graphNodes[ihel * ndiagrams + 1174];
      gpuDiagram( &graph, &graphExec, &node1175, &node1174, diagram1175, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1176 = graphNodes[ihel * ndiagrams + 1175];
      gpuDiagram( &graph, &graphExec, &node1176, &node1175, diagram1176, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1177 = graphNodes[ihel * ndiagrams + 1176];
      gpuDiagram( &graph, &graphExec, &node1177, &node1176, diagram1177, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1178 = graphNodes[ihel * ndiagrams + 1177];
      gpuDiagram( &graph, &graphExec, &node1178, &node1177, diagram1178, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1179 = graphNodes[ihel * ndiagrams + 1178];
      gpuDiagram( &graph, &graphExec, &node1179, &node1178, diagram1179, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1180 = graphNodes[ihel * ndiagrams + 1179];
      gpuDiagram( &graph, &graphExec, &node1180, &node1179, diagram1180, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1181 = graphNodes[ihel * ndiagrams + 1180];
      gpuDiagram( &graph, &graphExec, &node1181, &node1180, diagram1181, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1182 = graphNodes[ihel * ndiagrams + 1181];
      gpuDiagram( &graph, &graphExec, &node1182, &node1181, diagram1182, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1183 = graphNodes[ihel * ndiagrams + 1182];
      gpuDiagram( &graph, &graphExec, &node1183, &node1182, diagram1183, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1184 = graphNodes[ihel * ndiagrams + 1183];
      gpuDiagram( &graph, &graphExec, &node1184, &node1183, diagram1184, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1185 = graphNodes[ihel * ndiagrams + 1184];
      gpuDiagram( &graph, &graphExec, &node1185, &node1184, diagram1185, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1186 = graphNodes[ihel * ndiagrams + 1185];
      gpuDiagram( &graph, &graphExec, &node1186, &node1185, diagram1186, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1187 = graphNodes[ihel * ndiagrams + 1186];
      gpuDiagram( &graph, &graphExec, &node1187, &node1186, diagram1187, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1188 = graphNodes[ihel * ndiagrams + 1187];
      gpuDiagram( &graph, &graphExec, &node1188, &node1187, diagram1188, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1189 = graphNodes[ihel * ndiagrams + 1188];
      gpuDiagram( &graph, &graphExec, &node1189, &node1188, diagram1189, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1190 = graphNodes[ihel * ndiagrams + 1189];
      gpuDiagram( &graph, &graphExec, &node1190, &node1189, diagram1190, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1191 = graphNodes[ihel * ndiagrams + 1190];
      gpuDiagram( &graph, &graphExec, &node1191, &node1190, diagram1191, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1192 = graphNodes[ihel * ndiagrams + 1191];
      gpuDiagram( &graph, &graphExec, &node1192, &node1191, diagram1192, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1193 = graphNodes[ihel * ndiagrams + 1192];
      gpuDiagram( &graph, &graphExec, &node1193, &node1192, diagram1193, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1194 = graphNodes[ihel * ndiagrams + 1193];
      gpuDiagram( &graph, &graphExec, &node1194, &node1193, diagram1194, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1195 = graphNodes[ihel * ndiagrams + 1194];
      gpuDiagram( &graph, &graphExec, &node1195, &node1194, diagram1195, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1196 = graphNodes[ihel * ndiagrams + 1195];
      gpuDiagram( &graph, &graphExec, &node1196, &node1195, diagram1196, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1197 = graphNodes[ihel * ndiagrams + 1196];
      gpuDiagram( &graph, &graphExec, &node1197, &node1196, diagram1197, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1198 = graphNodes[ihel * ndiagrams + 1197];
      gpuDiagram( &graph, &graphExec, &node1198, &node1197, diagram1198, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1199 = graphNodes[ihel * ndiagrams + 1198];
      gpuDiagram( &graph, &graphExec, &node1199, &node1198, diagram1199, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1200 = graphNodes[ihel * ndiagrams + 1199];
      gpuDiagram( &graph, &graphExec, &node1200, &node1199, diagram1200, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1201 = graphNodes[ihel * ndiagrams + 1200];
      gpuDiagram( &graph, &graphExec, &node1201, &node1200, diagram1201, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1202 = graphNodes[ihel * ndiagrams + 1201];
      gpuDiagram( &graph, &graphExec, &node1202, &node1201, diagram1202, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1203 = graphNodes[ihel * ndiagrams + 1202];
      gpuDiagram( &graph, &graphExec, &node1203, &node1202, diagram1203, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1204 = graphNodes[ihel * ndiagrams + 1203];
      gpuDiagram( &graph, &graphExec, &node1204, &node1203, diagram1204, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1205 = graphNodes[ihel * ndiagrams + 1204];
      gpuDiagram( &graph, &graphExec, &node1205, &node1204, diagram1205, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1206 = graphNodes[ihel * ndiagrams + 1205];
      gpuDiagram( &graph, &graphExec, &node1206, &node1205, diagram1206, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1207 = graphNodes[ihel * ndiagrams + 1206];
      gpuDiagram( &graph, &graphExec, &node1207, &node1206, diagram1207, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1208 = graphNodes[ihel * ndiagrams + 1207];
      gpuDiagram( &graph, &graphExec, &node1208, &node1207, diagram1208, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1209 = graphNodes[ihel * ndiagrams + 1208];
      gpuDiagram( &graph, &graphExec, &node1209, &node1208, diagram1209, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1210 = graphNodes[ihel * ndiagrams + 1209];
      gpuDiagram( &graph, &graphExec, &node1210, &node1209, diagram1210, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1211 = graphNodes[ihel * ndiagrams + 1210];
      gpuDiagram( &graph, &graphExec, &node1211, &node1210, diagram1211, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1212 = graphNodes[ihel * ndiagrams + 1211];
      gpuDiagram( &graph, &graphExec, &node1212, &node1211, diagram1212, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1213 = graphNodes[ihel * ndiagrams + 1212];
      gpuDiagram( &graph, &graphExec, &node1213, &node1212, diagram1213, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1214 = graphNodes[ihel * ndiagrams + 1213];
      gpuDiagram( &graph, &graphExec, &node1214, &node1213, diagram1214, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1215 = graphNodes[ihel * ndiagrams + 1214];
      gpuDiagram( &graph, &graphExec, &node1215, &node1214, diagram1215, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1216 = graphNodes[ihel * ndiagrams + 1215];
      gpuDiagram( &graph, &graphExec, &node1216, &node1215, diagram1216, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1217 = graphNodes[ihel * ndiagrams + 1216];
      gpuDiagram( &graph, &graphExec, &node1217, &node1216, diagram1217, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1218 = graphNodes[ihel * ndiagrams + 1217];
      gpuDiagram( &graph, &graphExec, &node1218, &node1217, diagram1218, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1219 = graphNodes[ihel * ndiagrams + 1218];
      gpuDiagram( &graph, &graphExec, &node1219, &node1218, diagram1219, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1220 = graphNodes[ihel * ndiagrams + 1219];
      gpuDiagram( &graph, &graphExec, &node1220, &node1219, diagram1220, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1221 = graphNodes[ihel * ndiagrams + 1220];
      gpuDiagram( &graph, &graphExec, &node1221, &node1220, diagram1221, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1222 = graphNodes[ihel * ndiagrams + 1221];
      gpuDiagram( &graph, &graphExec, &node1222, &node1221, diagram1222, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1223 = graphNodes[ihel * ndiagrams + 1222];
      gpuDiagram( &graph, &graphExec, &node1223, &node1222, diagram1223, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1224 = graphNodes[ihel * ndiagrams + 1223];
      gpuDiagram( &graph, &graphExec, &node1224, &node1223, diagram1224, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1225 = graphNodes[ihel * ndiagrams + 1224];
      gpuDiagram( &graph, &graphExec, &node1225, &node1224, diagram1225, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1226 = graphNodes[ihel * ndiagrams + 1225];
      gpuDiagram( &graph, &graphExec, &node1226, &node1225, diagram1226, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1227 = graphNodes[ihel * ndiagrams + 1226];
      gpuDiagram( &graph, &graphExec, &node1227, &node1226, diagram1227, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1228 = graphNodes[ihel * ndiagrams + 1227];
      gpuDiagram( &graph, &graphExec, &node1228, &node1227, diagram1228, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1229 = graphNodes[ihel * ndiagrams + 1228];
      gpuDiagram( &graph, &graphExec, &node1229, &node1228, diagram1229, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1230 = graphNodes[ihel * ndiagrams + 1229];
      gpuDiagram( &graph, &graphExec, &node1230, &node1229, diagram1230, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1231 = graphNodes[ihel * ndiagrams + 1230];
      gpuDiagram( &graph, &graphExec, &node1231, &node1230, diagram1231, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1232 = graphNodes[ihel * ndiagrams + 1231];
      gpuDiagram( &graph, &graphExec, &node1232, &node1231, diagram1232, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1233 = graphNodes[ihel * ndiagrams + 1232];
      gpuDiagram( &graph, &graphExec, &node1233, &node1232, diagram1233, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1234 = graphNodes[ihel * ndiagrams + 1233];
      gpuDiagram( &graph, &graphExec, &node1234, &node1233, diagram1234, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1235 = graphNodes[ihel * ndiagrams + 1234];
      gpuDiagram( &graph, &graphExec, &node1235, &node1234, diagram1235, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1236 = graphNodes[ihel * ndiagrams + 1235];
      gpuDiagram( &graph, &graphExec, &node1236, &node1235, diagram1236, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1237 = graphNodes[ihel * ndiagrams + 1236];
      gpuDiagram( &graph, &graphExec, &node1237, &node1236, diagram1237, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1238 = graphNodes[ihel * ndiagrams + 1237];
      gpuDiagram( &graph, &graphExec, &node1238, &node1237, diagram1238, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1239 = graphNodes[ihel * ndiagrams + 1238];
      gpuDiagram( &graph, &graphExec, &node1239, &node1238, diagram1239, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
      gpuGraphNode_t& node1240 = graphNodes[ihel * ndiagrams + 1239];
      gpuDiagram( &graph, &graphExec, &node1240, &node1239, diagram1240, gpublocks, gputhreads, gpustream, wfs, jamps, channelIds, couplings, numerators, denominators );
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
