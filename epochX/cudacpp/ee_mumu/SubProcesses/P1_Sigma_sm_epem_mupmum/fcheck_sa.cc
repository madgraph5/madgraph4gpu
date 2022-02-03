#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <string>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "mgOnGpuConfig.h"

#include "Bridge.h"
#include "CPPProcess.h"
//#include "CrossSectionKernels.h"
#include "MatrixElementKernels.h"
#include "MemoryAccessMomenta.h"
#include "MemoryAccessRandomNumbers.h"
#include "MemoryAccessWeights.h"
#include "MemoryBuffers.h"
#include "RamboSamplingKernels.h"
#include "RandomNumberKernels.h"

using FORTRANFPTYPE = double;

void fsampler( FORTRANFPTYPE* fortranMomenta )
{
  // Namespaces for CUDA and C++
#ifdef __CUDACC__
  using namespace mg5amcGpu;
#else
  using namespace mg5amcCpu;
#endif
  
  // HARDCODED DEFAULTS
  constexpr bool verbose = false;
  constexpr int gpublocks = 4;
  constexpr int gputhreads = 32;
  constexpr int niter = 2;
  constexpr int ndim = gpublocks * gputhreads; // number of threads in one GPU grid
  constexpr int nevt = ndim; // number of events in one iteration == number of GPU threads

#ifndef __CUDACC__
#ifdef _OPENMP
  // Set OMP_NUM_THREADS equal to 1 if it is not yet set
  omp_set_num_threads( 1 ); // https://stackoverflow.com/a/22816325
#endif
#endif

#ifndef __CUDACC__
  // Fail gently and avoid "Illegal instruction (core dumped)" if the host does not support the requested AVX
  // [NB: this prevents a crash on pmpe04 but not on some github CI nodes]
  auto supportsAvx = [](){
#if defined __AVX512VL__
    bool ok = __builtin_cpu_supports( "avx512vl" );
    const std::string tag = "skylake-avx512 (AVX512VL)";
#elif defined __AVX2__
    bool ok = __builtin_cpu_supports( "avx2" );
    const std::string tag = "haswell (AVX2)";
#elif defined __SSE4_2__
#ifdef __PPC__
    // See https://gcc.gnu.org/onlinedocs/gcc/Basic-PowerPC-Built-in-Functions-Available-on-all-Configurations.html
    bool ok = __builtin_cpu_supports( "vsx" );
    const std::string tag = "powerpc vsx (128bit as in SSE4.2)";
#else
    bool ok = __builtin_cpu_supports( "sse4.2" );
    const std::string tag = "nehalem (SSE4.2)";
#endif
#else
    bool ok = true;
    const std::string tag = "none";
#endif
    if ( tag == "none" )
      std::cout << "INFO: The application does not require the host to support any AVX feature" << std::endl;
    else if ( ok )
      std::cout << "INFO: The application is built for " << tag << " and the host supports it" << std::endl;
    else
      std::cout << "ERROR! The application is built for " << tag << " but the host does not support it" << std::endl;
    return ok;
  };
  if ( ! supportsAvx() ) throw std::runtime_error( "SIMD mode is not supported" );
#endif

#ifdef __CUDACC__
  // === STEP 0 - INITIALISE
  checkCuda( cudaFree( 0 ) ); // SLOW!

  // --- Book the tear down at the end of main:
  constexpr bool debug = false;
  struct CudaTearDown {
    CudaTearDown(bool print) : _print(print) { }
    ~CudaTearDown() {
      //if ( _print ) std::cout << "Calling cudaDeviceReset()." << std::endl;
      checkCuda( cudaDeviceReset() ); // this is needed by cuda-memcheck --leak-check full
    }
    bool _print{false};
  } cudaTearDown(debug);
#endif

  // --- 0a. Initialise physics process

  // Create a process object
  CPPProcess process( niter, gpublocks, gputhreads, verbose );

  // Read param_card and set parameters
  process.initProc("../../Cards/param_card.dat");
  constexpr fptype energy = 1500; // historical default, Ecms = 1500 GeV = 1.5 TeV (above the Z peak)

  // --- 0b. Allocate memory structures

  // Memory buffers for random numbers
#ifndef __CUDACC__
  HostBufferRandomNumbers hstRnarray( nevt );
#else
  PinnedHostBufferRandomNumbers hstRnarray( nevt );
  DeviceBufferRandomNumbers devRnarray( nevt );
#endif

  // Memory buffers for momenta
#ifndef __CUDACC__
  HostBufferMomenta hstMomenta( nevt );
#else
  PinnedHostBufferMomenta hstMomenta( nevt );
  DeviceBufferMomenta devMomenta( nevt );
#endif

  // Memory buffers for sampling weights
#ifndef __CUDACC__
  HostBufferWeights hstWeights( nevt );
#else
  PinnedHostBufferWeights hstWeights( nevt );
  DeviceBufferWeights devWeights( nevt );
#endif

  // --- 0c. Create curand or common generator

  // Allocate the appropriate RandomNumberKernel
  std::unique_ptr<RandomNumberKernelBase> prnk;
  prnk.reset( new CommonRandomNumberKernel( hstRnarray ) );

  // --- 0c. Create rambo sampling kernel
  std::unique_ptr<SamplingKernelBase> prsk;
  prsk.reset( new RamboSamplingKernelHost( energy, hstRnarray, hstMomenta, hstWeights, nevt ) );

  /*
  // Memory buffers for matrix elements
#ifndef __CUDACC__
  HostBufferMatrixElements hstMatrixElements( nevt );
#else
  PinnedHostBufferMatrixElements hstMatrixElements( nevt );
  DeviceBufferMatrixElements devMatrixElements( nevt );
#endif
  */

  /*
  // --- 0c. Create cross section kernel [keep this in 0c for the moment]
  EventStatistics hstStats;
  CrossSectionKernelHost xsk( hstWeights, hstMatrixElements, hstStats, nevt );
  */

  // **************************************
  // *** START MAIN LOOP ON #ITERATIONS ***
  // **************************************

  for (int iiter = 0; iiter < niter; ++iiter)
  {
    //std::cout << "Iteration #" << iiter+1 << " of " << niter << std::endl;

    // === STEP 1 OF 3

    // --- 1a. Seed curand generator (to get same results on host and device)
    // [NB This should not be necessary using the host API: "Generation functions
    // can be called multiple times on the same generator to generate successive
    // blocks of results. For pseudorandom generators, multiple calls to generation
    // functions will yield the same result as a single call with a large size."]
    const unsigned long long seed = 20200805;
    prnk->seedGenerator( seed+iiter );

    // --- 1b. Generate all relevant numbers to build nevt events (i.e. nevt phase space points) on the host
    prnk->generateRnarray();
    //std::cout << "Got random numbers" << std::endl;

    // === STEP 2 OF 3

    // --- 2a. Fill in momenta of initial state particles on the device
    prsk->getMomentaInitial();
    //std::cout << "Got initial momenta" << std::endl;

    // --- 2b. Fill in momenta of final state particles using the RAMBO algorithm on the device
    // (i.e. map random numbers to final-state particle momenta for each of nevt events)
    prsk->getMomentaFinal();
    //std::cout << "Got final momenta" << std::endl;

    // === STEP 3 OF 3
    // Evaluate matrix elements for all nevt events
    // 0d. For Bridge only, transpose C2F [renamed as 0d: this is not initialisation, but I want it out of the ME timers (#371)]
    // 0e. (Only on the first iteration) Get good helicities [renamed as 0e: this IS initialisation!]
    // 3a. Evaluate MEs on the device (include transpose F2C for Bridge)
    // 3b. Copy MEs back from device to host

    // --- 0d. TransC2F
    hst_transposeMomentaC2F( hstMomenta.data(), fortranMomenta, nevt );
  }  
}
