#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <unistd.h>

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"

#ifdef __CUDACC__
#include "grambo.cu"
#else
#include "rambo.h"
#endif

#include "CPPProcess.h"
#include "timermap.h"

bool is_number(const char *s) {
  const char *t = s;
  while (*t != '\0' && isdigit(*t))
    ++t;
  return (int)strlen(s) == t - s;
}

int usage(char* argv0, int ret = 1) {
  std::cout << "Usage: " << argv0
            << " [--verbose|-v] [--debug|-d] [--performance|-p]"
            << " [#gpuBlocksPerGrid #gpuThreadsPerBlock] #iterations" << std::endl << std::endl;
  std::cout << "The number of events per iteration is #gpuBlocksPerGrid * #gpuThreadsPerBlock" << std::endl;
  std::cout << "(also in CPU/C++ code, where only the product of these two parameters counts)" << std::endl;
  return ret;
}

int main(int argc, char **argv)
{
  // READ COMMAND LINE ARGUMENTS
  bool verbose = false;
  bool debug = false;
  bool perf = false;
  int niter = 0;
  int gpublocks = 1;
  int gputhreads = 32;
  int numvec[3] = {0,0,0};
  int nnum = 0;
  for (int argn = 1; argn < argc; ++argn) {
    if (strcmp(argv[argn], "--verbose") == 0 || strcmp(argv[argn], "-v") == 0)
      verbose = true;
    else if (strcmp(argv[argn], "--debug") == 0 ||
             strcmp(argv[argn], "-d") == 0)
      debug = true;
    else if (strcmp(argv[argn], "--performance") == 0 ||
             strcmp(argv[argn], "-p") == 0)
      perf = true;
    else if (is_number(argv[argn]) && nnum<3)
      numvec[nnum++] = atoi(argv[argn]);
    else
      return usage(argv[0]);
  }

  if (nnum == 3) {
    gpublocks = numvec[0];
    gputhreads = numvec[1];
    niter = numvec[2];
  } else if (nnum == 1) {
    niter = numvec[0];
  } else {
    return usage(argv[0]);
  }

  if (niter == 0)
    return usage(argv[0]);

  const int neppR = mgOnGpu::neppR; // ASA layout: constant at compile-time
  if ( gputhreads%neppR != 0 )
  {
    std::cout << "ERROR! #threads/block should be a multiple of " << neppR << std::endl;
    return usage(argv[0]);
  }

  const int neppM = mgOnGpu::neppM; // ASA layout: constant at compile-time
  if ( gputhreads%neppM != 0 )
  {
    std::cout << "ERROR! #threads/block should be a multiple of " << neppM << std::endl;
    return usage(argv[0]);
  }

  using mgOnGpu::ntpbMAX;
  if ( gputhreads > ntpbMAX )
  {
    std::cout << "ERROR! #threads/block should be <= " << ntpbMAX << std::endl;
    return usage(argv[0]);
  }

  const int ndim = gpublocks * gputhreads; // number of threads in one GPU grid
  const int nevt = ndim; // number of events in one iteration == number of GPU threads
  const int nevtALL = niter*nevt; // total number of ALL events in all iterations

  if (verbose)
    std::cout << "# iterations: " << niter << std::endl;

  // *** START THE NEW TIMERS ***
  mgOnGpu::TimerMap timermap;

  // === STEP 0 - INITIALISE

#ifdef __CUDACC__
  // --- 00. Initialise cuda (call cudaFree to ease cuda profile analysis)
  const std::string cdfrKey = "00 CudaFree";
  timermap.start( cdfrKey );
  //std::cout << "Calling cudaFree... " << std::endl;
  checkCuda( cudaFree( 0 ) ); // SLOW!
  //std::cout << "Calling cudaFree... done" << std::endl;
#endif

  // --- 0a. Initialise physics process
  const std::string procKey = "0a ProcInit";
  timermap.start( procKey );

  // Create a process object
#ifdef __CUDACC__
  gProc::CPPProcess process( niter, gpublocks, gputhreads, verbose );
#else
  Proc::CPPProcess process( niter, gpublocks, gputhreads, verbose );
#endif

  // Read param_card and set parameters
  process.initProc("../../Cards/param_card.dat");

  const fptype energy = 1500;
  const int meGeVexponent = -(2 * process.nexternal - 8);

  // --- 0b. Allocate memory structures
  const std::string alloKey = "0b MemAlloc";
  timermap.start( alloKey );

  // Memory structures for random numbers, momenta, matrix elements and weights on host and device
  using mgOnGpu::np4;
  using mgOnGpu::nparf;
  using mgOnGpu::npar;
  const int nRnarray = np4*nparf*nevt; // (NB: ASA layout with nevt=npagR*neppR events per iteration)
#ifdef __CUDACC__
  const int nbytesRnarray = nRnarray * sizeof(fptype);
  fptype* devRnarray = 0; // AOSOA[npagR][nparf][np4][neppR] (NB: nevt=npagR*neppR)
  checkCuda( cudaMalloc( &devRnarray, nbytesRnarray ) );
#if defined MGONGPU_CURAND_ONHOST
  fptype* hstRnarray = 0; // AOSOA[npagR][nparf][np4][neppR] (NB: nevt=npagR*neppR)
  checkCuda( cudaMallocHost( &hstRnarray, nbytesRnarray ) );
#endif
#else
  fptype* hstRnarray = 0; // AOSOA[npagR][nparf][np4][neppR] (NB: nevt=npagR*neppR)
  hstRnarray = new fptype[nRnarray]();
#endif

  const int nMomenta = np4*npar*nevt; // (NB: nevt=npagM*neppM for ASA layouts)
  fptype* hstMomenta = 0; // AOSOA[npagM][npar][np4][neppM] (previously was: lp)
#ifdef __CUDACC__
  const int nbytesMomenta = nMomenta * sizeof(fptype);
  checkCuda( cudaMallocHost( &hstMomenta, nbytesMomenta ) );
  fptype* devMomenta = 0; // (previously was: allMomenta)
  checkCuda( cudaMalloc( &devMomenta, nbytesMomenta ) );
#else
  hstMomenta = new fptype[nMomenta]();
#endif

#if defined __CUDACC__ && defined MGONGPU_WFMEM_GLOBAL
  using mgOnGpu::nwf;
  using mgOnGpu::nw6;
  const int nAllWFs = nwf * nw6 * nevt;
  const int nbytesAllWFs = nAllWFs * sizeof(cxtype);
  cxtype* devAllWFs = 0; // AOSOA[nblk][nparf][np4][ntpb] (NB: nevt=nblk*ntpb)
  checkCuda( cudaMalloc( &devAllWFs, nbytesAllWFs ) );
#endif

  const int nWeights = nevt;
  fptype* hstWeights = 0; // (previously was: meHostPtr)
#ifdef __CUDACC__
  const int nbytesWeights = nWeights * sizeof(fptype);
  checkCuda( cudaMallocHost( &hstWeights, nbytesWeights ) );
  fptype* devWeights = 0; // (previously was: meDevPtr)
  checkCuda( cudaMalloc( &devWeights, nbytesWeights ) );
#else
  hstWeights = new fptype[nWeights]();
#endif

  const int nMEs = nevt;
  fptype* hstMEs = 0; // (previously was: meHostPtr)
#ifdef __CUDACC__
  const int nbytesMEs = nMEs * sizeof(fptype);
  checkCuda( cudaMallocHost( &hstMEs, nbytesMEs ) );
  fptype* devMEs = 0; // (previously was: meDevPtr)
  checkCuda( cudaMalloc( &devMEs, nbytesMEs ) );
#else
  hstMEs = new fptype[nMEs]();
#endif

  double* rambtimes = new double[niter]();
  double* wavetimes = new double[niter]();
  fptype* matrixelementALL = new fptype[nevtALL](); // FIXME: assume process.nprocesses == 1

#ifdef __CUDACC__
#if defined MGONGPU_WFMEM_SHARED
  const int nbytesSharedSK = gProc::sigmakin_sharedmem_nbytes(gputhreads);
#endif
#endif

  // --- 0c. Create curand generator
  const std::string cgenKey = "0c GenCreat";
  timermap.start( cgenKey );
  curandGenerator_t rnGen;
#ifdef __CUDACC__
  grambo2toNm0::createGenerator( &rnGen );
#else
  rambo2toNm0::createGenerator( &rnGen );
#endif

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
    const std::string sgenKey = "1a GenSeed ";
    timermap.start( sgenKey );
    const unsigned long long seed = 20200805;
#ifdef __CUDACC__
    grambo2toNm0::seedGenerator( rnGen, seed+iiter );
#else
    rambo2toNm0::seedGenerator( rnGen, seed+iiter );
#endif

    // --- 1b. Generate all relevant numbers to build nevt events (i.e. nevt phase space points) on the host
    const std::string rngnKey = "1b GenRnGen";
    timermap.start( rngnKey );
#ifdef __CUDACC__
#if defined MGONGPU_CURAND_ONDEVICE
    grambo2toNm0::generateRnarray( rnGen, devRnarray, nevt );
#elif defined MGONGPU_CURAND_ONHOST
    grambo2toNm0::generateRnarray( rnGen, hstRnarray, nevt );
#endif
#else
    rambo2toNm0::generateRnarray( rnGen, hstRnarray, nevt );
#endif
    //std::cout << "Got random numbers" << std::endl;

#ifdef __CUDACC__
#if defined MGONGPU_CURAND_ONHOST
    // --- 1c. Copy rnarray from host to device
    const std::string htodKey = "1c CpHTDrnd";
    timermap.start( htodKey );
    checkCuda( cudaMemcpy( devRnarray, hstRnarray, nbytesRnarray, cudaMemcpyHostToDevice ) );
#endif
#endif

    // === STEP 2 OF 3
    // Fill in particle momenta for each of nevt events on the device

    // *** START THE OLD-STYLE TIMER FOR RAMBO ***
    double rambtime = 0;

    // --- 2a. Fill in momenta of initial state particles on the device
    const std::string riniKey = "2a RamboIni";
    timermap.start( riniKey );
#ifdef __CUDACC__
    grambo2toNm0::getMomentaInitial<<<gpublocks, gputhreads>>>( energy, devMomenta );
#else
    rambo2toNm0::getMomentaInitial( energy, hstMomenta, nevt );
#endif
    //std::cout << "Got initial momenta" << std::endl;

    // --- 2b. Fill in momenta of final state particles using the RAMBO algorithm on the device
    // (i.e. map random numbers to final-state particle momenta for each of nevt events)
    const std::string rfinKey = "2b RamboFin";
    rambtime += timermap.start( rfinKey );
#ifdef __CUDACC__
    grambo2toNm0::getMomentaFinal<<<gpublocks, gputhreads>>>( energy, devRnarray, devMomenta, devWeights );
#else
    rambo2toNm0::getMomentaFinal( energy, hstRnarray, hstMomenta, hstWeights, nevt );
#endif
    //std::cout << "Got final momenta" << std::endl;

#ifdef __CUDACC__
    // --- 2c. CopyDToH Weights
    const std::string cwgtKey = "2c CpDTHwgt";
    rambtime += timermap.start( cwgtKey );
    checkCuda( cudaMemcpy( hstWeights, devWeights, nbytesWeights, cudaMemcpyDeviceToHost ) );

    // --- 2d. CopyDToH Momenta
    const std::string cmomKey = "2d CpDTHmom";
    rambtime += timermap.start( cmomKey );
    checkCuda( cudaMemcpy( hstMomenta, devMomenta, nbytesMomenta, cudaMemcpyDeviceToHost ) );
#endif

    // *** STOP THE OLD-STYLE TIMER FOR RAMBO ***
    rambtime += timermap.stop();

    // === STEP 3 OF 3
    // Evaluate matrix elements for all nevt events
    // 3a. Evaluate MEs on the device
    // 3b. Copy MEs back from device to host

    // *** START THE OLD TIMER FOR WAVEFUNCTIONS ***
    double wavetime = 0;

    // --- 3a. SigmaKin
    const std::string skinKey = "3a SigmaKin";
    timermap.start( skinKey );
#ifdef __CUDACC__
#if defined MGONGPU_WFMEM_GLOBAL
    gProc::sigmaKin<<<gpublocks, gputhreads>>>(devMomenta, devMEs, devAllWFs);
#elif defined MGONGPU_WFMEM_SHARED
    gProc::sigmaKin<<<gpublocks, gputhreads, nbytesSharedSK>>>(devMomenta, devMEs);
#else
    gProc::sigmaKin<<<gpublocks, gputhreads>>>(devMomenta, devMEs);
#endif
    checkCuda( cudaPeekAtLastError() );
#else
    Proc::sigmaKin(hstMomenta, hstMEs, nevt);
#endif

#ifdef __CUDACC__
    // --- 3b. CopyDToH MEs
    const std::string cmesKey = "3b CpDTHmes";
    wavetime += timermap.start( cmesKey );
    checkCuda( cudaMemcpy( hstMEs, devMEs, nbytesMEs, cudaMemcpyDeviceToHost ) );
#endif

    // *** STOP THE OLD TIMER FOR WAVEFUNCTIONS ***
    wavetime += timermap.stop();

    // === STEP 4 FINALISE LOOP
    // --- 4a Dump within the loop
    const std::string loopKey = "4a DumpLoop";
    timermap.start(loopKey);
    rambtimes[iiter] = rambtime;
    wavetimes[iiter] = wavetime;

    if (verbose)
    {
      std::cout << "***********************************" << std::endl
                << "Iteration #" << iiter+1 << " of " << niter << std::endl;
      if (perf) std::cout << "Wave function time: " << wavetime << std::endl;
    }

    if (verbose || perf)
    {
      for (int ievt = 0; ievt < nevt; ++ievt) // Loop over all events in this iteration
      {
        const int ipagM = ievt/neppM; // #eventpage in this iteration
        const int ieppM = ievt%neppM; // #event in the current eventpage in this iteration
        if (verbose)
        {
          std::cout << "Momenta:" << std::endl;
          for (int ipar = 0; ipar < npar; ipar++)
          {
            std::cout << std::setw(4) << ipar + 1
                      << setiosflags(std::ios::scientific) << std::setw(14)
                      << hstMomenta[ipagM*npar*np4*neppM + ipar*neppM*np4 + 0*neppM + ieppM] // AOSOA[ipagM][ipar][0][ieppM]
                      << setiosflags(std::ios::scientific) << std::setw(14)
                      << hstMomenta[ipagM*npar*np4*neppM + ipar*neppM*np4 + 1*neppM + ieppM] // AOSOA[ipagM][ipar][1][ieppM]
                      << setiosflags(std::ios::scientific) << std::setw(14)
                      << hstMomenta[ipagM*npar*np4*neppM + ipar*neppM*np4 + 2*neppM + ieppM] // AOSOA[ipagM][ipar][2][ieppM]
                      << setiosflags(std::ios::scientific) << std::setw(14)
                      << hstMomenta[ipagM*npar*np4*neppM + ipar*neppM*np4 + 3*neppM + ieppM] // AOSOA[ipagM][ipar][3][ieppM]
                      << std::endl;
          }
          std::cout << std::string(80, '-') << std::endl;
        }

        // Display matrix elements
        // FIXME: assume process.nprocesses == 1
        {
          if (verbose)
            std::cout << " Matrix element = "
              //   << setiosflags(ios::fixed) << setprecision(17)
                      << hstMEs[ievt] << " GeV^" << meGeVexponent << std::endl;
          if (perf)
            matrixelementALL[iiter*nevt + ievt] = hstMEs[ievt];
        }

        if (verbose)
          std::cout << std::string(80, '-') << std::endl;
      }
    }
    else if ( !debug )
    {
      std::cout << ".";
    }
  }

  // **************************************
  // *** END MAIN LOOP ON #ITERATIONS ***
  // **************************************

  // === STEP 9 FINALISE
  // --- 9a Dump after the loop
  const std::string dumpKey = "9a DumpAll ";
  timermap.start(dumpKey);

  if (!(verbose || debug || perf))
  {
    std::cout << std::endl;
  }

  if (perf)
  {
    double sumrtim = 0;
    double sqsrtim = 0;
    double minrtim = rambtimes[0];
    double maxrtim = rambtimes[0];
    for ( int iiter = 0; iiter < niter; ++iiter )
    {
      sumrtim += rambtimes[iiter];
      sqsrtim += rambtimes[iiter]*rambtimes[iiter];
      minrtim = std::min( minrtim, rambtimes[iiter] );
      maxrtim = std::max( maxrtim, rambtimes[iiter] );
    }
    //double meanrtim = sumrtim / niter; // unused
    //double stdrtim = std::sqrt( sqsrtim / niter - meanrtim * meanrtim ); // unused

    double sumwtim = 0;
    double sqswtim = 0;
    double minwtim = wavetimes[0];
    double maxwtim = wavetimes[0];
    for ( int iiter = 0; iiter < niter; ++iiter )
    {
      sumwtim += wavetimes[iiter];
      sqswtim += wavetimes[iiter]*wavetimes[iiter];
      minwtim = std::min( minwtim, wavetimes[iiter] );
      maxwtim = std::max( maxwtim, wavetimes[iiter] );
    }
    double meanwtim = sumwtim / niter;
    double stdwtim = std::sqrt( sqswtim / niter - meanwtim * meanwtim );

    int nnan = 0;
    double sumelem = 0;
    double sqselem = 0;
    double minelem = matrixelementALL[0];
    double maxelem = matrixelementALL[0];
    for ( int ievtALL = 0; ievtALL < nevtALL; ++ievtALL )
    {
      if ( isnan( matrixelementALL[ievtALL] ) )
      {
        if ( debug ) // only printed out with "-p -d" (matrixelementALL is not filled without -p)
          std::cout << "WARNING! ME[" << ievtALL << "} is nan" << std::endl;
        nnan++;
        continue;
      }
      sumelem += matrixelementALL[ievtALL];
      sqselem += matrixelementALL[ievtALL]*matrixelementALL[ievtALL];
      minelem = std::min( minelem, (double)matrixelementALL[ievtALL] );
      maxelem = std::max( maxelem, (double)matrixelementALL[ievtALL] );
    }
    double meanelem = sumelem / ( nevtALL - nnan );
    double stdelem = std::sqrt( sqselem / ( nevtALL - nnan) - meanelem * meanelem );

    std::cout << "***************************************" << std::endl
              << "NumIterations             = " << niter << std::endl
              << "NumThreadsPerBlock        = " << gputhreads << std::endl
              << "NumBlocksPerGrid          = " << gpublocks << std::endl
              << "---------------------------------------" << std::endl
#if defined MGONGPU_FPTYPE_DOUBLE
              << "FP precision              = DOUBLE (nan=" << nnan << ")" << std::endl
#elif defined MGONGPU_FPTYPE_FLOAT
              << "FP precision              = FLOAT (nan=" << nnan << ")" << std::endl
#endif
#ifdef __CUDACC__
#if defined MGONGPU_CXTYPE_CUCOMPLEX
              << "Complex type              = CUCOMPLEX" << std::endl
#elif defined MGONGPU_CXTYPE_THRUST
              << "Complex type              = THRUST::COMPLEX" << std::endl
#endif
#else
              << "Complex type              = STD::COMPLEX" << std::endl
#endif
              << "RanNumb memory layout     = AOSOA[" << neppR << "]" << std::endl
              << "Momenta memory layout     = AOSOA[" << neppM << "]" 
              << ( neppM == 1 ? " == AOS" : "" ) << std::endl
#ifdef __CUDACC__
#if defined MGONGPU_WFMEM_LOCAL
              << "Wavefunction GPU memory   = LOCAL" << std::endl
#elif defined MGONGPU_WFMEM_GLOBAL
              << "Wavefunction GPU memory   = GLOBAL" << std::endl
#elif defined MGONGPU_WFMEM_SHARED
              << "Wavefunction GPU memory   = SHARED (" << nbytesSharedSK/sizeof(char) << " bytes)" << std::endl
#endif
#endif
#ifdef __CUDACC__
#if defined MGONGPU_CURAND_ONDEVICE
              << "Curand generation         = DEVICE (CUDA code)" << std::endl
#elif defined MGONGPU_CURAND_ONHOST
              << "Curand generation         = HOST (CUDA code)" << std::endl
#endif
#else
              << "Curand generation         = HOST (C++ code)" << std::endl
#endif
              << "---------------------------------------" << std::endl
              << "NumberOfEntries           = " << niter << std::endl
              << std::scientific
              << "TotalTimeInWaveFuncs      = " << sumwtim << " sec" << std::endl
              << "MeanTimeInWaveFuncs       = " << meanwtim << " sec" << std::endl
              << "StdDevTimeInWaveFuncs     = " << stdwtim << " sec" << std::endl
              << "MinTimeInWaveFuncs        = " << minwtim << " sec" << std::endl
              << "MaxTimeInWaveFuncs        = " << maxwtim << " sec" << std::endl
              << "---------------------------------------" << std::endl
      //<< "ProcessID:                = " << getpid() << std::endl
      //<< "NProcesses                = " << process.nprocesses << std::endl
              << "TotalEventsComputed       = " << nevtALL << std::endl
              << "RamboEventsPerSec         = " << nevtALL/sumrtim << " sec^-1" << std::endl
              << "MatrixElemEventsPerSec    = " << nevtALL/sumwtim << " sec^-1" << std::endl;

    std::cout << "***************************************" << std::endl
              << "NumMatrixElements(notNan) = " << nevtALL - nnan << std::endl
              << std::scientific
              << "MeanMatrixElemValue       = " << meanelem << " GeV^" << meGeVexponent << std::endl
              << "StdErrMatrixElemValue     = " << stdelem/sqrt(nevtALL) << " GeV^" << meGeVexponent << std::endl
              << "StdDevMatrixElemValue     = " << stdelem << " GeV^" << meGeVexponent << std::endl
              << "MinMatrixElemValue        = " << minelem << " GeV^" << meGeVexponent << std::endl
              << "MaxMatrixElemValue        = " << maxelem << " GeV^" << meGeVexponent << std::endl;
  }

  // --- 9b. Destroy curand generator
  const std::string dgenKey = "9b GenDestr";
  timermap.start( dgenKey );
#ifdef __CUDACC__
  grambo2toNm0::destroyGenerator( rnGen );
#else
  rambo2toNm0::destroyGenerator( rnGen );
#endif

  // --- 9c Free memory structures
  const std::string freeKey = "9c MemFree ";
  timermap.start( freeKey );

#ifdef __CUDACC__
  checkCuda( cudaFreeHost( hstMEs ) );
  checkCuda( cudaFreeHost( hstWeights ) );
  checkCuda( cudaFreeHost( hstMomenta ) );
#if defined MGONGPU_CURAND_ONHOST
  checkCuda( cudaFreeHost( hstRnarray ) );
#endif
  checkCuda( cudaFree( devMEs ) );
#if defined MGONGPU_WFMEM_GLOBAL
  checkCuda( cudaFree( devAllWFs ) );
#endif
  checkCuda( cudaFree( devWeights ) );
  checkCuda( cudaFree( devMomenta ) );
  checkCuda( cudaFree( devRnarray ) );
#else
  delete[] hstMEs;
  delete[] hstWeights;
  delete[] hstMomenta;
  delete[] hstRnarray;
#endif

  delete[] rambtimes;
  delete[] wavetimes;
  delete[] matrixelementALL;

#ifdef __CUDACC__
  // --- 9d. Finalise cuda
  const std::string cdrsKey = "9d CudReset";
  timermap.start( cdrsKey );
  checkCuda( cudaDeviceReset() ); // this is needed by cuda-memcheck --leak-check full
#endif

  // *** STOP THE NEW TIMERS ***
  timermap.stop();
  if (perf)
  {
    std::cout << "***************************************" << std::endl;
    timermap.dump();
    std::cout << "***************************************" << std::endl;
  }

  //std::cout << "ALL OK" << std::endl;
}
