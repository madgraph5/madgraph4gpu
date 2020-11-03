#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
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
            << " [--verbose|-v] [--debug|-d] [--performance|-p] [--json|-j]"
            << " [#gpuBlocksPerGrid #gpuThreadsPerBlock] #iterations" << std::endl << std::endl;
  std::cout << "The number of events per iteration is #gpuBlocksPerGrid * #gpuThreadsPerBlock" << std::endl;
  std::cout << "(also in CPU/C++ code, where only the product of these two parameters counts)" << std::endl << std::endl;
  std::cout << "Summary stats are always computed: '-p' and '-j' only control their printout" << std::endl;
  std::cout << "The '-d' flag only controls if nan's emit warnings" << std::endl;
  return ret;
}

int main(int argc, char **argv)
{
  // READ COMMAND LINE ARGUMENTS
  bool verbose = false;
  bool debug = false;
  bool perf = false;
  bool json = false;
  int niter = 0;
  int gpublocks = 1;
  int gputhreads = 32;
  int jsondate = 0;
  int jsonrun = 0;
  int numvec[5] = {0,0,0,0,0};
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
    else if (strcmp(argv[argn], "--json") == 0 ||
             strcmp(argv[argn], "-j") == 0)
      json = true;
    else if (is_number(argv[argn]) && nnum<5)
      numvec[nnum++] = atoi(argv[argn]);
    else
      return usage(argv[0]);
  }

  if (nnum == 3 || nnum == 5) {
    gpublocks = numvec[0];
    gputhreads = numvec[1];
    niter = numvec[2];
    if (nnum == 5){
      jsondate = numvec[3];
      jsonrun = numvec[4];
    }
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
    std::cout << "ERROR! #threads/block should be a multiple of neppR=" << neppR << std::endl;
    return usage(argv[0]);
  }

  const int neppM = mgOnGpu::neppM; // ASA layout: constant at compile-time
  if ( gputhreads%neppM != 0 )
  {
    std::cout << "ERROR! #threads/block should be a multiple of neppM=" << neppM << std::endl;
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
  const fptype energy = 1500; // historical default, Ecms = 1500 GeV = 1.5 TeV (above the Z peak)
  //const fptype energy = 91.2; // Ecms = 91.2 GeV (Z peak)
  //const fptype energy = 0.100; // Ecms = 100 MeV (well below the Z peak, pure em scattering)
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

#ifdef __CUDACC__
  using mgOnGpu::ncomb;
  const int nbytesIsGoodHel = ncomb * sizeof(bool);
  bool* hstIsGoodHel = 0;
  checkCuda( cudaMallocHost( &hstIsGoodHel, nbytesIsGoodHel ) );
  bool* devIsGoodHel = 0;
  checkCuda( cudaMalloc( &devIsGoodHel, nbytesIsGoodHel ) );
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
  fptype* weightALL = new fptype[nevtALL]();

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
    // 3a. (Only on the first iteration) Get good helicities
    // 3b. Evaluate MEs on the device
    // 3c. Copy MEs back from device to host

    // --- 3a. SGoodHel
#ifdef __CUDACC__
    if ( iiter == 0 )
    {
      const std::string ghelKey = "3a SGoodHel";
      timermap.start( ghelKey );
      // ... 3a1. Compute good helicity mask on the device
      gProc::sigmaKin_getGoodHel<<<gpublocks, gputhreads>>>(devMomenta, devIsGoodHel);
      checkCuda( cudaPeekAtLastError() );
      // ... 3a2. Copy back good helicity mask to the host
      checkCuda( cudaMemcpy( hstIsGoodHel, devIsGoodHel, nbytesIsGoodHel, cudaMemcpyDeviceToHost ) );
      // ... 3a3. Copy back good helicity list to constant memory on the device
      gProc::sigmaKin_setGoodHel(hstIsGoodHel);
    }
#endif

    // *** START THE OLD TIMER FOR WAVEFUNCTIONS ***
    double wavetime = 0;

    // --- 3b. SigmaKin
    const std::string skinKey = "3b SigmaKin";
    timermap.start( skinKey );
#ifdef __CUDACC__
#ifndef MGONGPU_NSIGHT_DEBUG
    gProc::sigmaKin<<<gpublocks, gputhreads>>>(devMomenta, devMEs);
#else
    gProc::sigmaKin<<<gpublocks, gputhreads, ntpbMAX*sizeof(float)>>>(devMomenta, devMEs);
#endif
    checkCuda( cudaPeekAtLastError() );
#else
    Proc::sigmaKin(hstMomenta, hstMEs, nevt);
#endif

#ifdef __CUDACC__
    // --- 3c. CopyDToH MEs
    const std::string cmesKey = "3c CpDTHmes";
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

    for (int ievt = 0; ievt < nevt; ++ievt) // Loop over all events in this iteration
    {
      if (verbose)
      {
        // Display momenta
        const int ipagM = ievt/neppM; // #eventpage in this iteration
        const int ieppM = ievt%neppM; // #event in the current eventpage in this iteration
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
        // Display matrix elements
        std::cout << " Matrix element = "
          //   << setiosflags(ios::fixed) << setprecision(17)
                  << hstMEs[ievt] << " GeV^" << meGeVexponent << std::endl; // FIXME: assume process.nprocesses == 1
        std::cout << std::string(80, '-') << std::endl;
      }
      // Fill the arrays with ALL MEs and weights
      matrixelementALL[iiter*nevt + ievt] = hstMEs[ievt]; // FIXME: assume process.nprocesses == 1
      weightALL[iiter*nevt + ievt] = hstWeights[ievt];
    }

    if (!(verbose || debug || perf))
    {
      std::cout << ".";
    }

  }

  // **************************************
  // *** END MAIN LOOP ON #ITERATIONS ***
  // **************************************

  // === STEP 8 ANALYSIS
  // --- 8a Analysis: compute stats after the loop
  const std::string statKey = "8a CompStat";
  timermap.start(statKey);

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
  double sumweig = 0;
  double sqsweig = 0;
  double minweig = weightALL[0];
  double maxweig = weightALL[0];
  for ( int ievtALL = 0; ievtALL < nevtALL; ++ievtALL )
  {
    if ( std::isnan( matrixelementALL[ievtALL] ) )
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
    sumweig += weightALL[ievtALL];
    sqsweig += weightALL[ievtALL]*weightALL[ievtALL];
    minweig = std::min( minweig, (double)weightALL[ievtALL] );
    maxweig = std::max( maxweig, (double)weightALL[ievtALL] );
  }
  double meanelem = sumelem / ( nevtALL - nnan );
  double stdelem = std::sqrt( sqselem / ( nevtALL - nnan) - meanelem * meanelem );
  double meanweig = sumweig / ( nevtALL - nnan );
  double stdweig = std::sqrt( sqsweig / ( nevtALL - nnan) - meanweig * meanweig );

  // === STEP 9 FINALISE
  // --- 9a. Destroy curand generator
  const std::string dgenKey = "9a GenDestr";
  timermap.start( dgenKey );
#ifdef __CUDACC__
  grambo2toNm0::destroyGenerator( rnGen );
#else
  rambo2toNm0::destroyGenerator( rnGen );
#endif

  // --- 9b Free memory structures
  const std::string freeKey = "9b MemFree ";
  timermap.start( freeKey );

#ifdef __CUDACC__
  checkCuda( cudaFreeHost( hstMEs ) );
  checkCuda( cudaFreeHost( hstIsGoodHel ) );
  checkCuda( cudaFreeHost( hstWeights ) );
  checkCuda( cudaFreeHost( hstMomenta ) );
#if defined MGONGPU_CURAND_ONHOST
  checkCuda( cudaFreeHost( hstRnarray ) );
#endif
  checkCuda( cudaFree( devMEs ) );
  checkCuda( cudaFree( devIsGoodHel ) );
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
  delete[] weightALL;

#ifdef __CUDACC__
  // --- 9c. Finalise cuda
  const std::string cdrsKey = "9c CudReset";
  timermap.start( cdrsKey );
  checkCuda( cudaDeviceReset() ); // this is needed by cuda-memcheck --leak-check full
#endif

  // --- 9d Dump to screen
  const std::string dumpKey = "9d DumpScrn";
  timermap.start(dumpKey);

  if (!(verbose || debug || perf))
  {
    std::cout << std::endl;
  }

  if (perf)
  {
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
              << "RanNumb memory layout     = AOSOA[" << neppR << "]"
              << ( neppR == 1 ? " == AOS" : "" ) << std::endl
              << "Momenta memory layout     = AOSOA[" << neppM << "]"
              << ( neppM == 1 ? " == AOS" : "" ) << std::endl
#ifdef __CUDACC__
              << "Wavefunction GPU memory   = LOCAL" << std::endl
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
              << "MaxMatrixElemValue        = " << maxelem << " GeV^" << meGeVexponent << std::endl
              << "MeanWeight                = " << meanweig << std::endl
              << "StdErrWeight              = " << stdweig/sqrt(nevtALL) << std::endl
              << "StdDevWeight              = " << stdweig << std::endl
              << "MinWeight                 = " << minweig << std::endl
              << "MaxWeight                 = " << maxweig << std::endl;
  }

  // --- 9e Dump to json
  const std::string jsonKey = "9e DumpJson";
  timermap.start(jsonKey);

  if(json)
  {
    std::string jsonFileName = std::to_string(jsondate) + "-perf-test-run" + std::to_string(jsonrun) + ".json";
    jsonFileName = "./perf/data/" + jsonFileName;

    //Checks if file exists
    std::ifstream fileCheck;
    bool fileExists = false;
    fileCheck.open(jsonFileName);
    if(fileCheck){
      fileExists = true;
      fileCheck.close();
    }

    std::ofstream jsonFile;
    jsonFile.open(jsonFileName, std::ios_base::app);
    if(!fileExists){
      jsonFile << "[" << std::endl;
    }
    else{
      //deleting the last bracket and outputting a ", "
      std::string temp = "truncate -s-1 " + jsonFileName;
      const char *command = temp.c_str();
      system(command);
      jsonFile << ", " << std::endl;
    }

    jsonFile << "{" << std::endl
             << "\"NumIterations\": " << niter << ", " << std::endl
             << "\"NumThreadsPerBlock\": " << gputhreads << ", " << std::endl
             << "\"NumBlocksPerGrid\": " << gpublocks << ", " << std::endl
#if defined MGONGPU_FPTYPE_DOUBLE
             << "\"FP precision\": "
             << "\"DOUBLE (nan=" << nnan << ")\"," << std::endl
#elif defined MGONGPU_FPTYPE_FLOAT
             << "\"FP precision\": " << "FLOAT (nan=" << nnan << ")," << std::endl
#endif
             << "\"Complex type\": "
#ifdef __CUDACC__
#if defined MGONGPU_CXTYPE_CUCOMPLEX
             << "\"CUCOMPLEX\"," << std::endl
#elif defined MGONGPU_CXTYPE_THRUST
             << "\"THRUST::COMPLEX\"," << std::endl
#endif
#else
             << "\"STD::COMPLEX\"," << std::endl
#endif
             << "\"RanNumb memory layout\": " << "\"AOSOA[" << neppR << "]\""
             << ( neppR == 1 ? " == AOS" : "" ) << ", " << std::endl
             << "\"Momenta memory layout\": " << "\"AOSOA[" << neppM << "]\""
             << ( neppM == 1 ? " == AOS" : "" ) << ", " << std::endl
#ifdef __CUDACC__
             << "\"Wavefunction GPU memory\": " << "\"LOCAL\"," << std::endl
#endif
             << "\"Curand generation\": "
#ifdef __CUDACC__
#if defined MGONGPU_CURAND_ONDEVICE
             << "\"DEVICE (CUDA code)\"," << std::endl
#elif defined MGONGPU_CURAND_ONHOST
             << "\"HOST (CUDA code)\"," << std::endl
#endif
#else
             << "\"HOST (C++ code)\"," << std::endl
#endif
             << "\"NumberOfEntries\": " << niter << "," << std::endl
      //<< std::scientific // Not sure about this
             << "\"TotalTimeInWaveFuncs\": "
             << "\"" << std::to_string(sumwtim) << " sec\"," << std::endl
             << "\"MeanTimeInWaveFuncs\": "
             << "\"" << std::to_string(meanwtim) << " sec\"," << std::endl
             << "\"StdDevTimeInWaveFuncs\": "
             << "\"" << std::to_string(stdwtim) << " sec\"," << std::endl
             << "\"MinTimeInWaveFuncs\": "
             << "\"" << std::to_string(minwtim) << " sec\"," << std::endl
             << "\"MaxTimeInWaveFuncs\": "
             << "\"" << std::to_string(maxwtim) << " sec\"," << std::endl
      //<< "ProcessID:                = " << getpid() << std::endl
      //<< "NProcesses                = " << process.nprocesses << std::endl
             << "\"TotalEventsComputed\": " << nevtALL << "," << std::endl
             << "\"RamboEventsPerSec\": "
             << "\"" << std::to_string(nevtALL/sumrtim) << " sec^-1\"," << std::endl
             << "\"MatrixElemEventsPerSec\": "
             << "\"" << std::to_string(nevtALL/sumwtim) << " sec^-1\"," << std::endl
             << "\"NumMatrixElements(notNan)\": " << nevtALL - nnan << "," << std::endl
             << std::scientific
             << "\"MeanMatrixElemValue\": "
             << "\"" << std::to_string(meanelem) << " GeV^"
             << std::to_string(meGeVexponent) << "\"," << std::endl
             << "\"StdErrMatrixElemValue\": "
             << "\"" << std::to_string(stdelem/sqrt(nevtALL)) << " GeV^"
             << std::to_string(meGeVexponent) << "\"," << std::endl
             << "\"StdDevMatrixElemValue\": "
             << "\"" << std::to_string(stdelem)
             << " GeV^" << std::to_string(meGeVexponent) << "\"," << std::endl
             << "\"MinMatrixElemValue\": "
             << "\"" << std::to_string(minelem) << " GeV^"
             << std::to_string(meGeVexponent) << "\"," << std::endl
             << "\"MaxMatrixElemValue\": "
             << "\"" << std::to_string(maxelem) << " GeV^"
             << std::to_string(meGeVexponent) <<  "\"," << std::endl;

    timermap.dump(jsonFile, true); // NB For the active json timer this dumps a partial total

    jsonFile << "}" << std::endl << "]";
    jsonFile.close();
  }

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
