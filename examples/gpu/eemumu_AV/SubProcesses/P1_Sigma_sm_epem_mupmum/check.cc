#include <algorithm> // perf stats
#include <cstring>
#include <iomanip>
#include <iostream>
#include <numeric> // perf stats
#include <unistd.h>

#include "mgOnGpuConfig.h"

#ifdef __CUDACC__
#include "grambo2toNm0.cu"
#else
#include "rambo2toNm0.h"
#endif

#include "CPPProcess.h"
#include "timermap.h"

bool is_number(const char *s) {
  const char *t = s;
  while (*t != '\0' && isdigit(*t))
    ++t;
  return strlen(s) == t - s;
}

int usage(char* argv0, int ret = 1) {
  std::cout << "Usage: " << argv0
            << " [--verbose|-v] [--debug|-d] [--performance|-p]"
            << " [#gpuBlocksPerGrid #gpuThreadsPerBlock] #iterations" << std::endl;
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

  const int ndim = gpublocks * gputhreads; // number of events (threads) in one iteration
#if defined MGONGPU_LAYOUT_ASA
  using mgOnGpu::nepp;
  if ( gputhreads%nepp != 0 )
  {
    std::cout << "ERROR! #threads/block should be a multiple of " << nepp << std::endl;
    return usage(argv[0]);
  }
#endif

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
  gpuErrchk3( cudaFree( 0 ) ); // SLOW!
  //std::cout << "Calling cudaFree... done" << std::endl;
#endif

  // --- 0a. Initialise physics process
  const std::string procKey = "0a ProcInit";
  timermap.start( procKey );

  // Create a process object
#ifdef __CUDACC__
  gProc::CPPProcess process(niter, gpublocks, gputhreads, verbose, debug);
#else
  Proc::CPPProcess process(niter, gpublocks, gputhreads, verbose, debug);
#endif

  // Read param_card and set parameters
  process.initProc("../../Cards/param_card.dat");

  const double energy = 1500;
  const int meGeVexponent = -(2 * process.nexternal - 8);

  // --- 0b. Allocate memory structures
  const std::string alloKey = "0b MemAlloc";
  timermap.start( alloKey );

  // Memory structures for random numbers, momenta, matrix elements and weights on host and device
  const int npar = process.nexternal; // for this process (eemumu): npar=4 (e+, e-, mu+, mu-)
  const int nparf = npar - process.ninitial; // for this process (eemumu): nparf=2 (mu+, mu-)
  const int np4 = 4; // dimension of 4-momenta (E,px,py,pz): copy all of them from rambo

  const int nbytesRnarray = np4*nparf*ndim * sizeof(double); // (NB: ndim=npag*nepp for ASA layouts)
#ifdef __CUDACC__
  double* devRnarray = 0; // AOSOA[npag][nparf][np4][nepp] (NB: ndim=npag*nepp)
  gpuErrchk3( cudaMalloc( &devRnarray, nbytesRnarray ) );
#if defined MGONGPU_CURAND_ONHOST
  double* hstRnarray = 0; // AOSOA[npag][nparf][np4][nepp] (NB: ndim=npag*nepp)
  gpuErrchk3( cudaMallocHost( &hstRnarray, nbytesRnarray ) );
#endif
#else
  double* hstRnarray = 0; // AOSOA[npag][nparf][np4][nepp] (NB: ndim=npag*nepp)
  gpuErrchk3( cudaMallocHost( &hstRnarray, nbytesRnarray ) );
#endif

  const int nbytesMomenta = np4*npar*ndim * sizeof(double); // (NB: ndim=npag*nepp for ASA layouts)
#if defined MGONGPU_LAYOUT_ASA
  double* hstMomenta = 0; // AOSOA[npag][npar][np4][nepp] (previously was: lp)
#elif defined MGONGPU_LAYOUT_SOA
  double* hstMomenta = 0; // SOA[npar][np4][ndim] (previously was: lp)
#elif defined MGONGPU_LAYOUT_AOS
  double* hstMomenta = 0; // AOS[ndim][npar][np4] (previously was: lp)
#endif
  gpuErrchk3( cudaMallocHost( &hstMomenta, nbytesMomenta ) );
  double* devMomenta = 0; // (previously was: allMomenta)
  gpuErrchk3( cudaMalloc( &devMomenta, nbytesMomenta ) );

  const int nbytesWeights = ndim * sizeof(double); //  (NB: ndim=npag*nepp for ASA layouts)
  double* hstWeights = 0; // (previously was: meHostPtr)
  gpuErrchk3( cudaMallocHost( &hstWeights, nbytesWeights ) );
  double* devWeights = 0; // (previously was: meDevPtr)
  gpuErrchk3( cudaMalloc( &devWeights, nbytesWeights ) );

  const int nbytesMEs = ndim * sizeof(double); //  (NB: ndim=npag*nepp for ASA layouts)
  double* hstMEs = 0; // (previously was: meHostPtr)
  gpuErrchk3( cudaMallocHost( &hstMEs, nbytesMEs ) );
  double* devMEs = 0; // (previously was: meDevPtr)
  gpuErrchk3( cudaMalloc( &devMEs, nbytesMEs ) );

  float* wavetimes = new float[niter]();
  double* matrixelementvector = new double[niter * ndim * process.nprocesses]();

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

    // --- 1b. Generate all relevant numbers to build ndim events (i.e. ndim phase space points) on the host
    const std::string rngnKey = "1b GenRnGen";
    timermap.start( rngnKey );
#ifdef __CUDACC__
#if defined MGONGPU_CURAND_ONDEVICE
    grambo2toNm0::generateRnArray( rnGen, devRnarray, ndim );
#elif defined MGONGPU_CURAND_ONHOST
    grambo2toNm0::generateRnArray( rnGen, hstRnarray, ndim );
#endif
#else
    rambo2toNm0::generateRnArray( rnGen, hstRnarray, ndim );
#endif
    //std::cout << "Got random numbers" << std::endl;

#ifdef __CUDACC__
#if defined MGONGPU_CURAND_ONHOST
    // --- 1c. Copy rnarray from host to device
    const std::string htodKey = "1c CpHTDrnd";
    timermap.start( htodKey );
    gpuErrchk3( cudaMemcpy( devRnarray, hstRnarray, nbytesRnarray, cudaMemcpyHostToDevice ) );
#endif
#endif

    // === STEP 2 OF 3
    // Fill in particle momenta for each of ndim events on the device

    // --- 2a. Fill in momenta of initial state particles on the device
    const std::string riniKey = "2a RamboIni";
    timermap.start( riniKey );
#ifdef __CUDACC__
    grambo2toNm0::getMomentaInitial<<<gpublocks, gputhreads>>>( energy, devMomenta, ndim );
#else
    rambo2toNm0::getMomentaInitial( energy, hstMomenta, ndim );
#endif
    //std::cout << "Got initial momenta" << std::endl;

    // --- 2b. Fill in momenta of final state particles using the RAMBO algorithm on the device
    // (i.e. map random numbers to final-state particle momenta for each of ndim events)
    const std::string rfinKey = "2b RamboFin";
    timermap.start( rfinKey );
#ifdef __CUDACC__
    grambo2toNm0::getMomentaFinal<<<gpublocks, gputhreads>>>( energy, devRnarray, devMomenta, devWeights, ndim );
#else
    rambo2toNm0::getMomentaFinal( energy, hstRnarray, hstMomenta, hstWeights, ndim );
#endif
    //std::cout << "Got final momenta" << std::endl;

    // --- 2c. CopyDToH Weights
    const std::string cwgtKey = "2c CpDTHwgt";
    timermap.start( cwgtKey );
    gpuErrchk3( cudaMemcpy( hstWeights, devWeights, nbytesWeights, cudaMemcpyDeviceToHost ) );

    // --- 2d. CopyDToH Momenta
    const std::string cmomKey = "2d CpDTHmom";
    timermap.start( cmomKey );
    gpuErrchk3( cudaMemcpy( hstMomenta, devMomenta, nbytesMomenta, cudaMemcpyDeviceToHost ) );

    // === STEP 3 OF 3
    // Evaluate matrix elements for all ndim events
    // 3a. Evaluate MEs on the device
    // 3b. Copy MEs back from device to host

    // *** START THE OLD TIMER ***
    float gputime = 0;

    // --- 3a. SigmaKin
    const std::string skinKey = "3a SigmaKin";
    timermap.start( skinKey );
#ifdef __CUDACC__
    gProc::sigmaKin<<<gpublocks, gputhreads>>>(devMomenta,  devMEs);//, debug, verbose);
#else
    Proc::sigmaKin<<<gpublocks, gputhreads>>>(devMomenta,  devMEs);//, debug, verbose);
#endif
    gpuErrchk3( cudaPeekAtLastError() );

    // --- 3b. CopyDToH MEs
    const std::string cmesKey = "3b CpDTHmes";
    gputime += timermap.start( cmesKey );
    gpuErrchk3( cudaMemcpy( hstMEs, devMEs, nbytesMEs, cudaMemcpyDeviceToHost ) );

    // *** STOP THE OLD TIMER ***
    gputime += timermap.stop();

    // === STEP 4 FINALISE LOOP
    // --- 4a Dump within the loop
    const std::string loopKey = "4a DumpLoop";
    timermap.start(loopKey);
    wavetimes[iiter] = gputime;

    if (verbose)
    {
      std::cout << "***********************************" << std::endl
                << "Iteration #" << iiter+1 << " of " << niter << std::endl;
      if (perf) std::cout << "Wave function time: " << gputime << std::endl;
    }

    if (verbose || perf)
    {
      for (int idim = 0; idim < ndim; ++idim) // Loop over all events in this iteration
      {
#if defined MGONGPU_LAYOUT_ASA
        const int ipag = idim/nepp; // #eventpage in this iteration
        const int iepp = idim%nepp; // #event in the current eventpage in this iteration
#endif
        if (verbose)
        {
          std::cout << "Momenta:" << std::endl;
          for (int ipar = 0; ipar < npar; ipar++)
          {
#if defined MGONGPU_LAYOUT_ASA
            std::cout << std::setw(4) << ipar + 1
                      << setiosflags(std::ios::scientific) << std::setw(14)
                      << hstMomenta[ipag*npar*np4*nepp + ipar*nepp*np4 + 0*nepp + iepp] // AOSOA[ipag][ipar][0][iepp]
                      << setiosflags(std::ios::scientific) << std::setw(14)
                      << hstMomenta[ipag*npar*np4*nepp + ipar*nepp*np4 + 1*nepp + iepp] // AOSOA[ipag][ipar][1][iepp]
                      << setiosflags(std::ios::scientific) << std::setw(14)
                      << hstMomenta[ipag*npar*np4*nepp + ipar*nepp*np4 + 2*nepp + iepp] // AOSOA[ipag][ipar][2][iepp]
                      << setiosflags(std::ios::scientific) << std::setw(14)
                      << hstMomenta[ipag*npar*np4*nepp + ipar*nepp*np4 + 3*nepp + iepp] // AOSOA[ipag][ipar][3][iepp]
                      << std::endl;
#elif defined MGONGPU_LAYOUT_SOA
            std::cout << std::setw(4) << ipar + 1
                      << setiosflags(std::ios::scientific) << std::setw(14)
                      << hstMomenta[ipar*ndim*np4 + 0*ndim + idim] // SOA[ipar][0][idim]
                      << setiosflags(std::ios::scientific) << std::setw(14)
                      << hstMomenta[ipar*ndim*np4 + 1*ndim + idim] // SOA[ipar][1][idim]
                      << setiosflags(std::ios::scientific) << std::setw(14)
                      << hstMomenta[ipar*ndim*np4 + 2*ndim + idim] // SOA[ipar][2][idim]
                      << setiosflags(std::ios::scientific) << std::setw(14)
                      << hstMomenta[ipar*ndim*np4 + 3*ndim + idim] // SOA[ipar][3][idim]
                      << std::endl;
#elif defined MGONGPU_LAYOUT_AOS
            std::cout << std::setw(4) << ipar + 1
                      << setiosflags(std::ios::scientific) << std::setw(14)
                      << hstMomenta[idim*npar*np4 + ipar*np4 + 0] // AOS[idim][ipar][0]
                      << setiosflags(std::ios::scientific) << std::setw(14)
                      << hstMomenta[idim*npar*np4 + ipar*np4 + 1] // AOS[idim][ipar][1]
                      << setiosflags(std::ios::scientific) << std::setw(14)
                      << hstMomenta[idim*npar*np4 + ipar*np4 + 2] // AOS[idim][ipar][2]
                      << setiosflags(std::ios::scientific) << std::setw(14)
                      << hstMomenta[idim*npar*np4 + ipar*np4 + 3] // AOS[idim][ipar][3]
                      << std::endl;
#endif
          }
          std::cout << std::string(80, '-') << std::endl;
        }

        // Display matrix elements
        // FIXME: assume process.nprocesses == 1
        {
          if (verbose)
            std::cout << " Matrix element = "
              //   << setiosflags(ios::fixed) << setprecision(17)
                      << hstMEs[idim] << " GeV^" << meGeVexponent << std::endl;
          if (perf)
            matrixelementvector[iiter*ndim + idim] = hstMEs[idim];
        }

        if (verbose)
          std::cout << std::string(80, '-') << std::endl;
      }
    }
    else if (!debug)
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

    float sum = 0;
    float sq_sum = 0;
    float mintime = wavetimes[0];
    float maxtime = wavetimes[0];
    for (int iiter = 0; iiter < niter; ++iiter)
    {
      sum += wavetimes[iiter];
      sq_sum += wavetimes[iiter]*wavetimes[iiter];
      mintime = std::min( mintime, wavetimes[iiter] );
      maxtime = std::max( mintime, wavetimes[iiter] );
    }
    float mean = sum / niter;
    float stdev = std::sqrt( sq_sum / niter - mean * mean );

    int num_mes = niter*ndim;
    float sumelem = 0;
    float sqselem = 0;
    float minelem = matrixelementvector[0];
    float maxelem = matrixelementvector[0];
    for (int imes = 0; imes < num_mes; ++imes)
    {
      sumelem += matrixelementvector[imes];
      sqselem += matrixelementvector[imes]*matrixelementvector[imes];
      minelem = std::min( mintime, (float)matrixelementvector[imes] );
      maxelem = std::max( mintime, (float)matrixelementvector[imes] );
    }
    float meanelem = sumelem / num_mes;
    float stdelem = std::sqrt( sqselem / num_mes - meanelem * meanelem );

    std::cout << "***********************************" << std::endl
              << "NumIterations         = " << niter << std::endl
              << "NumThreadsPerBlock    = " << gputhreads << std::endl
              << "NumBlocksPerGrid      = " << gpublocks << std::endl
              << "-----------------------------------" << std::endl
#if defined MGONGPU_LAYOUT_ASA
              << "Momenta memory layout = AOSOA" << std::endl
#elif defined MGONGPU_LAYOUT_SOA
              << "Momenta memory layout = SOA" << std::endl
#elif defined MGONGPU_LAYOUT_AOS
              << "Momenta memory layout = AOS" << std::endl
#endif
#ifdef __CUDACC__
#if defined MGONGPU_CURAND_ONDEVICE
              << "Curand generation     = DEVICE" << std::endl
#elif defined MGONGPU_CURAND_ONHOST
              << "Curand generation     = HOST" << std::endl
#endif
#else
              << "Curand generation     = HOST (C++ code)" << std::endl
#endif
              << "-----------------------------------" << std::endl
              << "NumberOfEntries       = " << niter << std::endl
              << std::scientific
              << "TotalTimeInWaveFuncs  = " << sum << " sec" << std::endl
              << "MeanTimeInWaveFuncs   = " << mean << " sec" << std::endl
              << "StdDevTimeInWaveFuncs = " << stdev << " sec" << std::endl
              << "MinTimeInWaveFuncs    = " << mintime << " sec" << std::endl
              << "MaxTimeInWaveFuncs    = " << maxtime << " sec" << std::endl
              << "-----------------------------------" << std::endl
              << "ProcessID:            = " << getpid() << std::endl
              << "NProcesses            = " << process.nprocesses << std::endl
              << "NumMatrixElements     = " << num_mes << std::endl
              << "MatrixElementsPerSec  = " << num_mes/sum << " sec^-1" << std::endl;

    std::cout << "***********************************" << std::endl
              << "NumMatrixElements     = " << num_mes << std::endl
              << std::scientific
              << "MeanMatrixElemValue   = " << meanelem << " GeV^" << meGeVexponent << std::endl
              << "StdErrMatrixElemValue = " << stdelem/sqrt(num_mes) << " GeV^" << meGeVexponent << std::endl
              << "StdDevMatrixElemValue = " << stdelem << " GeV^" << meGeVexponent << std::endl
              << "MinMatrixElemValue    = " << minelem << " GeV^" << meGeVexponent << std::endl
              << "MaxMatrixElemValue    = " << maxelem << " GeV^" << meGeVexponent << std::endl;
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

  gpuErrchk3( cudaFreeHost( hstMEs ) );
  gpuErrchk3( cudaFreeHost( hstWeights ) );
  gpuErrchk3( cudaFreeHost( hstMomenta ) );

  gpuErrchk3( cudaFree( devMEs ) );
  gpuErrchk3( cudaFree( devWeights ) );
  gpuErrchk3( cudaFree( devMomenta ) );

#ifdef __CUDACC__
  gpuErrchk3( cudaFree( devRnarray ) );
#if defined MGONGPU_CURAND_ONHOST
  gpuErrchk3( cudaFreeHost( hstRnarray ) );
#endif
#else
  gpuErrchk3( cudaFreeHost( hstRnarray ) );
#endif

  delete[] wavetimes;
  delete[] matrixelementvector;

  // --- 9d. Finalise cuda
  const std::string cdrsKey = "9d CudReset";
  timermap.start( cdrsKey );
  gpuErrchk3( cudaDeviceReset() ); // this is needed by cuda-memcheck --leak-check full

  // *** STOP THE NEW TIMERS ***
  timermap.stop();
  if (perf)
  {
    std::cout << "***********************************" << std::endl;
    timermap.dump();
    std::cout << "***********************************" << std::endl;
  }

  //std::cout << "ALL OK" << std::endl;
}
