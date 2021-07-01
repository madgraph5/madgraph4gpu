#include <algorithm> // perf stats
#include <cstring>
#include <iomanip>
#include <iostream>
#include <numeric> // perf stats
#include <unistd.h>
#include <vector>

#include "CPPProcess.h"
//#include "HelAmps_sm.h"

#include "rambo.h"
#include "random_generator.h"
#include "timermap.h"

#include "Kokkos_Core.hpp"

#define TIMERTYPE std::chrono::high_resolution_clock

bool is_number(const char *s) {
  const char *t = s;
  while (*t != '\0' && isdigit(*t))
    ++t;
  return strlen(s) == t - s;
}

int usage(char* argv0, int ret = 1) {
  std::cout << "Usage: " << argv0 
            << " [--verbose|-v] [--debug|-d] [--performance|-p]"
            << " #events_per_iter #iterations" << std::endl;
  return ret;
}

int main(int argc, char **argv) {
  bool verbose = false, debug = false, perf = false;
  int num_iter = 0, events_per_iter = 0;
  std::vector<int> numvec;
  Timer<TIMERTYPE> timer;
  std::vector<float> wavetimes;


  for (int argn = 1; argn < argc; ++argn) {
    if (strcmp(argv[argn], "--verbose") == 0 || strcmp(argv[argn], "-v") == 0)
      verbose = true;
    else if (strcmp(argv[argn], "--debug") == 0 ||
             strcmp(argv[argn], "-d") == 0)
      debug = true;
    else if (strcmp(argv[argn], "--performance") == 0 ||
             strcmp(argv[argn], "-p") == 0)
      perf = true;
    else if (is_number(argv[argn]))
      numvec.push_back(atoi(argv[argn]));
    else
      return usage(argv[0]);
  }
  int veclen = numvec.size();
  if (veclen == 2) {
    events_per_iter = numvec[0];
    num_iter = numvec[1];
  } else {
    return usage(argv[0]);
  }

  if (num_iter == 0)
    return usage(argv[0]);

  Kokkos::initialize(argc, argv);
  if (verbose)
    std::cout << "# iterations: " << num_iter << std::endl;
  { // start Kokkos View space
    Kokkos::Timer ktimer;
    // *** START THE NEW TIMERS ***
    MG5::TimerMap timermap;

    // --- 0a. Initialise physics process
    const std::string procKey = "0a ProcInit";
    timermap.start( procKey );

    // Create a process object
    CPPProcess<Kokkos::DefaultExecutionSpace> process(num_iter, events_per_iter);

    // Read param_card and set parameters
    process.initProc("../../Cards/param_card.dat");

    double energy = 1500;

    int meGeVexponent = -(2 * process.nexternal - 8);
    
    double nevtALL = events_per_iter * num_iter;

    // --- 0b. Allocate memory structures
    const std::string alloKey = "0b MemAlloc";
    timermap.start( alloKey );

    Kokkos::View<double*,Kokkos::DefaultExecutionSpace> meDevPtr("meDevPtr",events_per_iter);
    auto meHostPtr = Kokkos::create_mirror_view(meDevPtr);

    Kokkos::View<double,Kokkos::DefaultExecutionSpace> d_wgt("d_wgt");
    auto h_wgt = Kokkos::create_mirror_view(d_wgt);

    const int nprocesses = 1; // TODO: hardcoded value
    Kokkos::View<double**,Kokkos::DefaultExecutionSpace> matrix_element("matrix_element",events_per_iter,nprocesses);
    Kokkos::View<double**,Kokkos::DefaultExecutionSpace> random_numbers(Kokkos::ViewAllocateWithoutInitializing("rns"),events_per_iter,4*(process.nexternal - process.ninitial));
    Kokkos::View<double***,Kokkos::DefaultExecutionSpace> p(Kokkos::ViewAllocateWithoutInitializing("p"),events_per_iter,process.nexternal,4);
    get_initial_momenta(p,process.nexternal,energy,process.cmME,events_per_iter);
    auto h_p = Kokkos::create_mirror_view(p);

    // used in calculate_wavefunctions
    Kokkos::View<Kokkos::complex<double>**,Kokkos::DefaultExecutionSpace> amp("amp",events_per_iter,2);
    Kokkos::View<Kokkos::complex<double>***,Kokkos::DefaultExecutionSpace> w("w",events_per_iter,5,6);


    std::vector<double> matrixelementvector;
    std::vector<double> genrtimes;
    std::vector<double> rambtimes;
    std::vector<double> wavetimes;
    std::vector<double> matrixelementALL;
    std::vector<double> weightALL;
    Kokkos::Timer sigmaKin_timer;
    for (int x = 0; x < num_iter; ++x) {
      // printf("iter %d of %d\n",x,num_iter);

      // get random numbers
      double genrtime = 0;
      const std::string sgenKey = "1 GenRnGen";
      timermap.start( sgenKey );
      fill_random_numbers_2d(random_numbers,events_per_iter,4*(process.nexternal - process.ninitial));
      genrtime += timermap.stop();

      double rambtime = 0;
      const std::string riniKey = "2a RamboIni";
      timermap.start( riniKey );
      // Get input particle momenta
      Kokkos::deep_copy(p,h_p);


      const std::string rfinKey = "2b RamboFin";
      rambtime += timermap.start( rfinKey );
      // Get phase space point
      get_final_momenta(process.ninitial, process.nexternal, energy, process.cmME, p, random_numbers, d_wgt, events_per_iter);
      

      const std::string cwgtKey = "2c CpDTHwgt";
      rambtime += timermap.start( cwgtKey );
      Kokkos::deep_copy(h_wgt,d_wgt);

      const std::string cmomKey = "2d CpDTHmom";
      rambtime += timermap.start( cmomKey );
      auto h_p = Kokkos::create_mirror_view(p);
      Kokkos::deep_copy(h_p,p);

      rambtime += timermap.stop();

      if (perf) {
        timer.Start();
      }
      // Evaluate matrix element
      // later process.sigmaKin(ncomb, goodhel, ntry, sum_hel, ngood, igood,
      // jhel);
      double wavetime = 0;
      const std::string skinKey = "3a SigmaKin";
      timermap.start( skinKey );
      // sigmaKin_timer.reset();
      sigmaKin(p, meDevPtr, process.cHel, process.cIPD, process.cIPC, amp, w, matrix_element, events_per_iter);//, debug, verbose);
      // printf("sigmaKin_timer: %f\n",sigmaKin_timer.seconds());
      wavetime += timermap.stop();
      

      // auto hp = Kokkos::create_mirror_view(p);
      // Kokkos::deep_copy(hp,p);

      const std::string cmesKey = "3b CpDTHmes";
      timermap.start( cmesKey );
      Kokkos::deep_copy(meHostPtr,meDevPtr);

      const std::string loopKey = "4a DumpLoop";
      timermap.start(loopKey);
      genrtimes.push_back(genrtime);
      rambtimes.push_back(rambtime);
      wavetimes.push_back(wavetime);

      if (verbose){
        std::cout << "***********************************" << std::endl
                  << "Iteration #" << x+1 << " of " << num_iter << std::endl;

        if (perf) std::cout << "Wave function time: " << wavetime << std::endl;
      }

      if (verbose || perf) {

        for (int d = 0; d < events_per_iter; ++d) {

          // Display matrix elements
          for (int i = 0; i < process.nprocesses; i++) {
            if (verbose)
              std::cout << " Matrix element = "
                        //	 << setiosflags(ios::fixed) << setprecision(17)
                        << meHostPtr(i*1 + d) << " GeV^" << meGeVexponent << std::endl;
            if (perf)
              matrixelementvector.push_back(meHostPtr(i*1 + d));
              weightALL.push_back(h_wgt());
          }

          if (verbose)
            std::cout << std::string(80, '-') << std::endl;
        }
      } else if (!debug) {
        std::cout << ".";
      }

    } // end for num_iter

    // **************************************
    // *** END MAIN LOOP ON #ITERATIONS ***
    // **************************************

    // === STEP 8 ANALYSIS
    // --- 8a Analysis: compute stats after the loop
    const std::string statKey = "8a CompStat";
    timermap.start(statKey);

  double sumgtim = 0;
  double sqsgtim = 0;
  double mingtim = genrtimes[0];
  double maxgtim = genrtimes[0];
  for ( int iiter = 0; iiter < num_iter; ++iiter )
  {
    sumgtim += genrtimes[iiter];
    sqsgtim += genrtimes[iiter]*genrtimes[iiter];
    mingtim = std::min( mingtim, genrtimes[iiter] );
    maxgtim = std::max( maxgtim, genrtimes[iiter] );
  }

  double sumrtim = 0;
  double sqsrtim = 0;
  double minrtim = rambtimes[0];
  double maxrtim = rambtimes[0];
  for ( int iiter = 0; iiter < num_iter; ++iiter )
  {
    sumrtim += rambtimes[iiter];
    sqsrtim += rambtimes[iiter]*rambtimes[iiter];
    minrtim = std::min( minrtim, rambtimes[iiter] );
    maxrtim = std::max( maxrtim, rambtimes[iiter] );
  }

  double sumwtim = 0;
  double sqswtim = 0;
  double minwtim = wavetimes[0];
  double maxwtim = wavetimes[0];
  for ( int iiter = 0; iiter < num_iter; ++iiter )
  {
    sumwtim += wavetimes[iiter];
    sqswtim += wavetimes[iiter]*wavetimes[iiter];
    minwtim = std::min( minwtim, wavetimes[iiter] );
    maxwtim = std::max( maxwtim, wavetimes[iiter] );
  }
  double meanwtim = sumwtim / num_iter;
  double stdwtim = std::sqrt( sqswtim / num_iter - meanwtim * meanwtim );

  int nnan = 0;
  double minelem = matrixelementvector[0];
  double maxelem = matrixelementvector[0];
  double minweig = weightALL[0];
  double maxweig = weightALL[0];
  for ( int ievtALL = 0; ievtALL < nevtALL; ++ievtALL )
  {
    // Compute min/max
    if ( std::isnan( matrixelementvector[ievtALL] ) )
    {
      if ( debug ) // only printed out with "-p -d" (matrixelementALL is not filled without -p)
        std::cout << "WARNING! ME[" << ievtALL << "} is nan" << std::endl;
      nnan++;
      continue;
    }
    minelem = std::min( minelem, (double)matrixelementvector[ievtALL] );
    maxelem = std::max( maxelem, (double)matrixelementvector[ievtALL] );
    minweig = std::min( minweig, (double)weightALL[ievtALL] );
    maxweig = std::max( maxweig, (double)weightALL[ievtALL] );
  }
  double sumelemdiff = 0;
  double sumweigdiff = 0;
  for ( int ievtALL = 0; ievtALL < nevtALL; ++ievtALL )
  {
    // Compute mean from the sum of diff to min
    if ( std::isnan( matrixelementvector[ievtALL] ) ) continue;
    sumelemdiff += ( matrixelementvector[ievtALL] - minelem );
    sumweigdiff += ( weightALL[ievtALL] - minweig );
  }
  double meanelem = minelem + sumelemdiff / ( nevtALL - nnan );
  double meanweig = minweig + sumweigdiff / ( nevtALL - nnan );
  double sqselemdiff = 0;
  double sqsweigdiff = 0;
  for ( int ievtALL = 0; ievtALL < nevtALL; ++ievtALL )
  {
    // Compute stddev from the squared sum of diff to mean
    if ( std::isnan( matrixelementvector[ievtALL] ) ) continue;
    sqselemdiff += std::pow( matrixelementvector[ievtALL] - meanelem, 2 );
    sqsweigdiff += std::pow( weightALL[ievtALL] - meanweig, 2 );
  }
  double stdelem = std::sqrt( sqselemdiff / ( nevtALL - nnan ) );
  double stdweig = std::sqrt( sqsweigdiff / ( nevtALL - nnan ) );

  // === STEP 9 FINALISE
  // --- 9b Dump to screen
  const std::string dumpKey = "9b DumpScrn";
  timermap.start(dumpKey);

  if (!(verbose || debug || perf))
  {
    std::cout << std::endl;
  }

  if (perf)
  {
    std::cout << "***********************************************************************" << std::endl
              << "events_per_iter           = " << events_per_iter << std::endl
              << "num_iterations              = " << num_iter << std::endl
              << "-----------------------------------------------------------------------" << std::endl
              << "Complex type               = KOKKOS::COMPLEX" << std::endl
              << "RanNumb memory layout      = AOSOA[ KOKKOS ]" << std::endl
              << "Momenta memory layout      = AOSOA[ KOKKOS ]" << std::endl
              << "Random number generation   = KOKKOS (C++ code)" << std::endl
              << "-----------------------------------------------------------------------" << std::endl
              << "NumberOfEntries            = " << num_iter << std::endl
              << std::scientific // fixed format: affects all floats (default precision: 6)
              << "TotalTime[Rnd+Rmb+ME] (123)= ( " << sumgtim+sumrtim+sumwtim << std::string(16, ' ') << " )  sec" << std::endl
              << "TotalTime[Rambo+ME]    (23)= ( " << sumrtim+sumwtim << std::string(16, ' ') << " )  sec" << std::endl
              << "TotalTime[RndNumGen]    (1)= ( " << sumgtim << std::string(16, ' ') << " )  sec" << std::endl
              << "TotalTime[Rambo]        (2)= ( " << sumrtim << std::string(16, ' ') << " )  sec" << std::endl
              << "TotalTime[MatrixElems]  (3)= ( " << sumwtim << std::string(16, ' ') << " )  sec" << std::endl
              << "MeanTimeInMatrixElems      = ( " << meanwtim << std::string(16, ' ') << " )  sec" << std::endl
              << "[Min,Max]TimeInMatrixElems = [ " << minwtim
              << " ,  " << maxwtim << " ]  sec" << std::endl
      //<< "StdDevTimeInWaveFuncs      = ( " << stdwtim << std::string(16, ' ') << " )  sec" << std::endl
              << "-----------------------------------------------------------------------" << std::endl
      //<< "ProcessID:                 = " << getpid() << std::endl
      //<< "NProcesses                 = " << process.nprocesses << std::endl
              << "TotalEventsComputed        = " << nevtALL << std::endl
              << "EvtsPerSec[Rnd+Rmb+ME](123)= ( " << nevtALL/(sumgtim+sumrtim+sumwtim)
              << std::string(16, ' ') << " )  sec^-1" << std::endl
              << "EvtsPerSec[Rmb+ME]     (23)= ( " << nevtALL/(sumrtim+sumwtim)
              << std::string(16, ' ') << " )  sec^-1" << std::endl
      //<< "EvtsPerSec[RndNumbGen]   (1)= ( " << nevtALL/sumgtim
      //<< std::string(16, ' ') << " )  sec^-1" << std::endl
      //<< "EvtsPerSec[Rambo]        (2)= ( " << nevtALL/sumrtim
      //<< std::string(16, ' ') << " )  sec^-1" << std::endl
              << "EvtsPerSec[MatrixElems] (3)= ( " << nevtALL/sumwtim
              << std::string(16, ' ') << " )  sec^-1" << std::endl
              << std::defaultfloat; // default format: affects all floats
    std::cout << "***********************************************************************" << std::endl
              << "NumMatrixElements(notNan)  = " << nevtALL - nnan << std::endl
              << std::scientific // fixed format: affects all floats (default precision: 6)
              << "MeanMatrixElemValue        = ( " << meanelem
              << " +- " << stdelem/sqrt(nevtALL - nnan) << " )  GeV^" << meGeVexponent << std::endl // standard error
              << "[Min,Max]MatrixElemValue   = [ " << minelem
              << " ,  " << maxelem << " ]  GeV^" << meGeVexponent << std::endl
              << "StdDevMatrixElemValue      = ( " << stdelem << std::string(16, ' ') << " )  GeV^" << meGeVexponent << std::endl
              << "MeanWeight                 = ( " << meanweig
              << " +- " << stdweig/sqrt(nevtALL - nnan) << " )" << std::endl // standard error
              << "[Min,Max]Weight            = [ " << minweig
              << " ,  " << maxweig << " ]" << std::endl
              << "StdDevWeight               = ( " << stdweig << std::string(16, ' ') << " )" << std::endl
              << std::defaultfloat; // default format: affects all floats
  }

  // *** STOP THE NEW TIMERS ***
  timermap.stop();
  if (perf)
  {
    std::cout << "***********************************************************************" << std::endl;
    timermap.dump();
    std::cout << "***********************************************************************" << std::endl;
  }
  printf("total runtime: %f\n",ktimer.seconds());
  } // end Kokkos View Space
  Kokkos::finalize();
}
