#include <algorithm> // perf stats
#include <cstring>
#include <iomanip>
#include <iostream>
#include <numeric> // perf stats
#include <unistd.h>
#include <vector>

#include <nvToolsExt.h> 

#include "CPPProcess.h"
#include "random_generator.h"
#include "rambo.h"
#include "Kokkos_Core.hpp"
#include "CalcMean.h"

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

int main(int argc, char **argv) {
  bool verbose = false, debug = false, perf = false;
  int numiter = 0, league_size = 1, team_size = 1;
  std::vector<int> numvec;
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
  if (veclen == 3) {
    league_size = numvec[0];
    team_size = numvec[1];
    numiter = numvec[2];
  } else if (veclen == 1) {
    numiter = numvec[0];
  } else {
    return usage(argv[0]);
  }

  if (numiter == 0)
    return usage(argv[0]);

  Kokkos::initialize(argc, argv);
  if (verbose)
    std::cout << "# iterations: " << numiter << std::endl;
  { // start Kokkos View space
    // Create a process object
    Kokkos::Timer total_time;
    CPPProcess<Kokkos::DefaultExecutionSpace> process(numiter, league_size, team_size);

    // Read param_card and set parameters
    process.initProc("../../Cards/param_card.dat");

    double energy = 1500;

    int meGeVexponent = -(2 * process.nexternal - 8);
    const int events_per_iter = league_size * team_size;

    Kokkos::View<double*,Kokkos::DefaultExecutionSpace> meDevPtr("meDevPtr",events_per_iter*1);
    auto meHostPtr = Kokkos::create_mirror_view(meDevPtr);

    Kokkos::View<double,Kokkos::DefaultExecutionSpace> d_wgt("d_wgt");
    auto h_wgt = Kokkos::create_mirror_view(d_wgt);

    // const int nprocesses = 1; // TODO: hardcoded value
    Kokkos::View<double**,Kokkos::DefaultExecutionSpace> random_numbers(Kokkos::ViewAllocateWithoutInitializing("rns"),events_per_iter,4*(process.nexternal - process.ninitial));
    Kokkos::View<double***,Kokkos::DefaultExecutionSpace> p(Kokkos::ViewAllocateWithoutInitializing("p"),events_per_iter,process.nexternal,4);
    get_initial_momenta(p,process.nexternal,energy,process.cmME,league_size,team_size);
    auto h_p = Kokkos::create_mirror_view(p);

    // init random number generator pool
    auto rand_pool = init_random_generator();

    CalcMean ave_me;
    CalcMean tmr_rand;
    CalcMean tmr_dcpA;
    CalcMean tmr_mom;
    CalcMean tmr_skin;
    CalcMean tmr_dcpB;
    CalcMean tmr_dcpC;
    Kokkos::Timer lptimer;
    for (int x = 0; x < numiter; ++x) {
      // printf("iter %d of %d\n",x,numiter);
      // Get phase space point
      // Kokkos::Timer ptimer;
      nvtxRangePush("fill_random_numbers_2d");
      lptimer.reset();
      fill_random_numbers_2d(random_numbers,events_per_iter,4*(process.nexternal - process.ninitial), rand_pool, league_size, team_size);
      Kokkos::DefaultExecutionSpace().fence();
      tmr_rand.add_value(lptimer.seconds());
      nvtxRangePop();
      
      nvtxRangePush("get_initial_momenta");
      lptimer.reset();
      Kokkos::deep_copy(p,h_p);
      tmr_dcpA.add_value(lptimer.seconds());
      nvtxRangePop();
      
      nvtxRangePush("get_final_momenta");
      lptimer.reset();
      get_final_momenta(process.ninitial, process.nexternal, energy, process.cmME, p, random_numbers, d_wgt, league_size, team_size);
      Kokkos::DefaultExecutionSpace().fence();
      tmr_mom.add_value(lptimer.seconds());
      nvtxRangePop();
      
      nvtxRangePush("sigmaKin");
      lptimer.reset();
      sigmaKin(p, meDevPtr, process.cHel, process.cIPD, process.cIPC, league_size, team_size);//, debug, verbose);
      Kokkos::DefaultExecutionSpace().fence();
      tmr_skin.add_value(lptimer.seconds());
      nvtxRangePop();
      
      nvtxRangePush("copy_momenta_DtoH");
      lptimer.reset();
      Kokkos::deep_copy(h_p,p);
      tmr_dcpB.add_value(lptimer.seconds());
      nvtxRangePop();
      
      nvtxRangePush("copy_matrixEl_DtoH");
      lptimer.reset();
      Kokkos::deep_copy(meHostPtr,meDevPtr);
      tmr_dcpC.add_value(lptimer.seconds());
      nvtxRangePop();

      nvtxRangePush("reporting");
      if (verbose)
        std::cout << "***********************************" << std::endl
                  << "Iteration #" << x+1 << " of " << numiter << std::endl;

      if (verbose || perf) {

        for (int d = 0; d < events_per_iter; ++d) {

          if (verbose) {
            std::cout << "Momenta:" << std::endl;
            for (int i = 0; i < process.nexternal; i++)
              std::cout << std::setw(4) << i + 1
                        << setiosflags(std::ios::scientific) << std::setw(14)
                        << h_p(d,i,0) << setiosflags(std::ios::scientific)
                        << std::setw(14) << h_p(d,i,1)
                        << setiosflags(std::ios::scientific) << std::setw(14)
                        << h_p(d,i,2) << setiosflags(std::ios::scientific)
                        << std::setw(14) << h_p(d,i,3) << std::endl;
            std::cout << std::string(80, '-') << std::endl;
          }

          // Display matrix elements
          for (int i = 0; i < process.nprocesses; i++) {
            if (verbose)
              std::cout << " Matrix element = "
                        //	 << setiosflags(ios::fixed) << setprecision(17)
                        << meHostPtr(i*1 + d) << " GeV^" << meGeVexponent << std::endl;
            if (perf)
              ave_me.add_value(meHostPtr(i*1 + d));
          }

          if (verbose)
            std::cout << std::string(80, '-') << std::endl;
        }
      } else if (!debug) {
        std::cout << ".";
      }
      nvtxRangePop();

    } // end for numiter


    if (!(verbose || debug || perf)) {
      std::cout << std::endl;
    }

    if (perf) {
      

      std::cout << "***********************************" << std::endl
                << "NumIterations         = " << numiter << std::endl
                << "NumThreadsPerBlock    = " << team_size << std::endl
                << "NumBlocksPerGrid      = " << league_size << std::endl
                << "-----------------------------------" << std::endl
                << std::scientific
                << "TotalTimeInWaveFuncs  = " << tmr_skin.sum() << " sec" << std::endl
                << "MeanTimeInWaveFuncs   = " << tmr_skin.mean() << " sec" << std::endl
                << "StdDevTimeInWaveFuncs = " << tmr_skin.sigma() << " sec" << std::endl
                << "MinTimeInWaveFuncs    = " << tmr_skin.min() << " sec" << std::endl
                << "MaxTimeInWaveFuncs    = " << tmr_skin.max() << " sec" << std::endl
                << "-----------------------------------" << std::endl
                << "ProcessID:            = " << getpid() << std::endl
                << "NProcesses            = " << process.nprocesses << std::endl
                << "NumMatrixElements     = " << ave_me.n() << std::endl
                << "MatrixElementsPerSec  = " << ave_me.n()/tmr_skin.sum() << " sec^-1" << std::endl;

      std::cout << "***********************************" << std::endl
                << "NumMatrixElements     = " << ave_me.n() << std::endl
                << std::scientific
                << "MeanMatrixElemValue   = " << ave_me.mean() << " GeV^" << meGeVexponent << std::endl
                << "StdErrMatrixElemValue = " << ave_me.sigma()/sqrt(ave_me.n()) << " GeV^" << meGeVexponent << std::endl
                << "StdDevMatrixElemValue = " << ave_me.sigma() << " GeV^" << meGeVexponent << std::endl
                << "MinMatrixElemValue    = " << ave_me.min() << " GeV^" << meGeVexponent << std::endl
                << "MaxMatrixElemValue    = " << ave_me.max() << " GeV^" << meGeVexponent << std::endl;

      std::cout << "***********************************" << std::endl
                << "fill_random_numbers   = " << tmr_rand.mean() << " +/- " << tmr_rand.sigma() << " seconds\n"
                << "copy momenta          = " << tmr_dcpA.mean() << " +/- " << tmr_dcpA.sigma() << " seconds\n"
                << "get final momenta     = " << tmr_mom.mean()  << " +/- " << tmr_mom.sigma()  << " seconds\n"
                << "sigmaKin              = " << tmr_skin.mean() << " +/- " << tmr_skin.sigma() << " seconds\n"
                << "copy momenta          = " << tmr_dcpB.mean() << " +/- " << tmr_dcpB.sigma() << " seconds\n"
                << "copy matrix_element   = " << tmr_dcpC.mean() << " +/- " << tmr_dcpC.sigma() << " seconds\n";
    }

    printf("total time: %e\n",total_time.seconds());
  } // end Kokkos View Space
  Kokkos::finalize();
}
