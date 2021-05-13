#include <algorithm> // perf stats
#include <cstring>
#include <iomanip>
#include <iostream>
#include <numeric> // perf stats
#include <unistd.h>
#include <vector>
#include <fstream>

#include "CPPProcess.h"
#include "random_generator.h"
#include "rambo.h"
#include "Kokkos_Core.hpp"
#include "CalcMean.h"

#ifdef __CUDACC__
#include <nvToolsExt.h> 
#else

void nvtxRangePush(const char* text){
  return;
}

void nvtxRangePop(void){
  return;
}

#endif


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
    Kokkos::Timer lptimer;
    
    nvtxRangePush("0a_ProcInit");
    lptimer.reset();
    CPPProcess<Kokkos::DefaultExecutionSpace> process(numiter, league_size, team_size);

    // Read param_card and set parameters
    process.initProc("../../Cards/param_card.dat");

    double energy = 1500;

    int meGeVexponent = -(2 * process.nexternal - 8);

    auto time_procInit = lptimer.seconds();
    nvtxRangePop();

    nvtxRangePush("0b_MemAlloc");
    lptimer.reset();

    const int events_per_iter = league_size * team_size;

    Kokkos::View<double*,Kokkos::DefaultExecutionSpace> meDevPtr(Kokkos::ViewAllocateWithoutInitializing("meDevPtr"),events_per_iter*1);
    auto meHostPtr = Kokkos::create_mirror_view(meDevPtr);

    Kokkos::View<double,Kokkos::DefaultExecutionSpace> d_wgt(Kokkos::ViewAllocateWithoutInitializing("d_wgt"));
    auto h_wgt = Kokkos::create_mirror_view(d_wgt);

    // const int nprocesses = 1; // TODO: hardcoded value
    Kokkos::View<double**,Kokkos::DefaultExecutionSpace> random_numbers(Kokkos::ViewAllocateWithoutInitializing("rns"),events_per_iter,4*(process.nexternal - process.ninitial));
    
    Kokkos::View<double***,Kokkos::DefaultExecutionSpace> p(Kokkos::ViewAllocateWithoutInitializing("p"),events_per_iter,process.nexternal,4);
    auto h_p = Kokkos::create_mirror_view(p);
    
    Kokkos::View<int*,Kokkos::DefaultExecutionSpace> nGoodHel("nGoodHel",1);
    Kokkos::View<int*,Kokkos::DefaultExecutionSpace> iGoodHel("iGoodHel",process.ncomb);

    auto time_memAlloc = lptimer.seconds();
    nvtxRangePop();

    nvtxRangePush("0c_GenCreat");
    lptimer.reset();
    
    // init random number generator pool
    auto rand_pool = init_random_generator();

    auto time_genCreat = lptimer.seconds();
    nvtxRangePop();

    CalcMean<float,unsigned int> ave_me;
    CalcMean<float,unsigned int> tmr_rand;
    CalcMean<float,unsigned int> tmr_momini;
    CalcMean<float,unsigned int> tmr_momfin;
    CalcMean<float,unsigned int> tmr_cpyMom;
    CalcMean<float,unsigned int> tmr_cpyWgt;
    CalcMean<float,unsigned int> tmr_skin;
    CalcMean<float,unsigned int> tmr_cpyME;
    CalcMean<float,unsigned int> tmr_dumploop;
    CalcMean<float,unsigned int> tmr_iter;
    float time_SGoodHel = 0;
    for (int x = 0; x < numiter; ++x) {
      // printf("iter %d of %d\n",x,numiter);
      // Get phase space point
      Kokkos::Timer iter_timer;

      nvtxRangePush("1a_1b_1c_GenSeed");
      lptimer.reset();
      fill_random_numbers_2d(random_numbers,events_per_iter,4*(process.nexternal - process.ninitial), rand_pool, league_size, team_size);
      Kokkos::DefaultExecutionSpace().fence();
      tmr_rand.add_value(lptimer.seconds());
      nvtxRangePop();
      
      nvtxRangePush("2a_RamboIni");
      lptimer.reset();
      get_initial_momenta(p,process.nexternal,energy,process.cmME,league_size,team_size);
      tmr_momini.add_value(lptimer.seconds());
      nvtxRangePop();
      
      nvtxRangePush("2b_RamboFin");
      lptimer.reset();
      get_final_momenta(process.ninitial, process.nexternal, energy, process.cmME, p, random_numbers, d_wgt, league_size, team_size);
      Kokkos::DefaultExecutionSpace().fence();
      tmr_momfin.add_value(lptimer.seconds());
      nvtxRangePop();

      nvtxRangePush("2c_CpDTHwgt");
      lptimer.reset();
      Kokkos::deep_copy(h_wgt,d_wgt);
      tmr_cpyWgt.add_value(lptimer.seconds());
      nvtxRangePop();

      nvtxRangePush("2d_CpDTHmom");
      lptimer.reset();
      Kokkos::deep_copy(h_p,p);
      tmr_cpyMom.add_value(lptimer.seconds());
      nvtxRangePop();
      
      if(x == 0){
        nvtxRangePush("0d_SGoodHel");
        lptimer.reset();
        sigmaKin_setup(p, process.cHel, process.cIPD, process.cIPC, iGoodHel, nGoodHel, process.ncomb, league_size, team_size);
        time_SGoodHel = lptimer.seconds();
        nvtxRangePop();
      }
      
      nvtxRangePush("3a_SigmaKin");
      lptimer.reset();
      sigmaKin(p, meDevPtr, process.cHel, process.cIPD, process.cIPC, iGoodHel, nGoodHel, process.ncomb, league_size, team_size);//, debug, verbose);
      Kokkos::DefaultExecutionSpace().fence();
      tmr_skin.add_value(lptimer.seconds());
      nvtxRangePop();

      nvtxRangePush("3b_CpDTHmes");
      lptimer.reset();
      Kokkos::deep_copy(meHostPtr,meDevPtr);
      tmr_cpyME.add_value(lptimer.seconds());
      nvtxRangePop();

      nvtxRangePush("4a_DumpLoop");
      lptimer.reset();
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
      tmr_dumploop.add_value(lptimer.seconds());
      tmr_iter.add_value(iter_timer.seconds());
    } // end for numiter


    if (!(verbose || debug || perf)) {
      std::cout << std::endl;
    }
    nvtxRangePush("8a_9a_DumpStat");
    lptimer.reset();
    std::ofstream fout("output_data.json");
    fout << "{" << std::endl;
    if (perf) {
      

      printf("**********************************************************************\n");
      printf("NumIterations               = %8d\n",numiter);
      printf("NumThreadsPerBlock          = %8d\n",team_size);
      printf("NumBlocksPerGrid            = %8d\n",league_size);
      printf("----------------------------------------------------------------------\n");
      printf("FP Precision                = DOUBLE\n");
      printf("Complex type                = KOKKOS::COMPLEX\n");
      printf("Random number generator     = Kokkos Device Side\n");
      printf("----------------------------------------------------------------------\n");
      printf("TotalTimeInWaveFuncs        = %8.6f sec\n",tmr_skin.sum());
      printf("MeanTimeInWaveFuncs         = %8.6f sec\n",tmr_skin.mean());
      printf("StdDevTimeInWaveFuncs       = %8.6f sec\n",tmr_skin.sigma());
      printf("MinTimeInWaveFuncs          = %8.6f sec\n",tmr_skin.min());
      printf("MaxTimeInWaveFuncs          = %8.6f sec\n",tmr_skin.max());
      printf("----------------------------------------------------------------------\n");
      printf("NProcesses                  =  %8d \n",process.nprocesses);
      printf("NumMatrixElements           =  %8d \n",ave_me.n());
      auto me_per_sec = ave_me.n() / tmr_skin.sum();
      printf("EventsPerSec[MatrixEls]     =  %8.6e / sec\n",me_per_sec);
      fout << "\"me_per_sec\": " << me_per_sec << "," << std::endl;

      printf("**********************************************************************\n");
      printf("NumMatrixElements           =  %8d \n",ave_me.n());
      fout << "\"num_me\": " << ave_me.n() << "," << std::endl;
      printf("MeanMatrixElemValue         = (%8.6e +/- %8.6e ) GeV^%d \n",ave_me.mean(),ave_me.sigma()/sqrt(ave_me.n()),meGeVexponent);
      fout << "\"mean_me\": " << ave_me.mean() << "," << std::endl;
      fout << "\"me_error\": " << ave_me.sigma()/sqrt(ave_me.n()) << "," << std::endl;
      printf("[Min,Max]MatrixElemValue    = (%8.6e +/- %8.6e ) GeV^%d \n",ave_me.min(),ave_me.max(),meGeVexponent);
      printf("StdDevMatrixElemValue       =  %8.6e GeV^%d \n",ave_me.sigma(),meGeVexponent);
      fout << "\"me_stdev\": " << ave_me.sigma() << "," << std::endl;

      printf("**********************************************************************\n");
      printf("0a_ProcInit           = %8.6f seconds\n",time_procInit);
      printf("0b_MemAlloc           = %8.6f seconds\n",time_memAlloc);
      printf("0c_GenCreat           = %8.6f seconds\n",time_genCreat);
      printf("0d_SGoodHel           = %8.6f seconds\n",time_SGoodHel);

      printf("1a_1b_1c_GenSeed      = %8.6f +/- %8.6f seconds\n",tmr_rand.mean(),tmr_rand.sigma());
      printf("2a_RamboIni           = %8.6f +/- %8.6f seconds\n",tmr_momini.mean(),tmr_momini.sigma());
      printf("2b_RamboFin           = %8.6f +/- %8.6f seconds\n",tmr_momfin.mean(),tmr_momfin.sigma());
      printf("2c_CpDTHwgt           = %8.6f +/- %8.6f seconds\n",tmr_cpyWgt.mean(),tmr_cpyWgt.sigma());
      printf("2d_CpDTHmom           = %8.6f +/- %8.6f seconds\n",tmr_cpyMom.mean(),tmr_cpyMom.sigma());
      printf("3a_SigmaKin           = %8.6f +/- %8.6f seconds\n",tmr_skin.mean(),tmr_skin.sigma());
      printf("3b_CpDTHmes           = %8.6f +/- %8.6f seconds\n",tmr_cpyME.mean(),tmr_cpyME.sigma());
      printf("4a_DumpLoop           = %8.6f +/- %8.6f seconds\n",tmr_dumploop.mean(),tmr_dumploop.sigma());

      printf("8a_9a_DumpStat        = %8.6f seconds\n",lptimer.seconds());
      printf("**********************************************************************\n");

    }
    printf("iteration time        = %10.3f +/- %10.3f seconds\n",tmr_iter.mean(),tmr_iter.sigma());
    printf("total time            = %10.3f seconds\n",total_time.seconds());
    fout << "\"total_time\": " << total_time.seconds() << std::endl << "}\n";
  } // end Kokkos View Space
  Kokkos::finalize();
}
