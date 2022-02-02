#include <algorithm> // perf stats
#include <cstring>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <numeric> // perf stats
#include <unistd.h>
#include <vector>

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

int main(int argc, char **argv) {
  clock_t start = clock(), end;
  bool verbose = false, debug = false, perf = false, json = false;
  int numiter = 0, league_size = 1, team_size = 1;
  std::vector<int> numvec;
  std::vector<float> wavetimes;
  // int jsondate = 0;
  int jsonrun = 0;


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
    else if (is_number(argv[argn]))
      numvec.push_back(atoi(argv[argn]));
    else
      return usage(argv[0]);
  }
  int veclen = numvec.size();
  if (veclen == 3 || veclen == 5) {
    league_size = numvec[0];
    team_size = numvec[1];
    numiter = numvec[2];
    if (veclen == 5){
      // jsondate = numvec[3];
      jsonrun = numvec[4];
    }
  } else if (veclen == 1) {
    numiter = numvec[0];
  } else {
    return usage(argv[0]);
  }

  if (numiter == 0)
    return usage(argv[0]);

  end = clock();
  double cmdline_parse_sec = ((double) (end - start)) / CLOCKS_PER_SEC;
  std::cout << "command line parsing: " << cmdline_parse_sec << " seconds\n";
  start = clock();

  Kokkos::initialize(argc, argv);

  end = clock();
  double kokkos_init_sec = ((double) (end - start)) / CLOCKS_PER_SEC;
  std::cout << "Kokkkos initialize: " << kokkos_init_sec << " seconds\n";
  start = clock();

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
    constexpr double energy = 1500;
    const int meGeVexponent = -(2 * process.nexternal - 8);

    auto time_procInit = lptimer.seconds();
    nvtxRangePop();

    nvtxRangePush("0b_MemAlloc");
    lptimer.reset();

    const int events_per_iter = league_size * team_size;

    Kokkos::View<double*,Kokkos::DefaultExecutionSpace> d_me(Kokkos::ViewAllocateWithoutInitializing("d_me"),events_per_iter*process.nprocesses);
    auto h_me = Kokkos::create_mirror_view(d_me);

    Kokkos::View<double*,Kokkos::DefaultExecutionSpace> d_wgt(Kokkos::ViewAllocateWithoutInitializing("d_wgt"),events_per_iter*process.nprocesses);
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
    CalcMean<float,unsigned int> ave_weight;
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
      sigmaKin(p, d_me, process.cHel, process.cIPD, process.cIPC, iGoodHel, nGoodHel, process.ncomb, league_size, team_size);//, debug, verbose);
      Kokkos::DefaultExecutionSpace().fence();
      tmr_skin.add_value(lptimer.seconds());
      nvtxRangePop();

      nvtxRangePush("3b_CpDTHmes");
      lptimer.reset();
      Kokkos::deep_copy(h_me,d_me);
      tmr_cpyME.add_value(lptimer.seconds());
      nvtxRangePop();

      nvtxRangePush("4a_DumpLoop");
      lptimer.reset();
      if (verbose)
        std::cout << "***********************************" << std::endl
                  << "Iteration #" << x+1 << " of " << numiter << std::endl;

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
                      // << setiosflags(ios::fixed) << setprecision(17)
                      << h_me(i*1 + d) << " GeV^" << meGeVexponent << std::endl;

          ave_me.add_value(h_me(i*1 + d));
          ave_weight.add_value(h_wgt(i*1 + d));

        }

        if (verbose)
          std::cout << std::string(80, '-') << std::endl;
      }

      if (!(verbose || debug || perf))
        std::cout << ".";
      nvtxRangePop();
      tmr_dumploop.add_value(lptimer.seconds());
      tmr_iter.add_value(iter_timer.seconds());
    } // end for numiter
    if (!(verbose || debug || perf))
      std::cout << std::endl;

    nvtxRangePush("8a_9a_DumpStat");
    lptimer.reset();
    int nevtALL = numiter*events_per_iter;
    // timer sums
    double tmr_sum_me = tmr_skin.sum() + tmr_cpyME.sum();
    double tmr_sum_rmb = tmr_momini.sum() + tmr_momfin.sum() + tmr_cpyWgt.sum() + tmr_cpyMom.sum();
    double tmr_sum_rnd = tmr_rand.sum();
    double tmr_sum_rmb_me = tmr_sum_me + tmr_sum_rmb;
    double tmr_sum_rnd_rmb_me = tmr_sum_me + tmr_sum_rmb + tmr_sum_rnd;

    if (perf) {
      printf("**********************************************************************\n");
      printf("NumBlocksPerGrid            = %8d\n",league_size);
      printf("NumThreadsPerBlock          = %8d\n",team_size);
      printf("NumIterations               = %8d\n",numiter);
      printf("----------------------------------------------------------------------\n");
      printf("FP Precision                = DOUBLE\n");
#ifdef THRUST_COMPLEX
      printf("Complex type                = THRUST::COMPLEX\n");
#else
      printf("Complex type                = KOKKOS::COMPLEX\n");
#endif
      printf("Random number generator     = Kokkos Device Side\n");
      printf("----------------------------------------------------------------------\n");
      printf("NumberOfEntries             = %8d\n",numiter);
      printf("TotalTime[Rnd+Rmb+ME] (123) = ( %.6e ) sec\n",tmr_sum_rnd_rmb_me);
      printf("TotalTime[Rambo+ME]    (23) = ( %.6e ) sec\n",tmr_sum_rmb_me);
      printf("TotalTime[RndNumGen]    (1) = ( %.6e ) sec\n",tmr_sum_rnd);
      printf("TotalTime[Rambo]        (2) = ( %.6e ) sec\n",tmr_sum_rmb);
      printf("TotalTime[MatrixElems]  (3) = ( %.6e ) sec\n",tmr_sum_me);
      printf("MeanTimeInMatrixElems       = ( %.6e ) sec\n",tmr_skin.mean()+tmr_cpyME.mean());
      printf("[Min,Max]TimeInMatrixElems  = [ %.6e , %.6e ] sec\n",tmr_skin.min()+tmr_cpyME.min(),tmr_skin.max()+tmr_cpyME.max());

      printf("----------------------------------------------------------------------\n");
      printf("TotalEventsComputed         = %8d\n",nevtALL);
      printf("EvtsPerSec[Rnd+Rmb+ME](123) = ( %.6e ) sec^-1 \n",nevtALL/tmr_sum_rnd_rmb_me);
      printf("EvtsPerSec[Rmb+ME]     (23) = ( %.6e ) sec^-1 \n",nevtALL/tmr_sum_rmb_me);
      printf("EvtsPerSec[MatrixElems] (3) = ( %.6e ) sec^-1 \n",nevtALL/tmr_sum_me);

      printf("**********************************************************************\n");
      printf("NumMatrixElements(notNan)   = %8d\n",nevtALL);
      printf("MeanMatrixElemValue         = ( %.6e +- %.6e ) GeV^%d\n",ave_me.mean(),ave_me.sigma(),meGeVexponent);
      printf("[Min,Max]MatrixElemValue    = [ %.6e , %.6e ]  GeV^%d\n",ave_me.min(),ave_me.max(),meGeVexponent);
      printf("StdDevMatrixElemValue       = ( %.6e ) GeV^%d\n",ave_me.sigma()/sqrt(nevtALL),meGeVexponent);
      printf("MeanWeight                  = ( %.6e +- %.6e ) GeV^%d\n",ave_weight.mean(),ave_weight.sigma(),meGeVexponent);
      printf("[Min,Max]Weight             = [ %.6e , %.6e ]  GeV^%d\n",ave_weight.min(),ave_weight.max(),meGeVexponent);
      printf("StdDevWeight                = ( %.6e ) GeV^%d\n",ave_weight.sigma()/sqrt(nevtALL),meGeVexponent);

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

    if(json){
      std::stringstream json_fn;
      json_fn << "./perf/data/" << league_size << "-" << team_size << "-" << numiter 
              << "-perf-test-run" << jsonrun << ".json";

      std::ofstream fout(json_fn.str());
      fout << "[{\n";
      fout << "  \"NumberOfEntries\": "         << numiter << ",\n";
      fout << "  \"NumThreadsPerBlock\": "      << team_size << ",\n";
      fout << "  \"NumBlocksPerGrid\": "        << league_size << ",\n";
      fout << "  \"FP precision\": \"DOUBLE\",\n";
#ifdef THRUST_COMPLEX
      fout << "  \"Complex type\": \"THRUST::COMPLEX\",\n";
#else
      fout << "  \"Complex type\": \"KOKKOS::COMPLEX\",\n";
#endif
      fout << "  \"TotalTimeInWaveFuncs\": "    << std::scientific << tmr_skin.sum()+tmr_cpyME.sum() << ",\n";
      fout << "  \"MeanTimeInWaveFuncs\": "     << tmr_skin.mean()+tmr_cpyME.mean() << ",\n";
      fout << "  \"StdDevTimeInWaveFuncs\": "   << tmr_skin.sigma()+tmr_cpyME.sigma() << ",\n";
      fout << "  \"TotalEventsComputed\": "     << nevtALL << ",\n";
      fout << "  \"MinTimeInWaveFuncs\": "      << tmr_skin.min()+tmr_cpyME.min() << ",\n";
      fout << "  \"MaxTimeInWaveFuncs\": "      << tmr_skin.max()+tmr_cpyME.max() << ",\n";
      fout << "  \"RamboEventsPerSec\": "       << nevtALL/tmr_sum_rmb << ",\n";
      fout << "  \"MatrixElemEventsPerSec\": "  << nevtALL/tmr_sum_me << ",\n";
      fout << "  \"MeanMatrixElemValue\": "     << ave_me.mean() << ",\n";
      fout << "  \"StdDevMatrixElemValue\": "   << ave_me.sigma() << ",\n";
      fout << "  \"StdErrMatrixElemValue\": "   << ave_me.sigma()/sqrt(nevtALL) << ",\n";
      fout << "  \"MinMatrixElemValue\": "      << ave_me.min() << ",\n";
      fout << "  \"MaxMatrixElemValue\": "      << ave_me.max() << ",\n";
      fout << "  \"MatrixElemUnits\": "         << " \"GeV^" << meGeVexponent << "\",\n";
      fout << "  \"rateUnits\": "               << " \"sec^-1\",\n";
      fout << "  \"periodUnits\": "               << " \"sec\"\n";
      fout << "}]";
      fout.close();
    }

    printf("iteration time        = %10.3f +/- %10.3f seconds\n",tmr_iter.mean(),tmr_iter.sigma());
    printf("total time            = %10.3f seconds\n",total_time.seconds());
  } // end Kokkos View Space
  Kokkos::finalize();
}
