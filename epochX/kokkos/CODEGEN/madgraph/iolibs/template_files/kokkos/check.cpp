#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <string>
#include <unistd.h>

#include "Kokkos_Core.hpp"
#include "CPPProcess.h"
#include "mgKokkosTypes.h"
#include "random_generator.h"
#include "rambo.h"
#include "CalcMean.h"

bool is_number(const char *s) {
  const char *t = s;
  while (*t != '\0' && isdigit(*t))
    ++t;
  return (int)strlen(s) == t - s;
}

int usage(char* argv0, int ret = 1) {
  std::cout << "Usage: " << argv0
            << " [--verbose|-v] [--debug|-d] [--performance|-p] [--json|-j]"
            << " [#league-size #team-size] #iterations" << std::endl << std::endl;
  std::cout << "The number of events per iteration is #thread-teams * #team-size" << std::endl;
  std::cout << "(also in CPU/C++ code, where only the product of these two parameters counts)" << std::endl << std::endl;
  std::cout << "Summary stats are always computed: '-p' and '-j' only control their printout" << std::endl;
  std::cout << "The '-d' flag only controls if nan's emit warnings" << std::endl;
  return ret;
}

// #ifdef __CUDACC__
// template<typename T = fptype>
// struct CudaDevDeleter {
//   void operator()(T* mem) {
//     checkCuda( cudaFree( mem ) );
//   }
// };
// template<typename T = fptype>
// std::unique_ptr<T, CudaDevDeleter<T>> devMakeUnique(std::size_t N) {
//   T* tmp = nullptr;
//   checkCuda( cudaMalloc( &tmp, N * sizeof(T) ) );
//   return std::unique_ptr<T, CudaDevDeleter<T>>{ tmp };
// }
// template<typename T = fptype>
// struct CudaHstDeleter {
//   void operator()(T* mem) {
//     checkCuda( cudaFreeHost( mem ) );
//   }
// };
// template<typename T = fptype>
// std::unique_ptr<T[], CudaHstDeleter<T>> hstMakeUnique(std::size_t N) {
//   T* tmp = nullptr;
//   checkCuda( cudaMallocHost( &tmp, N * sizeof(T) ) );
//   return std::unique_ptr<T[], CudaHstDeleter<T>>{ tmp };
// };
// #else
// template<typename T = fptype>
// std::unique_ptr<T[]> hstMakeUnique(std::size_t N) { return std::unique_ptr<T[]>{ new T[N] }; };
// #endif

int main(int argc, char **argv)
{
  clock_t clock_start = clock(), clock_end = clock();
  // READ COMMAND LINE ARGUMENTS
  bool verbose = false;
  bool debug = false;
  bool perf = false;
  bool json = false;
  int niter = 1, league_size = 1, team_size = 1;
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
    league_size = numvec[0];
    team_size = numvec[1];
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

  clock_end = clock();
  double cmdline_parse_sec =  ((double) (clock_end - clock_start)) / CLOCKS_PER_SEC;
  std::cout << "command line parsing: " << cmdline_parse_sec << " seconds\n";
  clock_start = clock();

  // initialize Kokkos
  Kokkos::initialize(argc, argv);

  clock_end = clock();
  
  double kokkos_init_sec = ((double) (clock_end - clock_start)) / CLOCKS_PER_SEC;
  std::cout << "Kokkkos initialize: " << kokkos_init_sec << " seconds\n";

  // bracket Kokkos Views so that deconstructors are called
  // before Kokkos::finalize()
  {

    // const int neppR = mgKokkos::neppR; // ASA layout: constant at compile-time
    // if ( team_size%neppR != 0 )
    // {
    //   std::cout << "ERROR! #threads/block should be a multiple of neppR=" << neppR << std::endl;
    //   return usage(argv[0]);
    // }

    // const int neppM = mgKokkos::neppM; // ASA layout: constant at compile-time
    // if ( team_size%neppM != 0 )
    // {
    //   std::cout << "ERROR! #threads/block should be a multiple of neppM=" << neppM << std::endl;
    //   return usage(argv[0]);
    // }

    // using mgKokkos::ntpbMAX;
    // if ( team_size > ntpbMAX )
    // {
    //   std::cout << "ERROR! #threads/block should be <= " << ntpbMAX << std::endl;
    //   return usage(argv[0]);
    // }

    const int ndim = league_size * team_size; // number of threads in one GPU grid
    const int nevt = ndim; // number of events in one iteration == number of GPU threads
    const int nevtALL = niter*nevt; // total number of ALL events in all iterations

    if (verbose)
      std::cout << "# iterations: " << niter << std::endl;

    // *** START THE NEW TIMERS ***
    Kokkos::Timer total_timer,loop_timer,section_timer;

    // === STEP 0 - INITIALISE

    // --- 0a. Initialise physics process
    nvtxRangePush("0a_ProcInit");
    section_timer.reset();

    // Create a process object
    Proc::CPPProcess<Kokkos::DefaultExecutionSpace> process(niter, league_size, team_size);

    // Read param_card and set parameters
    process.initProc("../../Cards/param_card.dat");
    const double energy = 1500; // historical default, Ecms = 1500 GeV = 1.5 TeV (above the Z peak)
    //const fptype energy = 91.2; // Ecms = 91.2 GeV (Z peak)
    //const fptype energy = 0.100; // Ecms = 100 MeV (well below the Z peak, pure em scattering)
    const int meGeVexponent = -(2 * process.nexternal - 8);

    auto time_procInit = section_timer.seconds();
    nvtxRangePop();

    // --- 0b. Allocate memory structures
    nvtxRangePush("0b_MemAlloc");
    section_timer.reset();

    // random numbers (device-side only)
    Kokkos::View<double**,Kokkos::DefaultExecutionSpace> devRnarray(Kokkos::ViewAllocateWithoutInitializing("devRnarray"),nevt,4*(process.nexternal - process.ninitial));
    
    // momenta
    Kokkos::View<double***,Kokkos::DefaultExecutionSpace> devMomenta(Kokkos::ViewAllocateWithoutInitializing("devMomenta"),nevt,process.nexternal,4);
    auto hstMomenta = Kokkos::create_mirror_view(devMomenta);

    // matrix elements
    Kokkos::View<double*,Kokkos::DefaultExecutionSpace> devMEs(Kokkos::ViewAllocateWithoutInitializing("devMEs"),nevt*process.nprocesses);
    auto hstMEs = Kokkos::create_mirror_view(devMEs);

    // weights
    Kokkos::View<double*,Kokkos::DefaultExecutionSpace> devWeights(Kokkos::ViewAllocateWithoutInitializing("devWeights"),nevt*process.nprocesses);
    auto hstWeights = Kokkos::create_mirror_view(devWeights);

    // good helicity indices tracking (device-side only)
    Kokkos::View<int*,Kokkos::DefaultExecutionSpace> devNGoodHel("devNGoodHel",1);
    Kokkos::View<int*,Kokkos::DefaultExecutionSpace> devIsGoodHel("devIsGoodHel",process.ncomb);

    auto time_memAlloc = section_timer.seconds();
    nvtxRangePop();

    // --- 0c. Create curand or common generator
    nvtxRangePush("0c_GenCreat");
    section_timer.reset();

    // init random number generator pool
    auto rand_pool = init_random_generator();

    auto time_genCreat = section_timer.seconds();
    nvtxRangePop();


    // **************************************
    // *** START MAIN LOOP ON #ITERATIONS ***
    // **************************************

    // creating statistics counters for metrics of interest
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
    
    for (int iiter = 0; iiter < niter; ++iiter)
    {
      //std::cout << "Iteration #" << iiter+1 << " of " << niter << std::endl;
      
      // time each loop
      loop_timer.reset();
      

      // === STEP 1 OF 3

      nvtxRangePush("1a_1b_1c_RandGen");
      section_timer.reset();

      fill_random_numbers_2d(devRnarray,nevt,4*(process.nexternal - process.ninitial), rand_pool, league_size, team_size);
      Kokkos::DefaultExecutionSpace().fence();
      tmr_rand.add_value(section_timer.seconds());
      nvtxRangePop();

      // === STEP 2 OF 3
      // Fill in particle momenta for each of nevt events on the device

      // --- 2a. Fill in momenta of initial state particles on the device
      nvtxRangePush("2a_RamboIni");
      section_timer.reset();
      get_initial_momenta(devMomenta,process.nexternal,energy,process.cmME,league_size,team_size);
      tmr_momini.add_value(section_timer.seconds());
      nvtxRangePop();
      //std::cout << "Got initial momenta" << std::endl;

      // --- 2b. Fill in momenta of final state particles using the RAMBO algorithm on the device
      // (i.e. map random numbers to final-state particle momenta for each of nevt events)
      nvtxRangePush("2b_RamboFin");
      section_timer.reset();
      get_final_momenta(process.ninitial, process.nexternal, energy, process.cmME, devMomenta, devRnarray, devWeights, league_size, team_size);
      Kokkos::DefaultExecutionSpace().fence();
      tmr_momfin.add_value(section_timer.seconds());
      nvtxRangePop();
      //std::cout << "Got final momenta" << std::endl;

      nvtxRangePush("2c_CpDTHwgt");
      section_timer.reset();
      Kokkos::deep_copy(hstWeights,devWeights);
      Kokkos::DefaultExecutionSpace().fence();
      tmr_cpyWgt.add_value(section_timer.seconds());
      nvtxRangePop();

      nvtxRangePush("2d_CpDTHmom");
      section_timer.reset();
      Kokkos::deep_copy(hstMomenta,devMomenta);
      Kokkos::DefaultExecutionSpace().fence();
      tmr_cpyMom.add_value(section_timer.seconds());
      nvtxRangePop();

      // === STEP 3 OF 3
      // Evaluate matrix elements for all nevt events
      // 0d. (Only on the first iteration) Get good helicities [renamed as 0d: this is initialisation!]
      // 3a. Evaluate MEs on the device
      // 3b. Copy MEs back from device to host

      // --- 0d. SGoodHel
      if ( iiter == 0 )
      {
        nvtxRangePush("0d_SGoodHel");
        section_timer.reset();
        Proc::sigmaKin_setup(devMomenta, process.cHel, process.cIPD, process.cIPC, devIsGoodHel, devNGoodHel, process.ncomb, league_size, team_size);
        Kokkos::DefaultExecutionSpace().fence();
        time_SGoodHel = section_timer.seconds();
        nvtxRangePop();
      }

      // --- 3a. SigmaKin
       nvtxRangePush("3a_SigmaKin");
      section_timer.reset();
      Proc::sigmaKin(devMomenta, devMEs, process.cHel, process.cIPD, process.cIPC, devIsGoodHel, devNGoodHel, process.ncomb, league_size, team_size);//, debug, verbose);
      Kokkos::DefaultExecutionSpace().fence();
      tmr_skin.add_value(section_timer.seconds());
      nvtxRangePop();


      nvtxRangePush("3b_CpDTHmes");
      section_timer.reset();
      Kokkos::deep_copy(hstMEs,devMEs);
      Kokkos::DefaultExecutionSpace().fence();
      tmr_cpyME.add_value(section_timer.seconds());
      nvtxRangePop();


      // === STEP 4 FINALISE LOOP
      // --- 4a Dump within the loop
      nvtxRangePush("4a_DumpLoop");
      section_timer.reset();

      if (verbose)
      {
        std::cout << "***********************************************************************" << std::endl
                  << "Iteration #" << iiter+1 << " of " << niter << std::endl;
        if (perf) std::cout << "Wave function time: " << tmr_skin.mean() + tmr_cpyME.mean() << std::endl;
      }

      for (int ievt = 0; ievt < nevt; ++ievt) // Loop over all events in this iteration
      {
        if (verbose)
        {
          // Display momenta
          std::cout << "Momenta:" << std::endl;
          for (int i = 0; i < process.nexternal; i++)
            std::cout << std::setw(4) << i + 1
                      << setiosflags(std::ios::scientific) << std::setw(14)
                      << hstMomenta(ievt,i,0) << setiosflags(std::ios::scientific)
                      << std::setw(14) << hstMomenta(ievt,i,1)
                      << setiosflags(std::ios::scientific) << std::setw(14)
                      << hstMomenta(ievt,i,2) << setiosflags(std::ios::scientific)
                      << std::setw(14) << hstMomenta(ievt,i,3) << std::endl;
          std::cout << std::string(80, '-') << std::endl;
        }
        // Display matrix elements
        for (int i = 0; i < process.nprocesses; i++) {
          if (verbose)
            std::cout << " Matrix element = "
                      // << setiosflags(ios::fixed) << setprecision(17)
                      << hstMEs(i*1 + ievt) << " GeV^" << meGeVexponent << std::endl;
          
          ave_me.add_value(hstMEs(i*1 + ievt));
          ave_weight.add_value(hstWeights(i*1 + ievt));
        }

        if (verbose)
          std::cout << std::string(80, '-') << std::endl;
      } // end if ievt

      if (!(verbose || debug || perf))
          std::cout << ".";
      nvtxRangePop();
      tmr_dumploop.add_value(section_timer.seconds());
      tmr_iter.add_value(loop_timer.seconds());
    } // end for niter
    if (!(verbose || debug || perf))
      std::cout << std::endl;


    // **************************************
    // *** END MAIN LOOP ON #ITERATIONS ***
    // **************************************

    // === STEP 8 ANALYSIS
    // --- 8a Analysis: compute stats after the loop
    nvtxRangePush("8a_9a_DumpStat");
    section_timer.reset();
    
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
      printf("NumIterations               = %8d\n",niter);
      printf("----------------------------------------------------------------------\n");
      printf("FP Precision                = DOUBLE\n");
      printf("Complex type                = KOKKOS::COMPLEX\n");
      printf("Random number generator     = Kokkos Device Side\n");
      printf("----------------------------------------------------------------------\n");
      printf("NumberOfEntries             = %8d\n",niter);
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

      printf("8a_9a_DumpStat        = %8.6f seconds\n",section_timer.seconds());
      printf("**********************************************************************\n");
    }
    
    if(json){
      std::stringstream json_fn;
      json_fn << "./perf/data/" << league_size << "-" << team_size << "-" << niter 
              << "perf-test-run" << jsonrun << ".json";

      std::ofstream fout(json_fn.str());
      fout << "[{\n";
      fout << "  \"NumberOfEntries\": "         << niter << ",\n";
      fout << "  \"NumThreadsPerBlock\": "      << team_size << ",\n";
      fout << "  \"NumBlocksPerGrid\": "        << league_size << ",\n";
      fout << "  \"FP precision\": \"DOUBLE\",\n";
#ifdef THRUST_COMPLEX
      fout << "  \"Complex type\": \"thrust::complex\",\n";
#else
      fout << "  \"Complex type\": \"Kokkos::complex\",\n";
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

    printf("kokkos init time      = %10.3f\n",kokkos_init_sec);
    printf("iteration time        = %10.3f +/- %10.3f seconds\n",tmr_iter.mean(),tmr_iter.sigma());
    printf("total time            = %10.3f seconds\n",total_timer.seconds());
  } // end Kokkos View Space
  Kokkos::finalize();
}
