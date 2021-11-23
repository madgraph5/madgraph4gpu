#include "check.h"

#include <iomanip>
#include <thread>

void dumptime( const std::string& tag,
               const int nevtALL,
               const double sumgtim,
               const double sumrtim,
               const double sumwtim,
               const double sumw3atim )
{
  std::cout << std::string(SEP79, '-') << std::endl
            << tag << "TotalEventsComputed         = " << nevtALL << std::endl
            << std::scientific // fixed format: affects all floats (default precision: 6)
            << tag << "TotalTime[Rnd+Rmb+ME] (123) = ( "
            << sumgtim+sumrtim+sumwtim << std::string(16, ' ') << " )  sec" << std::endl
            << tag << "TotalTime[Rambo+ME]    (23) = ( "
            << sumrtim+sumwtim << std::string(16, ' ') << " )  sec" << std::endl
            << tag << "TotalTime[RndNumGen]    (1) = ( " << sumgtim << std::string(16, ' ') << " )  sec" << std::endl
            << tag << "TotalTime[Rambo]        (2) = ( " << sumrtim << std::string(16, ' ') << " )  sec" << std::endl
            << tag << "TotalTime[MatrixElems]  (3) = ( " << sumwtim << std::string(16, ' ') << " )  sec" << std::endl
            << tag << "TotalTime[MECalcOnly]  (3a) = ( " << sumw3atim << std::string(16, ' ') << " )  sec" << std::endl
            << std::defaultfloat // default format: affects all floats
            << std::string(SEP79, '-') << std::endl;
}

void dumptput( const std::string& tag,
               const int nevtALL,
               const double tputgrw,
               const double tputrw,
               const double tputw,
               const double tputw3a,
               const int nthreadsomp ) // >0 for CPU, 0 for GPU and HET
{
  if ( nthreadsomp > 0 )
  {
    std::cout << tag << "OMP threads / `nproc --all` = " << nthreadsomp << " / " << nprocall() << std::endl;
  }
  std::cout << tag << "TotalEventsComputed         = " << nevtALL << std::endl
            << std::scientific // fixed format: affects all floats (default precision: 6)
            << tag << "EvtsPerSec[Rnd+Rmb+ME](123) = ( " << tputgrw
            << std::string(16, ' ') << " )  sec^-1" << std::endl
            << tag << "EvtsPerSec[Rmb+ME]     (23) = ( " << tputrw
            << std::string(16, ' ') << " )  sec^-1" << std::endl
            << tag << "EvtsPerSec[MatrixElems] (3) = ( " << tputw
            << std::string(16, ' ') << " )  sec^-1" << std::endl
            << tag << "EvtsPerSec[MECalcOnly] (3a) = ( " << tputw3a
            << std::string(16, ' ') << " )  sec^-1" << std::endl
            << std::defaultfloat // default format: affects all floats
            << std::string(SEP79, '*') << std::endl;
}

// This is compiled using g++, and then linked using nvcc with objects compiled using nvcc or g++
int main( int argc, char **argv )
{
#ifdef _OPENMP
  int nthreadsomp = check_omp_threads(); // this is always > 0
#else
  int nthreadsomp = 1; // this is always > 0
#endif
  std::string gpuOut;
  std::vector<double> gpuStats;
  int gpuStatus;
  //int gpuMult = 1; // GPU processes same #events as CPU, with same random seeds
  int gpuMult = 192;
  //int gpuMult = 144;
  //int gpuMult = 72;
  if ( gpuMult > 1 )
  {
    // GPU processes processes gpuMult x #events as CPU, with different random seeds
#if defined __AVX512F__
    // Reduce the number of GPU iterations if SIMD is used on the CPU
    gpuMult /= 4;
#elif defined __AVX2__
    // Reduce the number of GPU iterations if SIMD is used on the CPU
    gpuMult /= 4;
#elif defined __SSE4_2__
    // Reduce the number of GPU iterations if SIMD is used on the CPU
    gpuMult /= 2;
#endif
    // Reduce the number of GPU iterations if several OMP threads are used on the CPU
    gpuMult /= nthreadsomp;
  }
  std::string gpuTag = "(GPU) ";
  std::thread gpuThread( [&]{ gpuStatus = gcheck( argc, argv, gpuOut, gpuStats, gpuTag, gpuMult ); });

  std::string cpuOut;
  std::vector<double> cpuStats;
  int cpuStatus;
  std::string cpuTag = "(CPU) ";
  std::thread cpuThread( [&]{ cpuStatus = check( argc, argv, cpuOut, cpuStats, cpuTag ); });

  gpuThread.join();
  std::cout << gpuOut;

  cpuThread.join();
  std::cout << cpuOut;

  if ( !gpuStats.empty() && !cpuStats.empty() )
  {
    int gpuNevtALL = (int)(gpuStats[0]);
    double gpuSumgtim = gpuStats[1];
    double gpuSumrtim = gpuStats[2];
    double gpuSumwtim = gpuStats[3];
    double gpuSumw3atim = gpuStats[4];
    
    int cpuNevtALL = (int)(cpuStats[0]);
    double cpuSumgtim = cpuStats[1];
    double cpuSumrtim = cpuStats[2];
    double cpuSumwtim = cpuStats[3];
    double cpuSumw3atim = cpuStats[4];

    dumptime( gpuTag, gpuNevtALL, gpuSumgtim, gpuSumrtim, gpuSumwtim, gpuSumw3atim );
    double gpuTputgrw = gpuNevtALL/(gpuSumgtim+gpuSumrtim+gpuSumwtim);
    double gpuTputrw  = gpuNevtALL/(gpuSumrtim+gpuSumwtim);
    double gpuTputw   = gpuNevtALL/(gpuSumwtim);
    double gpuTputw3a = gpuNevtALL/(gpuSumw3atim);
    dumptput( gpuTag, gpuNevtALL, gpuTputgrw, gpuTputrw, gpuTputw, gpuTputw3a, 0 );

    dumptime( cpuTag, cpuNevtALL, cpuSumgtim, cpuSumrtim, cpuSumwtim, cpuSumw3atim );
    double cpuTputgrw = cpuNevtALL/(cpuSumgtim+cpuSumrtim+cpuSumwtim);
    double cpuTputrw  = cpuNevtALL/(cpuSumrtim+cpuSumwtim);
    double cpuTputw   = cpuNevtALL/(cpuSumwtim);
    double cpuTputw3a = cpuNevtALL/(cpuSumw3atim);
    dumptput( cpuTag, cpuNevtALL, cpuTputgrw, cpuTputrw, cpuTputw, cpuTputw3a, nthreadsomp );
  
    std::string hetTag = "(HET) ";
    int hetNevtALL = gpuNevtALL+cpuNevtALL;
    double hetTputgrw = gpuTputgrw+cpuTputgrw;
    double hetTputrw  = gpuTputrw+cpuTputrw;
    double hetTputw   = gpuTputw+cpuTputw;
    double hetTputw3a = gpuTputw3a+cpuTputw3a;
    dumptput( hetTag, hetNevtALL, hetTputgrw, hetTputrw, hetTputw, hetTputw3a, 0 );
  }

  if ( gpuStatus != 0 ) return 1;
  if ( cpuStatus != 0 ) return 2;
  return 0;
}
