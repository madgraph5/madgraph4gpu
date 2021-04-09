#include "check.h"

#include <iomanip>
#include <thread>

void dumptime( const std::string& tag,
               const int nevtALL,
               const double sumgtim,
               const double sumrtim,
               const double sumwtim )
{
  std::cout << "-----------------------------------------------------------------------------" << std::endl
            << tag << "TotalEventsComputed         = " << nevtALL << std::endl
            << std::scientific // fixed format: affects all floats (default precision: 6)
            << tag << "TotalTime[Rnd+Rmb+ME] (123) = ( "
            << sumgtim+sumrtim+sumwtim << std::string(16, ' ') << " )  sec" << std::endl
            << tag << "TotalTime[Rambo+ME]    (23) = ( "
            << sumrtim+sumwtim << std::string(16, ' ') << " )  sec" << std::endl
            << tag << "TotalTime[RndNumGen]    (1) = ( " << sumgtim << std::string(16, ' ') << " )  sec" << std::endl
            << tag << "TotalTime[Rambo]        (2) = ( " << sumrtim << std::string(16, ' ') << " )  sec" << std::endl
            << tag << "TotalTime[MatrixElems]  (3) = ( " << sumwtim << std::string(16, ' ') << " )  sec" << std::endl
            << std::defaultfloat // default format: affects all floats
            << "-----------------------------------------------------------------------------" << std::endl;
}

void dumptput( const std::string& tag,
               const int nevtALL,
               const double tputgrw,
               const double tputrw,
               const double tputw,
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
            << std::defaultfloat // default format: affects all floats
            << "*****************************************************************************" << std::endl;
}

// This is compiled using g++, and then linked using nvcc with objects compiled using nvcc or g++
int main( int argc, char **argv )
{
  std::string gpuOut;
  std::vector<double> gpuStats;
  int gpuStatus;
  //int gpuMult = 1; // GPU processes same #events as CPU, with same random seeds
  int gpuMult = 72; // GPU processes 72x #events as CPU, with different random seeds
  int nthreadsomp = check_omp_threads(); // this is always > 0
  gpuMult /= nthreadsomp; // reduce the number of GPU iterations if several OMP threads are used on the CPU
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

  int gpuNevtALL = (int)(gpuStats[0]);
  double gpuSumgtim = gpuStats[1];
  double gpuSumrtim = gpuStats[2];
  double gpuSumwtim = gpuStats[3];

  dumptime( gpuTag, gpuNevtALL, gpuSumgtim, gpuSumrtim, gpuSumwtim );
  double gpuTputgrw = gpuNevtALL/(gpuSumgtim+gpuSumrtim+gpuSumwtim);
  double gpuTputrw  = gpuNevtALL/(gpuSumrtim+gpuSumwtim);
  double gpuTputw   = gpuNevtALL/(gpuSumwtim);
  dumptput( gpuTag, gpuNevtALL, gpuTputgrw, gpuTputrw, gpuTputw, 0 );

  int cpuNevtALL = (int)(cpuStats[0]);
  double cpuSumgtim = cpuStats[1];
  double cpuSumrtim = cpuStats[2];
  double cpuSumwtim = cpuStats[3];

  dumptime( cpuTag, cpuNevtALL, cpuSumgtim, cpuSumrtim, cpuSumwtim );
  double cpuTputgrw = cpuNevtALL/(cpuSumgtim+cpuSumrtim+cpuSumwtim);
  double cpuTputrw  = cpuNevtALL/(cpuSumrtim+cpuSumwtim);
  double cpuTputw   = cpuNevtALL/(cpuSumwtim);
  dumptput( cpuTag, cpuNevtALL, cpuTputgrw, cpuTputrw, cpuTputw, nthreadsomp );

  std::string hetTag = "(HET) ";
  int hetNevtALL = gpuNevtALL+cpuNevtALL;
  double hetTputgrw = gpuTputgrw+cpuTputgrw;
  double hetTputrw  = gpuTputrw+cpuTputrw;
  double hetTputw   = gpuTputw+cpuTputw;
  dumptput( hetTag, hetNevtALL, hetTputgrw, hetTputrw, hetTputw, 0 );

  if ( gpuStatus != 0 ) return 1;
  if ( cpuStatus != 0 ) return 2;
  return 0;
}
