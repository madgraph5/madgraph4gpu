#include "check.h"

#include <iomanip>
#include <thread>

// This is compiled with g++ and linked with objects compiled with nvcc or g++
int main( int argc, char **argv )
{
  std::string gpuOut;
  std::vector<double> gpuStats;
  int gpuStatus;
  //int gpuMult = 1; // GPU processes same #events as CPU, with same random seeds
  int gpuMult = 100; // GPU processes 100x #events as CPU, with different random seeds
  std::thread gpuThread( [&]{ gpuStatus = gcheck( argc, argv, gpuOut, gpuStats, "(GPU) ", gpuMult ); });

  std::string cpuOut;
  std::vector<double> cpuStats;
  int cpuStatus;
  std::thread cpuThread( [&]{ cpuStatus = check( argc, argv, cpuOut, cpuStats, "(CPU) " ); });

  gpuThread.join();
  std::cout << gpuOut;
  //for ( auto& stat : gpuStats ) std::cout << stat << std::endl;

  cpuThread.join();
  std::cout << cpuOut;
  //for ( auto& stat : cpuStats ) std::cout << stat << std::endl;

  int nevtALL = (int)(cpuStats[0]+gpuStats[0]);
  double sumgtim = cpuStats[1]+gpuStats[1];
  double sumrtim = cpuStats[2]+gpuStats[2];
  double sumwtim = cpuStats[3]+gpuStats[3];

  std::string tag = "(HET) ";
  std::cout << "----------------------------------------------------------------------------" << std::endl
            << tag << "TotalEventsComputed        = " << nevtALL << std::endl
            << std::scientific // fixed format: affects all floats (default precision: 6)
            << tag << "TotalTime[Rnd+Rmb+ME] (123)= ( "
            << sumgtim+sumrtim+sumwtim << std::string(16, ' ') << " )  sec" << std::endl
            << tag << "TotalTime[Rambo+ME]    (23)= ( "
            << sumrtim+sumwtim << std::string(16, ' ') << " )  sec" << std::endl
            << tag << "TotalTime[RndNumGen]    (1)= ( " << sumgtim << std::string(16, ' ') << " )  sec" << std::endl
            << tag << "TotalTime[Rambo]        (2)= ( " << sumrtim << std::string(16, ' ') << " )  sec" << std::endl
            << tag << "TotalTime[MatrixElems]  (3)= ( " << sumwtim << std::string(16, ' ') << " )  sec" << std::endl
            << "----------------------------------------------------------------------------" << std::endl
            << tag << "EvtsPerSec[Rnd+Rmb+ME](123)= ( " << nevtALL/(sumgtim+sumrtim+sumwtim)
            << std::string(16, ' ') << " )  sec^-1" << std::endl
            << tag << "EvtsPerSec[Rmb+ME]     (23)= ( " << nevtALL/(sumrtim+sumwtim)
            << std::string(16, ' ') << " )  sec^-1" << std::endl
            << tag << "EvtsPerSec[MatrixElems] (3)= ( " << nevtALL/sumwtim
            << std::string(16, ' ') << " )  sec^-1" << std::endl
            << std::defaultfloat // default format: affects all floats
            << "----------------------------------------------------------------------------" << std::endl;

  if ( gpuStatus != 0 ) return 1;
  if ( cpuStatus != 0 ) return 2;
  return 0;
}
