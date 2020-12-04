#include "check.h"

#include <iomanip>
#include <thread>

// This is compiled with g++ and linked with objects compiled with nvcc or g++
int main( int argc, char **argv )
{
  std::string gpuOut;
  std::vector<double> gpuStats;
  int gpuStatus;
  std::thread gpuThread( [&]{ gpuStatus = gcheck( argc, argv, gpuOut, gpuStats ); });

  std::string cpuOut;
  std::vector<double> cpuStats;
  int cpuStatus;
  std::thread cpuThread( [&]{ cpuStatus = check( argc, argv, cpuOut, cpuStats ); });

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
  
  std::cout << "-----------------------------------------------------------------------" << std::endl
            << std::scientific // fixed format: affects all floats (default precision: 6)
            << "TotalTime[Rnd+Rmb+ME] (123)= ( " << sumgtim+sumrtim+sumwtim << std::string(16, ' ') << " )  sec" << std::endl
            << "TotalTime[Rambo+ME]    (23)= ( " << sumrtim+sumwtim << std::string(16, ' ') << " )  sec" << std::endl
            << "TotalTime[RndNumGen]    (1)= ( " << sumgtim << std::string(16, ' ') << " )  sec" << std::endl
            << "TotalTime[Rambo]        (2)= ( " << sumrtim << std::string(16, ' ') << " )  sec" << std::endl
            << "TotalTime[MatrixElems]  (3)= ( " << sumwtim << std::string(16, ' ') << " )  sec" << std::endl
            << "-----------------------------------------------------------------------" << std::endl
            << "TotalEventsComputed        = " << nevtALL << std::endl
            << "EvtsPerSec[Rnd+Rmb+ME](123)= ( " << nevtALL/(sumgtim+sumrtim+sumwtim)
            << std::string(16, ' ') << " )  sec^-1" << std::endl
            << "EvtsPerSec[Rmb+ME]     (23)= ( " << nevtALL/(sumrtim+sumwtim)
            << std::string(16, ' ') << " )  sec^-1" << std::endl
            << "EvtsPerSec[MatrixElems] (3)= ( " << nevtALL/sumwtim
            << std::string(16, ' ') << " )  sec^-1" << std::endl
            << std::defaultfloat // default format: affects all floats  
            << "-----------------------------------------------------------------------" << std::endl;

  if ( gpuStatus != 0 ) return 1;
  if ( cpuStatus != 0 ) return 2;
  return 0;
}
