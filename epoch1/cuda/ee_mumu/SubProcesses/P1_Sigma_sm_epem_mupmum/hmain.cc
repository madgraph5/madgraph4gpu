#include <iostream>
#include <string>
#include <thread>
#include <vector>

int gcheck( int argc, char **argv, std::string& out, std::vector<double>& stats ); // from gcheck.cu (compiled with nvcc)
int check( int argc, char **argv, std::string& out, std::vector<double>& stats ); // from check.cc (compiled with g++)

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

  if ( gpuStatus != 0 ) return 1;
  if ( cpuStatus != 0 ) return 2;
  return 0;
}
