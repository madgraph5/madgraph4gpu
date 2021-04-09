#include "check.h"

// This is built with nvcc and linked with objects compiled with nvcc or g++
int main( int argc, char **argv )
{
  std::string gpuOut;
  std::vector<double> gpuStats;
  int gpuStatus = gcheck( argc, argv, gpuOut, gpuStats, "(GPU) " );
  std::cout << gpuOut;
  //for ( auto& stat : gpuStats ) std::cout << stat << std::endl;
  return gpuStatus;
}
