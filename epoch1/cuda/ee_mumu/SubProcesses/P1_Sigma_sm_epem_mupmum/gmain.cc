#include "check.h"

// This is compiled using g++, and then linked using nvcc with objects compiled using nvcc or g++
int main( int argc, char **argv )
{
  std::string gpuOut;
  std::vector<double> gpuStats;
  int gpuStatus = gcheck( argc, argv, gpuOut, gpuStats, "(GPU) " );
  std::cout << gpuOut;
  //for ( auto& stat : gpuStats ) std::cout << stat << std::endl;
  return gpuStatus;
}
