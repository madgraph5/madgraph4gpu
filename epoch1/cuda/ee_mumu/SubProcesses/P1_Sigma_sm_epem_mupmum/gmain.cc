#include <iostream>
#include <string>
#include <vector>

int gcheck( int argc, char **argv, std::string& out, std::vector<double>& stats ); // from gcheck.cu (compiled with nvcc)
int check( int argc, char **argv, std::string& out, std::vector<double>& stats ); // from check.cc (compiled with g++)

// This is built with nvcc and linked with objects compiled with nvcc or g++
int main( int argc, char **argv )
{
  std::string gpuOut;
  std::vector<double> gpuStats;
  int gpuStatus = gcheck( argc, argv, gpuOut, gpuStats );
  std::cout << gpuOut;
  //for ( auto& stat : gpuStats ) std::cout << stat << std::endl;
  return gpuStatus;
}
