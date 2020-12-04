#include <iostream>
#include <string>
#include <vector>

int gcheck( int argc, char **argv, std::string& out, std::vector<double>& stats ); // from gcheck.cu (compiled with nvcc)
int check( int argc, char **argv, std::string& out, std::vector<double>& stats ); // from check.cc (compiled with g++)

// This is built with g++ and linked with objects compiled with g++
int main( int argc, char **argv )
{
  std::string cpuOut;
  std::vector<double> cpuStats;
  int cpuStatus = check( argc, argv, cpuOut, cpuStats );
  std::cout << cpuOut;
  //for ( auto& stat : cpuStats ) std::cout << stat << std::endl;
  return cpuStatus;
}
