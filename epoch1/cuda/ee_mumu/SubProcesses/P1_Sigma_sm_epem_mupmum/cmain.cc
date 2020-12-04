#include "check.h"

// This is built with g++ and linked with objects compiled with g++
int main( int argc, char **argv )
{
  std::string cpuOut;
  std::vector<double> cpuStats;
  int cpuStatus = check( argc, argv, cpuOut, cpuStats, "(CPU) " );
  std::cout << cpuOut;
  //for ( auto& stat : cpuStats ) std::cout << stat << std::endl;
  return cpuStatus;
}
