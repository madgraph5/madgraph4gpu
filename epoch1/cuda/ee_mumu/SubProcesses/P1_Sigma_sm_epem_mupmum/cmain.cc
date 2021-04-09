#include "check.h"

// This is compiled using g++, and then linked using g++ with objects compiled using g++
int main( int argc, char **argv )
{
  std::string cpuOut;
  std::vector<double> cpuStats;
  int cpuStatus = check( argc, argv, cpuOut, cpuStats, "(CPU) " );
  std::cout << cpuOut;
  //for ( auto& stat : cpuStats ) std::cout << stat << std::endl;
  return cpuStatus;
}
