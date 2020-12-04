#include <iostream>
#include <string>

int check( int argc, char **argv, std::string& out ); // from check.cc (compiled with g++)
int gcheck( int argc, char **argv, std::string& out ); // from gcheck.cu (compiled with nvcc)

// This is built with g++ and linked with objects compiled with g++
int main( int argc, char **argv )
{
  std::string cpuOut;
  int cpuStatus = check( argc, argv, cpuOut );
  std::cout << cpuOut;
  return cpuStatus;
}
