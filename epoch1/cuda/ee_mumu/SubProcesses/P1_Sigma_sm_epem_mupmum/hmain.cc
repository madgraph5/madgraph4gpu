#include <iostream>
#include <string>

int check( int argc, char **argv, std::string& out ); // from check.cc (compiled with g++)
int gcheck( int argc, char **argv, std::string& out ); // from gcheck.cu (compiled with nvcc)

// This is compiled with g++ and linked with objects compiled with nvcc or g++
int main( int argc, char **argv )
{
  std::string gpuOut;
  int gpuStatus = gcheck( argc, argv, gpuOut );
  std::cout << gpuOut;
  if ( gpuStatus != 0 ) return 1;

  std::string cpuOut;
  int cpuStatus = check( argc, argv, cpuOut );
  std::cout << cpuOut;
  if ( cpuStatus != 0 ) return 2;

  return 0;
}
