#include <iostream>
#include <string>

int check( int argc, char **argv, std::string& out ); // from check.cc (compiled with g++)
int gcheck( int argc, char **argv, std::string& out ); // from gcheck.cu (compiled with nvcc)

// This is built with nvcc and linked with objects compiled with nvcc or g++
int main( int argc, char **argv )
{
  std::string gpuOut;
  int gpuStatus = gcheck( argc, argv, gpuOut );
  std::cout << gpuOut;
  return gpuStatus;
}
