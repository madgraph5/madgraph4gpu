#ifndef CHECK_H 
#define CHECK_H 1

#include <iostream>
#include <string>
#include <vector>

int gcheck( int argc, char **argv, std::string& out, std::vector<double>& stats ); // from gcheck.cu (compiled with nvcc)

int check( int argc, char **argv, std::string& out, std::vector<double>& stats ); // from check.cc (compiled with g++)

#endif // CHECK_H
