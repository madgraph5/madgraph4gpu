#ifndef CHECK_H 
#define CHECK_H 1

#include <iostream>
#include <string>
#include <vector>

// from gcheck.cu (compiled with nvcc)
int gcheck( int argc,
            char **argv,
            std::string& out,
            std::vector<double>& stats );

// from check.cc (compiled with g++)
int check( int argc,
           char **argv,
           std::string& out,
           std::vector<double>& stats );

#endif // CHECK_H
