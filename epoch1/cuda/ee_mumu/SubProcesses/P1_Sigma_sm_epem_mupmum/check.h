#ifndef CHECK_H 
#define CHECK_H 1

#include <iostream>
#include <string>
#include <vector>

// from gcheck.cu (compiled with nvcc)
int gcheck( int argc,
            char **argv,
            std::string& out,
            std::vector<double>& stats,
            const std::string& tag = "",
            const int niter_multiplier = 1 ); // only for the GPU

// from check.cc (compiled with g++)
int check( int argc,
           char **argv,
           std::string& out,
           std::vector<double>& stats,
           const std::string& tag = "" );

#endif // CHECK_H
