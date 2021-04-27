#ifndef CHECK_H 
#define CHECK_H 1

#include <iostream>
#include <string>
#include <thread>
#include <vector>

#define SEP79 79

// implemented in gheck.cu (compiled with nvcc)
// used by gmain/hmain.cc (compiled using g++)
int gcheck( int argc,
            char **argv,
            std::string& out,
            std::vector<double>& stats,
            const std::string& tag = "",
            const int niter_multiplier = 1 ); // only for the GPU

// implemented in check.cc (compiled with g++)
// used by cmain/hmain.cc (compiled using g++)
int check( int argc,
           char **argv,
           std::string& out,
           std::vector<double>& stats,
           const std::string& tag = "" );

// implemented in check.cc (compiled with g++)
// used by cmain/gmain/hmain.cc (compiled using g++)
#ifndef __CUDACC__
//#warning CUDA not defined
int check_omp_threads( bool debug = false ); // returns the number of OMP threads
#else
//#warning CUDA defined
#endif

// used by check.cc and hmain.cc (compiled using g++)
#ifndef __CUDACC__
inline int nprocall(){ return std::thread::hardware_concurrency(); }
#endif

#endif // CHECK_H
