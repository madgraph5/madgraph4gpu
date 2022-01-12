#ifndef MATRIXELEMENTKERNELDATA_H 
#define MATRIXELEMENTKERNELDATA_H 1

#include "mgOnGpuConfig.h"
#include "mgOnGpuFptypes.h"

#ifdef __CUDACC__
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
  // A struct encapsulating the input and output data needed by matrix element kernels
  struct MatrixElementKernelData
  {
    const fptype* allmomenta; // input: momenta[nevt*npar*4]
    fptype* allMEs;           // output: allMEs[nevt], |M|^2 running_sum_over_helicities
    const int nevt;           // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
    MatrixElementKernelData( const fptype* _allmomenta, fptype* _allMEs, const int _nevt )
      : allmomenta( _allmomenta ), allMEs( _allMEs ), nevt( _nevt ){}
  };
}
#endif // MATRIXELEMENTKERNELDATA_H
