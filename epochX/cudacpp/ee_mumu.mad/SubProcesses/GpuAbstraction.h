#ifndef MG5AMC_GPUABSTRACTION_H
#define MG5AMC_GPUABSTRACTION_H 1

#include <cassert>

/*
  ToDo:
      * Fix rpath in makefile when compiling with HIP
      * Fix warnings with improper hip function return code handling
*/


#ifdef __CUDACC__
  #define MGONGPUCPP_CUDACC 1
#endif

#ifdef __HIPCC__
  #include "hip/hip_runtime.h"
  #define MGONGPUCPP_HIPCC 1
#endif

#ifdef MGONGPUCPP_CUDACC

  // Defines correct compiler
  #define MGONGPUCPP_GPUIMPL __CUDACC__

  //--------------------------------------------------------------------------

  #define gpuError_t cudaError_t
  #define gpuPeekAtLastError cudaPeekAtLastError
  #define gpuGetErrorString cudaGetErrorString
  #define gpuSuccess cudaSuccess

  #define gpuMallocHost(ptr, size) checkGpu( cudaMallocHost(ptr, size) )
  #define gpuMalloc(ptr, size) checkGpu( cudaMalloc(ptr, size) )

  #define gpuMemcpy(dstData, srcData, srcBytes, func) checkGpu( cudaMemcpy(dstData, srcData, srcBytes, func) )
  #define gpuMemcpyHostToDevice cudaMemcpyHostToDevice
  #define gpuMemcpyDeviceToHost cudaMemcpyDeviceToHost
  #define gpuMemcpyToSymbol(type1, type2, size) checkGpu( cudaMemcpyToSymbol(type1, type2, size) )

  #define gpuFree(ptr) checkGpu( cudaFree(ptr) )
  #define gpuFreeHost(ptr) checkGpu( cudaFreeHost(ptr) )

  #define gpuSetDevice cudaSetDevice
  #define gpuDeviceSynchronize cudaDeviceSynchronize
  #define gpuDeviceReset cudaDeviceReset

  #define gpuLaunchKernel( kernel, blocks, threads, ...)                    kernel<<<blocks, threads>>> (__VA_ARGS__)
  #define gpuLaunchKernelSharedMem(kernel, blocks, threads, sharedMem, ...) kernel<<<blocks, threads, sharedMem>>>(__VA_ARGS__)

//--------------------------------------------------------------------------

#elif defined MGONGPUCPP_HIPCC

  // Defines correct compiler
  #define MGONGPUCPP_GPUIMPL __HCC__

  //--------------------------------------------------------------------------

  #define gpuError_t hipError_t
  #define gpuPeekAtLastError() checkGpu( hipPeekAtLastError() )
  #define gpuGetErrorString() checkGpu( hipGetErrorString() )
  #define gpuSuccess hipSuccess

  #define gpuMallocHost(ptr, size) checkGpu( hipHostMalloc(ptr, size) ) // HostMalloc better
  #define gpuMalloc(ptr, size) checkGpu( hipMalloc(ptr, size) )

  #define gpuMemcpy(dstData, srcData, srcBytes, func) checkGpu( hipMemcpy(dstData, srcData, srcBytes, func) )
  #define gpuMemcpyHostToDevice hipMemcpyHostToDevice
  #define gpuMemcpyDeviceToHost hipMemcpyDeviceToHost
  #define gpuMemcpyToSymbol(type1, type2, size) checkGpu( hipMemcpyToSymbol(type1, type2, size) )

  #define gpuFree(ptr) checkGpu( hipFree(ptr) )
  #define gpuFreeHost(ptr) checkGpu( hipHostFree(ptr) )

  #define gpuSetDevice hipSetDevice
  #define gpuDeviceSynchronize() checkGpu( hipDeviceSynchronize() )
  #define gpuDeviceReset() checkGpu( hipDeviceReset() )

  #define gpuLaunchKernel( kernel, blocks, threads, ...)                    kernel<<<blocks, threads>>> (__VA_ARGS__)
  #define gpuLaunchKernelSharedMem(kernel, blocks, threads, sharedMem, ...) kernel<<<blocks, threads, sharedMem>>>(__VA_ARGS__)

#endif

#endif // MG5AMC_GPUABSTRACTION_H