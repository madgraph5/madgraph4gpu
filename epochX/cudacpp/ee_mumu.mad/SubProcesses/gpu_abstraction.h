// gpu_abstraction.h
#pragma once

#if defined(__HIP_PLATFORM_HCC__)
    #include <hip/hip_runtime.h>

    #define gpuError_t hipError_t
    #define gpuMalloc hipMalloc
    #define gpuFree hipFree
    #define gpuMemcpy hipMemcpy
    // Add other necessary abstractions

#elif defined(__CUDACC__)

    #define gpuError_t cudaError_t
    #define gpuMalloc cudaMalloc
    #define gpuMallocHost cudaMallocHost
    #define gpuMemcpyHostToDevice cudaMemcpyHostToDevice
    #define gpuMemcpyDeviceToHost cudaMemcpyDeviceToHost
    #define gpuMemcpyToSymbol cudaMemcpyToSymbol
    #define gpuFree cudaFree
    #define gpuFreeHost cudaFreeHost
    #define gpuMemcpy cudaMemcpy
    // Add other necessary abstractions

#endif