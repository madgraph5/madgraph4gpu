// gpu_abstraction.h
#pragma once

#ifdef __CUDACC__
    // NVIDIA GPU using CUDA
    #include <CudaRuntime.h>

    #define gpuError_t cudaError_t
    #define gpuPeekAtLastError cudaPeekAtLastError

    #define gpuMallocHost(ptr, size) checkCuda( cudaMallocHost(ptr, size) )
    #define gpuMalloc(ptr, size) checkCuda( cudaMalloc(ptr, size) )

    #define gpuMemcpy(dstData, srcData, srcBytes, func) checkCuda( cudaMemcpy(dstData, srcData, srcBytes, func) )
    #define gpuMemcpyHostToDevice cudaMemcpyHostToDevice
    #define gpuMemcpyDeviceToHost cudaMemcpyDeviceToHost
    #define gpuMemcpyToSymbol(type1, type2, size) checkCuda( cudaMemcpyToSymbol(type1, type2, size) )

    #define gpuFree(ptr) checkCuda( cudaFree(ptr) )
    #define gpuFreeHost(ptr) checkCuda( cudaFreeHost(ptr) )

    #define gpuDeviceSynchronize cudaDeviceSynchronize
    #define gpuDeviceReset cudaDeviceReset

#elif defined(__HIP__)
    // AMD GPU using HIP
    #include <hip/hip_runtime.h>

    #define gpuError_t hipError_t
    #define gpuPeekAtLastError hipPeekAtLastError

    #define gpuMallocHost(ptr, size) checkHip( hipMallocHost(ptr, size) )
    #define gpuMalloc(ptr, size) checkHip( hipMalloc(ptr, size) )

    #define gpuMemcpy(dstData, srcData, srcBytes, func) checkHip( hipMemcpy(dstData, srcData, srcBytes, func) )
    #define gpuMemcpyHostToDevice hipMemcpyHostToDevice
    #define gpuMemcpyDeviceToHost hipMemcpyDeviceToHost
    #define gpuMemcpyToSymbol(type1, type2, size) checkHip( hipMemcpyToSymbol(type1, type2, size) )

    #define gpuFree(ptr) checkHip( hipFree(ptr) )
    #define gpuFreeHost(ptr) checkHip( hipFreeHost(ptr) )

    #define gpuDeviceSynchronize hipDeviceSynchronize
    #define gpuDeviceReset hipDeviceReset

#endif