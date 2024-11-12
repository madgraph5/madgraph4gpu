// Copyright (C) 2020-2024 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: J. Teig (Jul 2023) for the MG5aMC CUDACPP plugin.
// Further modified by: J. Teig, A. Valassi (2020-2024) for the MG5aMC CUDACPP plugin.

#ifndef MG5AMC_GPUABSTRACTION_H
#define MG5AMC_GPUABSTRACTION_H 1

#include <cassert>

//--------------------------------------------------------------------------

#ifdef __CUDACC__ // this must be __CUDACC__ (not MGONGPUCPP_GPUIMPL)

#define gpuError_t cudaError_t
#define gpuPeekAtLastError cudaPeekAtLastError
#define gpuGetErrorString cudaGetErrorString
#define gpuSuccess cudaSuccess

#define gpuMallocHost( ptr, size ) checkGpu( cudaMallocHost( ptr, size ) )
#define gpuMalloc( ptr, size ) checkGpu( cudaMalloc( ptr, size ) )

#define gpuMemcpy( dstData, srcData, srcBytes, func ) checkGpu( cudaMemcpy( dstData, srcData, srcBytes, func ) )
#define gpuMemset( data, value, bytes ) checkGpu( cudaMemset( data, value, bytes ) )
#define gpuMemcpyHostToDevice cudaMemcpyHostToDevice
#define gpuMemcpyDeviceToHost cudaMemcpyDeviceToHost
#define gpuMemcpyDeviceToDevice cudaMemcpyDeviceToDevice
#define gpuMemcpyToSymbol( type1, type2, size ) checkGpu( cudaMemcpyToSymbol( type1, type2, size ) )

#define gpuFree( ptr ) checkGpu( cudaFree( ptr ) )
#define gpuFreeHost( ptr ) checkGpu( cudaFreeHost( ptr ) )

#define gpuGetSymbolAddress( devPtr, symbol ) checkGpu( cudaGetSymbolAddress( devPtr, symbol ) )

#define gpuSetDevice cudaSetDevice
#define gpuDeviceSynchronize cudaDeviceSynchronize
#define gpuDeviceReset cudaDeviceReset

#define gpuLaunchKernel( kernel, blocks, threads, ... ) kernel<<<blocks, threads>>>( __VA_ARGS__ )
//#define gpuLaunchKernelSharedMem( kernel, blocks, threads, sharedMem, ... ) kernel<<<blocks, threads, sharedMem>>>( __VA_>
#define gpuLaunchKernelStream( kernel, blocks, threads, stream, ... ) kernel<<<blocks, threads, 0, stream>>>( __VA_ARGS__ )

#define gpuStream_t cudaStream_t
#define gpuStreamCreate cudaStreamCreate
#define gpuStreamDestroy cudaStreamDestroy

//--------------------------------------------------------------------------

#elif defined __HIPCC__

#define gpuError_t hipError_t
#define gpuPeekAtLastError hipPeekAtLastError
#define gpuGetErrorString hipGetErrorString
#define gpuSuccess hipSuccess

#define gpuMallocHost( ptr, size ) checkGpu( hipHostMalloc( ptr, size ) ) // HostMalloc better
#define gpuMalloc( ptr, size ) checkGpu( hipMalloc( ptr, size ) )

#define gpuMemcpy( dstData, srcData, srcBytes, func ) checkGpu( hipMemcpy( dstData, srcData, srcBytes, func ) )
#define gpuMemset( data, value, bytes ) checkGpu( hipMemset( data, value, bytes ) )
#define gpuMemcpyHostToDevice hipMemcpyHostToDevice
#define gpuMemcpyDeviceToHost hipMemcpyDeviceToHost
#define gpuMemcpyDeviceToDevice hipMemcpyDeviceToDevice
#define gpuMemcpyToSymbol( type1, type2, size ) checkGpu( hipMemcpyToSymbol( type1, type2, size ) )

#define gpuFree( ptr ) checkGpu( hipFree( ptr ) )
#define gpuFreeHost( ptr ) checkGpu( hipHostFree( ptr ) )

#define gpuGetSymbolAddress( devPtr, symbol ) checkGpu( hipGetSymbolAddress( devPtr, symbol ) )

#define gpuSetDevice hipSetDevice
#define gpuDeviceSynchronize hipDeviceSynchronize
#define gpuDeviceReset hipDeviceReset

#define gpuLaunchKernel( kernel, blocks, threads, ... ) kernel<<<blocks, threads>>>( __VA_ARGS__ )
//#define gpuLaunchKernelSharedMem( kernel, blocks, threads, sharedMem, ... ) kernel<<<blocks, threads, sharedMem>>>( __VA_>
#define gpuLaunchKernelStream( kernel, blocks, threads, stream, ... ) kernel<<<blocks, threads, 0, stream>>>( __VA_ARGS__ )

#define gpuStream_t hipStream_t
#define gpuStreamCreate hipStreamCreate
#define gpuStreamDestroy hipStreamDestroy

//--------------------------------------------------------------------------

#endif

#endif // MG5AMC_GPUABSTRACTION_H
