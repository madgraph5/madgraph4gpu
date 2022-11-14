/*
 *  Memory.h
 *
 *  Created on: 19.11.2020
 *      Author: shageboeck
 *    Modified: avalassi (vectors, alignment)
 */

#ifndef MEMORY_H
#define MEMORY_H 1

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"
#include "mgOnGpuVectors.h"

#include <memory>

#ifdef __CUDACC__

template<typename T = fptype>
struct CudaDevDeleter {
  void operator()(T* mem) {
    checkCuda( cudaFree( mem ) );
  }
};

template<typename T = fptype>
std::unique_ptr<T, CudaDevDeleter<T>> devMakeUnique(std::size_t N) {
  T* tmp = nullptr;
  checkCuda( cudaMalloc( &tmp, N * sizeof(T) ) );
  return std::unique_ptr<T, CudaDevDeleter<T>>{ tmp };
}

template<typename T = fptype>
struct CudaHstDeleter {
  void operator()(T* mem) {
    checkCuda( cudaFreeHost( mem ) );
  }
};

template<typename T = fptype>
std::unique_ptr<T[], CudaHstDeleter<T>> hstMakeUnique(std::size_t N) {
  T* tmp = nullptr;
  checkCuda( cudaMallocHost( &tmp, N * sizeof(T) ) );
  return std::unique_ptr<T[], CudaHstDeleter<T>>{ tmp };
};

#else

template<typename T = fptype>
struct CppHstDeleter {
  void operator()(T* mem) {
    ::operator delete( mem, std::align_val_t{ mgOnGpu::cppAlign } );
    //::operator delete( mem-1, std::align_val_t{ mgOnGpu::cppAlign } ); // TEST MISALIGNED ACCESS!
  }
};

template<typename T = fptype> inline
std::unique_ptr<T[], CppHstDeleter<T>> hstMakeUnique(std::size_t N) {
  // See https://www.bfilipek.com/2019/08/newnew-align.html#requesting-different-alignment
  return std::unique_ptr<T[], CppHstDeleter<T>>{ new( std::align_val_t{ mgOnGpu::cppAlign } ) T[N]() };
  //return std::unique_ptr<T[], CppHstDeleter<T>>{ new( std::align_val_t{ mgOnGpu::cppAlign } ) T[N+1]() + 1 }; // TEST MISALIGNED ACCESS!
};

/*
#ifdef MGONGPU_CPPSIMD

template<> inline
std::unique_ptr<fptype_v[], CppHstDeleter<fptype_v>> hstMakeUnique(std::size_t N) {
  return std::unique_ptr<fptype_v[], CppHstDeleter<fptype_v>>{ new( std::align_val_t{ mgOnGpu::cppAlign } ) fptype_v[N/neppV]() };
};

#endif
*/

#endif

#endif /* MEMORY_H */
