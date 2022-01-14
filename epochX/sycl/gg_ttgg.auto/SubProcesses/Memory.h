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

template<typename T = fptype>
struct CppHstDeleter {
  void operator()(T* mem) {
    ::operator delete( mem, std::align_val_t{ mgOnGpu::cppAlign } );
    //::operator delete( mem-1, std::align_val_t{ mgOnGpu::cppAlign } ); // TEST MISALIGNMENT!
  }
};

template<typename T = fptype> inline
std::unique_ptr<T[], CppHstDeleter<T>> hstMakeUnique(std::size_t N) {
  // See https://www.bfilipek.com/2019/08/newnew-align.html#requesting-different-alignment
  return std::unique_ptr<T[], CppHstDeleter<T>>{ new( std::align_val_t{ mgOnGpu::cppAlign } ) T[N]() };
  //return std::unique_ptr<T[], CppHstDeleter<T>>{ new( std::align_val_t{ mgOnGpu::cppAlign } ) T[N+1]() + 1 }; // TEST MISALIGNMENT!
};

#endif /* MEMORY_H */
