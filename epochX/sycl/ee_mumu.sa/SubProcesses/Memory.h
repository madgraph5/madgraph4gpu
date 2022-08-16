/*
 * memory.h
 *
 *  Created on: 19.11.2020
 *      Author: shageboeck
 */

#ifndef MEMORY_H
#define MEMORY_H 1

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"
#include "mgOnGpuVectors.h"

#include <memory>

template<typename T = fptype> inline
std::unique_ptr<T[]> hstMakeUnique(std::size_t N) { return std::unique_ptr<T[]>{ new T[N]() }; };

#endif /* MEMORY_H */
