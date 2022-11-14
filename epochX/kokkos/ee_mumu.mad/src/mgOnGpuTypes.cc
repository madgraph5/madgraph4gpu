#include "mgOnGpuTypes.h"

fptype fpmax(const fptype& L,const fptype& R){
#ifdef KOKKOS_ENABLE_CUDA
// CUDA DOUBLE
#ifdef MGONGPU_FPTYPE_DOUBLE
  return fmax(L, R);
#else
// CUDA FLOAT
  return fmaxf(L, R);
#endif
// C++
#else
  return std::max(L, R); 
#endif
}

fptype fpmin(const fptype& L,const fptype& R){
#ifdef KOKKOS_ENABLE_CUDA
// CUDA DOUBLE
#ifdef MGONGPU_FPTYPE_DOUBLE
  return fmin(L, R);
#else
// CUDA FLOAT
  return fminf(L, R);
#endif
// C++
#else
  return std::min(L, R); 
#endif
}

fptype fpsqrt(const fptype& V) {
#ifdef KOKKOS_ENABLE_CUDA
// CUDA DOUBLE
#ifdef MGONGPU_FPTYPE_DOUBLE
  return sqrt(V);
#else
// CUDA FLOAT
  return sqrtf(V);
#endif
// C++
#else
  return std::sqrt(V); 
#endif
}

fptype fpabs(const fptype& V) {
#ifdef KOKKOS_ENABLE_CUDA
// CUDA DOUBLE
#ifdef MGONGPU_FPTYPE_DOUBLE
  return fabs(V);
#else
// CUDA FLOAT
  return fabsf(V);
#endif
// C++
#else
  return std::fabs(V); 
#endif
}

int iabs(const int& V) {
#ifdef KOKKOS_ENABLE_CUDA
// CUDA INT
  return abs(V);
// C++
#else
  return std::abs(V); 
#endif
}
