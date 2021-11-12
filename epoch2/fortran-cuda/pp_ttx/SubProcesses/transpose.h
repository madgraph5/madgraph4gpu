#include <cassert>
#include <iostream>
#include <memory>

#include "CPPProcess.h"

/**
const int evnt_n = 4;  // the number of events
const int part_n = 4;  // number of in/out particles inside an event
const int mome_n = 3;  // number of momenta of one particle (usually 4)
const int strd_n = 2;  // stride length for aosoa data (# adjacent events)
const int array_bytes = evnt_n * part_n * mome_n * sizeof(T);
*/

const int gpublocks = 1024;
const int gputhreads = 256;

#ifdef __CUDACC__

#define checkCuda(code)                                                        \
  { assertCuda(code, __FILE__, __LINE__); }

/*
inline void assertCuda(cudaError_t code, const char *file, int line,
                       bool abort = true) {
  if (code != cudaSuccess) {
    printf("GPUassert: %s %s:%d\n", cudaGetErrorString(code), file, line);
    if (abort)
      assert(code == cudaSuccess);
  }
}*/

template <typename T> struct CudaDevDeleter {
  void operator()(T *mem) { checkCuda(cudaFree(mem)); }
};

template <typename T>
std::unique_ptr<T[], CudaDevDeleter<T>> devMakeUnique(int n_bytes) {
  T *tmp = nullptr;
  checkCuda(cudaMalloc(&tmp, n_bytes));
  return std::unique_ptr<T[], CudaDevDeleter<T>>{tmp};
}

template <typename T> struct CudaHstDeleter {
  void operator()(T *mem) { checkCuda(cudaFreeHost(mem)); }
};

template <typename T>
std::unique_ptr<T[], CudaHstDeleter<T>> hstMakeUnique(int n_bytes) {
  T *tmp = nullptr;
  checkCuda(cudaMallocHost(&tmp, n_bytes));
  return std::unique_ptr<T[], CudaHstDeleter<T>>{tmp};
};

template <typename T>
__global__ void dev_transpose(const T *inpArr, T *outArr, const int evnt_n,
                              const int part_n, const int mome_n,
                              const int strd_n) {

  int pos = blockDim.x * blockIdx.x + threadIdx.x;
  int arrlen = evnt_n * part_n * mome_n;

  if (pos < arrlen) {

    int page_i = pos / (strd_n * mome_n * part_n);
    int rest_1 = pos % (strd_n * mome_n * part_n);
    int part_i = rest_1 / (strd_n * mome_n);
    int rest_2 = rest_1 % (strd_n * mome_n);
    int mome_i = rest_2 / strd_n;
    int strd_i = rest_2 % strd_n;

    int inpos = (page_i * strd_n + strd_i) // event number
                    * (part_n * mome_n)    // event size (pos of event)
                + part_i * mome_n          // particle inside event
                + mome_i;                  // momentum inside particle

#ifdef DEBUG
    printf("opos:%d, ipos:%d, page_i:%d, strd_i:%d, part_i:%i, mome_i:%d\n",
           pos, inpos, page_i, strd_i, part_i, mome_i);
#endif

    outArr[pos] = inpArr[inpos];
  }
}

#endif // __CUDACC__

template <typename T> class Matrix {
public:
  Matrix(int evt, int par, int mom, int str);

  void fill(T *arr);
  void hst_transpose(T *arr);

private:
  int m_evnt;
  int m_part;
  int m_mome;
  int m_strd;
  int m_arrbytes;
};

/**
 *
 */
template <typename T>
Matrix<T>::Matrix(int evnt, int part, int mome, int strd)
    : m_evnt(evnt), m_part(part), m_mome(mome), m_strd(strd),
      m_arrbytes(m_evnt * m_part * m_mome * sizeof(T)) {}

/**
 *
 */
template <typename T> void Matrix<T>::fill(T *arr) {

  T(*aos)
  [m_part][m_mome] = (T(*)[m_part][m_mome])arr; // was -> m_hstInpArray.get();

  for (int i = 0; i < m_evnt; ++i) {
    for (int j = 0; j < m_part; ++j) {
      for (int k = 0; k < m_mome; ++k) {
        aos[i][j][k] = (i + 1) * 100 + (j + 1) * 10 + (k + 1);
      }
    }
  }

#ifdef DEBUG
  std::cout << std::string(80, '*') << std::endl;
  T *aos_p = (T *)arr; // was -> m_hstInpArray.get();
  for (int i = 0; i < m_evnt * m_part * m_mome; ++i) {
    if (i && i % m_mome == 0)
      std::cout << std::endl;
    if (i && i % (m_mome * m_part) == 0)
      std::cout << std::endl;
    std::cout << aos_p[i] << " ";
  }
  std::cout << std::endl;
#endif // DEBUG
}

/**
 *
 */
template <typename T> void Matrix<T>::hst_transpose(T *arr) {

  std::unique_ptr<T[], CudaHstDeleter<T>> hstInpArray =
      hstMakeUnique<T>(m_arrbytes);
  std::unique_ptr<T[], CudaDevDeleter<T>> devInpArray =
      devMakeUnique<T>(m_arrbytes);
  std::unique_ptr<T[], CudaHstDeleter<T>> hstOutArray =
      hstMakeUnique<T>(m_arrbytes);
  std::unique_ptr<T[], CudaDevDeleter<T>> devOutArray =
      devMakeUnique<T>(m_arrbytes);

  checkCuda(cudaMemcpy(devInpArray.get(), hstInpArray.get(), m_arrbytes,
                       cudaMemcpyHostToDevice));

  dev_transpose<<<gpublocks, gputhreads>>>(devInpArray.get(), devOutArray.get(),
                                           m_evnt, m_part, m_mome, m_strd);

  gProc::sigmaKin<<<gpublocks, gputhreads>>>(devInpArray.get(),
                                             devOutArray.get());

  checkCuda(cudaMemcpy(hstOutArray.get(), devOutArray.get(), m_arrbytes,
                       cudaMemcpyDeviceToHost));

#ifdef DEBUG
  std::cout << std::string(80, '*') << std::endl;
  T *aosoa_p = (T *)hstOutArray.get();
  for (int i = 0; i < evnt_n * part_n * mome_n; ++i) {
    if (i && i % m_strd == 0)
      std::cout << ", ";
    if (i && i % (m_mome * m_strd) == 0)
      std::cout << std::endl;
    if (i && i % (m_part * m_mome * m_strd) == 0)
      std::cout << std::endl;
    std::cout << aosoa_p[i] << " ";
  }
  std::cout << std::endl << std::string(80, '*') << std::endl;
#endif // DEBUG
}
