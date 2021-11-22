#include <cassert>
#include <iostream>
#include <memory>

#include "CPPProcess.h"
#include "Memory.h"

/**
const int evnt_n = 4;  // the number of events
const int part_n = 4;  // number of in/out particles inside an event
const int mome_n = 3;  // number of momenta of one particle (usually 4)
const int strd_n = 2;  // stride length for aosoa data (# adjacent events)
const int array_bytes = evnt_n * part_n * mome_n * sizeof(T);
*/

const int gpublocks = 1;   // 1024;
const int gputhreads = 16; // 256;

#ifdef __CUDACC__

#define checkCuda(code)                                                        \
  { assertCuda(code, __FILE__, __LINE__); }

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

    outArr[pos] = inpArr[inpos];
    //     if (pos == 0 and part_i == 0 and mome_i == 0)
    //       printf("tran: %d %d %d %d %f\n", blockIdx.x, threadIdx.x, part_i,
    //       mome_i,
    //              inpArr[inpos]);
    // #ifdef DEBUG

    // printf("opos:%d, ipos:%d, evt_i:%d, part_i:%i, mome_i:%d, page_i:%d, "
    //        "strd_i:%d, val:%f\n",
    //        pos, inpos, (page_i * strd_n + strd_i), part_i, mome_i, page_i,
    //        strd_i, (T)outArr[pos]);

    // #endif
  }
}

#endif // __CUDACC__

template <typename T> class Matrix {
public:
  Matrix(int evt, int par, int mom, int str, int ncomb);

  void fill(T *arr);
  void hst_transpose(T *momenta, double *mes);

private:
  int m_evnt;
  int m_part;
  int m_mome;
  int m_strd;
  int m_ncomb;
  int m_bts_momenta;
  int m_bts_goodhel;
  int m_bts_mes;
};

/**
 *
 */
template <typename T>
Matrix<T>::Matrix(int evnt, int part, int mome, int strd, int ncomb)
    : m_evnt(evnt), m_part(part), m_mome(mome), m_strd(strd), m_ncomb(ncomb),
      m_bts_momenta(m_evnt * m_part * m_mome * sizeof(T)),
      m_bts_goodhel(m_ncomb * sizeof(bool)), m_bts_mes(m_evnt * sizeof(T)) {
  gProc::CPPProcess process(1, gpublocks, gputhreads, false);
  process.initProc("../../Cards/param_card.dat");
}

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
template <typename T> void Matrix<T>::hst_transpose(T *momenta, double *mes) {

  auto devMomentaF2 = devMakeUnique<T>(m_evnt * m_part * m_mome);
  auto devMomentaC2 = devMakeUnique<T>(m_evnt * m_part * m_mome);
  auto hstIsGoodHel2 = hstMakeUnique<bool>(m_ncomb);
  auto devIsGoodHel2 = devMakeUnique<bool>(m_ncomb);
  auto hstMEs2 = hstMakeUnique<T>(m_evnt);
  auto devMEs2 = devMakeUnique<T>(m_evnt);

  checkCuda(cudaMemcpy(devMomentaF2.get(), momenta, m_bts_momenta,
                       cudaMemcpyHostToDevice));

  dev_transpose<<<gpublocks * 16, gputhreads>>>(
      devMomentaF2.get(), devMomentaC2.get(), m_evnt, m_part, m_mome, m_strd);

  // sr should go before //
  gProc::sigmaKin_getGoodHel<<<gpublocks, gputhreads>>>(
      devMomentaC2.get(), devMEs2.get(), devIsGoodHel2.get());

  checkCuda(cudaMemcpy(hstIsGoodHel2.get(), devIsGoodHel2.get(), m_bts_goodhel,
                       cudaMemcpyDeviceToHost));

  gProc::sigmaKin_setGoodHel(hstIsGoodHel2.get());

  gProc::sigmaKin<<<gpublocks, gputhreads>>>(devMomentaC2.get(), devMEs2.get());

  // sr copy directly into me
  checkCuda(cudaMemcpy(hstMEs2.get(), devMEs2.get(), m_bts_mes,
                       cudaMemcpyDeviceToHost));

  auto phstMEs = hstMEs2.get();
  size_t s = 16 * sizeof(double);
  memcpy(mes, hstMEs2.get(), s);
  // mes = phstMEs;

  // std::cout << "MEs: ";
  // for (int i = 0; i < 16; ++i) {
  //   std::cout << mes[i] << ", ";
  // }
  // std::cout << std::endl << std::endl;

  // std::cout << std::string(80, '*') << std::endl;
  // for (int i = 0; i < m_ncomb; ++i) {
  //   std::cout << i << ":" << hstIsGoodHel2[i] << " ";
  // }
  // std::cout << std::endl;

  // std::cout << std::string(80, '*') << std::endl;
  // for (int i = 0; i < m_evnt; ++i) {
  //   std::cout << i << ":" << hstMEs2[i] << " ";
  // }
  // std::cout << std::endl;

#ifdef DEBUG
  auto hstMomentaC2 = hstMakeUnique<T>(m_evnt * m_part * m_mome);
  checkCuda(cudaMemcpy(hstMomentaC2.get(), devMomentaC2.get(), m_bts_momenta,
                       cudaMemcpyDeviceToHost));

  std::cout << std::string(80, '*') << std::endl;
  T *aosoa_p = (T *)hstMomentaC2.get();
  for (int i = 0; i < m_evnt * m_part * m_mome; ++i) {
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
