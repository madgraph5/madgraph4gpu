#ifndef BRIDGE_H
#define BRIDGE_H

// includes from Cuda/C++ matrix element calculations
#include "CPPProcess.h"
#include "Memory.h"
#include "mgOnGpuTypes.h"

#include <cassert>
#include <cstring>
#include <iostream>
#include <memory>

// those should become fortran parameters passed in here
const int gpublocks = 1;   // 1024;
const int gputhreads = 16; // 256;

// Forward declare transposition kernel
#ifdef __CUDACC__

template <typename T>
__global__ void dev_transpose(const T *in, T *out, const int evt,
                              const int part, const int mome, const int strd);
#else

template <typename T>
void hst_transpose(const T *in, T *out, const int evt, const int part,
                   const int mome, const int strd);

#endif // __CUDACC__

// *****************************************************************************

/**
 * A templated class for calling the C++ / Cuda matrix element calculations of
 * the event generation workflow. The template parameter is used for the
 * precision of the calculations (float or double)
 *
 * The fortran momenta passed in are in the form of
 *   DOUBLE PRECISION P_MULTI(0:3, NEXTERNAL, NB_PAGE)
 * where the dimensions are <# momenta>, <# of particles>, <# events>
 */
template <typename T> class Bridge {
public:
  /**
   * class constructor
   *
   * @param evt number of events (NB_PAGE, vector.inc)
   * @param par number of particles / event (NEXTERNAL, nexternal.inc)
   * @param mom number of momenta / particle
   * @param str stride length
   * @param ncomb number of good helicities (ncomb, mgOnGpuConfig.h)
   */
  Bridge(int evt, int par, int mom, int str, int ncomb);

  /**
   * sequence to be executed for the Cuda matrix element calculation
   *
   * @param momenta memory address of the input 4-momenta
   * @param mes memory address of the output matrix elements
   */
  void gpu_sequence(T *momenta, double *mes);

  /**
   * sequence to be executed for the vectorized CPU matrix element calculation
   *
   * @param momenta memory address of the input 4-momenta
   * @param mes memory address of the output matrix elements
   */
  void cpu_sequence(T *momenta, double *mes);

private:
  int m_evt;                 ///< number of events
  int m_part;                ///< number of particles / event
  int m_mome;                ///< number of momenta / particle (usually 4)
  int m_strd;                ///< stride length of the AOSOA structure
  int m_ncomb;               ///< number of good helicities
  bool m_goodHelsCalculated; ///< have the good helicities been calculated?

#ifdef __CUDACC__
  typedef std::unique_ptr<bool[], CudaHstDeleter<bool>> CuBHPtr;
  typedef std::unique_ptr<bool, CudaDevDeleter<bool>> CuBDPtr;
  typedef std::unique_ptr<T[], CudaHstDeleter<T>> CuTHPtr;
  typedef std::unique_ptr<T, CudaDevDeleter<T>> CuTDPtr;
  CuTDPtr devMomentaF = devMakeUnique<T>(m_evt * m_part * m_mome);
  CuTDPtr devMomentaC = devMakeUnique<T>(m_evt * m_part * m_mome);
  CuBHPtr hstIsGoodHel = hstMakeUnique<bool>(m_ncomb);
  CuBDPtr devIsGoodHel = devMakeUnique<bool>(m_ncomb);
  CuTDPtr devMEs = devMakeUnique<T>(m_evt);
#else
  // This needs to be inside the function, why?
  // typedef std::unique_ptr<T[], CppHstDeleter<T>> CpTHPtr;
  // typedef std::unique_ptr<bool[], CppHstDeleter<bool>> CpBHPtr;
  // CpTHPtr hstMomenta = hstMakeUnique<T>(m_evt * m_part * m_mome);
  // CpBHPtr hstIsGoodHel2 = hstMakeUnique<bool>(m_ncomb);
  // CpTHPtr hstMEs = hstMakeUnique<T>(m_evt);
#endif
};

// *****************************************************************************

//
// Implementations of class Bridge member functions
//

template <typename T>
Bridge<T>::Bridge(int evnt, int part, int mome, int strd, int ncomb)
    : m_evt(evnt), m_part(part), m_mome(mome), m_strd(strd), m_ncomb(ncomb),
      m_goodHelsCalculated(false) {
#ifdef __CUDACC__
  gProc::CPPProcess process(1, gpublocks, gputhreads, false);
#else
  Proc::CPPProcess process(1, gpublocks, gputhreads, false);
#endif // __CUDACC__
  process.initProc("../../Cards/param_card.dat");
}

#ifdef __CUDACC__

template <typename T> void Bridge<T>::gpu_sequence(T *momenta, double *mes) {

  checkCuda(cudaMemcpy(devMomentaF.get(), momenta,
                       m_evt * m_part * m_mome * sizeof(T),
                       cudaMemcpyHostToDevice));

  dev_transpose<<<gpublocks * 16, gputhreads>>>(
      devMomentaF.get(), devMomentaC.get(), m_evt, m_part, m_mome, m_strd);

  if (!m_goodHelsCalculated) {
    gProc::sigmaKin_getGoodHel<<<gpublocks, gputhreads>>>(
        devMomentaC.get(), devMEs.get(), devIsGoodHel.get());

    checkCuda(cudaMemcpy(hstIsGoodHel.get(), devIsGoodHel.get(),
                         m_ncomb * sizeof(bool), cudaMemcpyDeviceToHost));

    gProc::sigmaKin_setGoodHel(hstIsGoodHel.get());

    m_goodHelsCalculated = true;
  }

  gProc::sigmaKin<<<gpublocks, gputhreads>>>(devMomentaC.get(), devMEs.get());

  checkCuda(
      cudaMemcpy(mes, devMEs.get(), m_evt * sizeof(T), cudaMemcpyDeviceToHost));
}

#else

template <typename T> void Bridge<T>::cpu_sequence(T *momenta, double *mes) {

  // should become class members...
  typedef std::unique_ptr<T[], CppHstDeleter<T>> CpTHPtr;
  typedef std::unique_ptr<bool[], CppHstDeleter<bool>> CpBHPtr;
  CpTHPtr hstMomenta = hstMakeUnique<T>(m_evt * m_part * m_mome);
  CpBHPtr hstIsGoodHel2 = hstMakeUnique<bool>(m_ncomb);
  CpTHPtr hstMEs = hstMakeUnique<T>(m_evt);

  // double(&hstMEs2)[m_evt] = reinterpret_cast<double(&)[m_evt]>(mes);

  hst_transpose(momenta, hstMomenta.get(), m_evt, m_part, m_mome, m_strd);

  if (!m_goodHelsCalculated) {
    Proc::sigmaKin_getGoodHel(hstMomenta.get(), hstMEs.get(),
                              hstIsGoodHel2.get(), m_evt);
    Proc::sigmaKin_setGoodHel(hstIsGoodHel2.get());
    m_goodHelsCalculated = true;
  }

  Proc::sigmaKin(hstMomenta.get(), hstMEs.get(), m_evt);

  memcpy(mes, hstMEs.get(), sizeof(T) * m_evt);
}

#endif // __CUDACC__
// *****************************************************************************

//
// Implementations of transposition functions
//
#ifdef __CUDACC__

/**
const int evnt_n = 4;  // the number of events
const int part_n = 4;  // number of in/out particles inside an event
const int mome_n = 3;  // number of momenta of one particle (usually 4)
const int strd_n = 2;  // stride length for aosoa data (# adjacent events)
const int array_bytes = evnt_n * part_n * mome_n * sizeof(T);
*/
template <typename T>
__global__ void dev_transpose(const T *in, T *out, const int evt,
                              const int part, const int mome, const int strd) {

  int pos = blockDim.x * blockIdx.x + threadIdx.x;
  int arrlen = evt * part * mome;

  if (pos < arrlen) {

    int page_i = pos / (strd * mome * part);
    int rest_1 = pos % (strd * mome * part);
    int part_i = rest_1 / (strd * mome);
    int rest_2 = rest_1 % (strd * mome);
    int mome_i = rest_2 / strd;
    int strd_i = rest_2 % strd;

    int inpos = (page_i * strd + strd_i) // event number
                    * (part * mome)      // event size (pos of event)
                + part_i * mome          // particle inside event
                + mome_i;                // momentum inside particle

    out[pos] = in[inpos];
  }
}

#else

template <typename T>
void hst_transpose(const T *in, T *out, const int evt, const int part,
                   const int mome, const int strd) {

  int arrlen = evt * part * mome;

  for (int pos = 0; pos < arrlen; ++pos) {

    int page_i = pos / (strd * mome * part);
    int rest_1 = pos % (strd * mome * part);
    int part_i = rest_1 / (strd * mome);
    int rest_2 = rest_1 % (strd * mome);
    int mome_i = rest_2 / strd;
    int strd_i = rest_2 % strd;

    int inpos = (page_i * strd + strd_i) // event number
                    * (part * mome)      // event size (pos of event)
                + part_i * mome          // particle inside event
                + mome_i;                // momentum inside particle

    out[pos] = in[inpos];
  }
}

#endif // __CUDACC__

// *****************************************************************************

//
// BACKUP
//

// debug

// std::cout << std::string(80, '*') << std::endl << "Momenta:" << std::endl;
// T *aosoa_p = (T *)hstMomenta.get();
// for (int i = 0; i < m_evt * m_part * m_mome; ++i) {
//   if (i && i % m_strd == 0)
//     std::cout << ", ";
//   if (i && i % (m_mome * m_strd) == 0)
//     std::cout << std::endl;
//   if (i && i % (m_part * m_mome * m_strd) == 0)
//     std::cout << std::endl;
//   std::cout << aosoa_p[i] << " ";
// }
// std::cout << std::endl << std::string(80, '*') << std::endl;

// template <typename T> void Matrix<T>::fill(T *arr) {
//
//   T(*aos)
//   [m_part][m_mome] = (T(*)[m_part][m_mome])arr; // was ->
//   m_hstInpArray.get();
//
//   for (int i = 0; i < m_evt; ++i) {
//     for (int j = 0; j < m_part; ++j) {
//       for (int k = 0; k < m_mome; ++k) {
//         aos[i][j][k] = (i + 1) * 100 + (j + 1) * 10 + (k + 1);
//       }
//     }
//   }
//
// #ifdef DEBUG
//   std::cout << std::string(80, '*') << std::endl;
//   T *aos_p = (T *)arr; // was -> m_hstInpArray.get();
//   for (int i = 0; i < m_evt * m_part * m_mome; ++i) {
//     if (i && i % m_mome == 0)
//       std::cout << std::endl;
//     if (i && i % (m_mome * m_part) == 0)
//       std::cout << std::endl;
//     std::cout << aos_p[i] << " ";
//   }
//   std::cout << std::endl;
// #endif // DEBUG
// }

#endif // BRIDGE_H
