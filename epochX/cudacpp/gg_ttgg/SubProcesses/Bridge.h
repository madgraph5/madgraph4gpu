#ifndef BRIDGE_H
#define BRIDGE_H 1

// Includes from Cuda/C++ matrix element calculations
#include "mgOnGpuConfig.h" // for mgOnGpu::npar, mgOnGpu::np4
#include "CPPProcess.h" // for CPPProcess
#include "MatrixElementKernels.h" // for MatrixElementKernelHost, MatrixElementKernelDevice
#include "MemoryAccessMomenta.h" // for MemoryAccessMomenta::neppM
#include "MemoryBuffers.h" // for HostBufferMomenta, DeviceBufferMomenta etc

#include <cassert>
#include <cmath>
#include <cstring>
#include <iostream>
#include <memory>

// Forward declare transposition kernels
// (Inverse transposition with template bool parameter F2C and "if constexpr" would require C++17, not available in cuda 11.1)

#ifdef __CUDACC__

template <typename T>
__global__
void dev_transposeMomentaF2C( const T *in, T *out, const int evt );

#endif // __CUDACC__

template <typename T>
void hst_transposeMomentaF2C( const T *in, T *out, const int evt );

template <typename T>
void hst_transposeMomentaC2F( const T *in, T *out, const int evt );

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
template <typename T> class Bridge
{
public:
  /**
   * class constructor
   *
   * @param evt number of events (NB_PAGE, vector.inc)                    [STRICTLY NEEDED]
   * @param par number of particles per event (NEXTERNAL, nexternal.inc)  [KEPT FOR SANITY CHECKS ONLY! remove it?]
   * @param mom number of momenta per particle                            [KEPT FOR SANITY CHECKS ONLY! remove it?]
   * @param str stride length                                             [IGNORED! REMOVE IT! fully encapsulated in Cuda/C++]
   * @param ncomb number of good helicities                               [IGNORED! REMOVE IT! fully encapsulated in Cuda/C++]
   */
  Bridge( int evt, int par, int mom, int /*dummy*/str=0, int /*dummy*/ncomb=0 );

#ifdef __CUDACC__
  /**
   * set the gpublocks and gputhreads for the gpusequence - throws if evnt != gpublocks*gputhreads
   * (this is needed for BridgeKernel tests rather than for actual production use in Fortran)
   *
   * @param gpublocks number of gpublocks
   * @param gputhreads number of gputhreads
   */
  void set_gpugrid( const int gpublocks, const int gputhreads );
#endif

  /**
   * sequence to be executed for the Cuda matrix element calculation
   *
   * @param momenta memory address of the input 4-momenta
   * @param mes memory address of the output matrix elements
   * @param goodHelOnly quit after computing good helicities
   */
  void gpu_sequence( const T *momenta, double *mes, const bool goodHelOnly=false );

  /**
   * sequence to be executed for the vectorized CPU matrix element calculation
   *
   * @param momenta memory address of the input 4-momenta
   * @param mes memory address of the output matrix elements
   * @param goodHelOnly quit after computing good helicities
   */
  void cpu_sequence( const T *momenta, double *mes, const bool goodHelOnly=false );

private:

  int m_evt;                 ///< number of events
  //int m_part;                ///< number of particles per event [NO LONGER NEEDED!]
  //int m_mome;                ///< number of momenta per particle (usually 4) [NO LONGER NEEDED!]
  //int m_strd;                ///< stride length of the AOSOA structure [NO LONGER NEEDED!]
  //int m_ncomb;               ///< number of good helicities [NO LONGER NEEDED!]
  bool m_goodHelsCalculated; ///< have the good helicities been calculated?

  int m_gputhreads;          ///< number of gpu threads (default set from number of events, can be modified)
  int m_gpublocks;           ///< number of gpu blocks (default set from number of events, can be modified)

#ifdef __CUDACC__
  mg5amcGpu::DeviceBufferMomenta m_devMomentaF;
  mg5amcGpu::DeviceBufferMomenta m_devMomentaC;
  mg5amcGpu::DeviceBufferMatrixElements m_devMEsC;
  mg5amcGpu::MatrixElementKernelDevice m_devMek;
#else
  mg5amcCpu::HostBufferMomenta m_hstMomentaC;
  mg5amcCpu::HostBufferMatrixElements m_hstMEsC;
  mg5amcCpu::MatrixElementKernelHost m_hstMek;
#endif

};

// *****************************************************************************

//
// Implementations of class Bridge member functions
//

template <typename T>
Bridge<T>::Bridge( int evnt, int part, int mome, int /*strd*/, int /*ncomb*/ )
  : m_evt( evnt )
    //, m_part( part )
    //, m_mome( mome )
    //, m_strd( strd )
    //, m_ncomb( ncomb )
  , m_goodHelsCalculated( false )
  , m_gputhreads( 256 ) // default number of gpu threads
  , m_gpublocks( ceil(double(m_evt)/m_gputhreads) ) // this ensures m_evt <= m_gpublocks*m_gputhreads
#ifdef __CUDACC__
  , m_devMomentaF( evnt )
  , m_devMomentaC( evnt )
  , m_devMEsC( evnt )
  , m_devMek( m_devMomentaC, m_devMEsC, m_gpublocks, m_gputhreads )
#else
  , m_hstMomentaC( evnt )
  , m_hstMEsC( evnt )
  , m_hstMek( m_hstMomentaC, m_hstMEsC, evnt )
#endif
{
  if ( part != mgOnGpu::npar ) throw std::runtime_error( "Bridge constructor: npar mismatch" );
  if ( mome != mgOnGpu::np4 ) throw std::runtime_error( "Bridge constructor: np4 mismatch" );
  std::cout << "WARNING! Instantiate Bridge (nevt=" << m_evt << ", gpublocks=" << m_gpublocks << ", gputhreads=" << m_gputhreads
            << ", gpublocks*gputhreads=" << m_gpublocks*m_gputhreads << ")" << std::endl;
#ifdef __CUDACC__
  mg5amcGpu::CPPProcess process( 1, m_gpublocks, m_gputhreads, false );
#else
  mg5amcCpu::CPPProcess process( 1, m_gpublocks, m_gputhreads, false );
#endif // __CUDACC__
  process.initProc( "../../Cards/param_card.dat" );
}

#ifdef __CUDACC__
template <typename T>
void Bridge<T>::set_gpugrid(const int gpublocks, const int gputhreads)
{
  if ( m_goodHelsCalculated )
    throw std::runtime_error( "Bridge: gpublocks and gputhreads cannot be set after calculating helicities" );
  if ( m_evt != gpublocks*gputhreads )
    throw std::runtime_error( "Bridge: gpublocks*gputhreads must equal m_evt in set_gpugrid" );
  m_gpublocks = gpublocks;
  m_gputhreads = gputhreads;
  std::cout << "WARNING! Set grid in Bridge (nevt=" << m_evt << ", gpublocks=" << m_gpublocks << ", gputhreads=" << m_gputhreads
            << ", gpublocks*gputhreads=" << m_gpublocks*m_gputhreads << ")" << std::endl;
  m_devMek.setGrid( m_gpublocks, m_gputhreads );
}
#endif

#ifdef __CUDACC__
template <typename T>
void Bridge<T>::gpu_sequence( const T *momenta, double *mes, const bool goodHelOnly )
{
  constexpr int neppM = MemoryAccessMomenta::neppM;
  if ( neppM == 1 ) // eventually move to "if constexpr" (need c++17, not available in cuda 11.1)
  {
    checkCuda( cudaMemcpy( m_devMomentaC.data(), momenta, m_devMomentaC.bytes(), cudaMemcpyHostToDevice ) );
  }
  else
  {
    checkCuda( cudaMemcpy( m_devMomentaF.data(), momenta, m_devMomentaF.bytes(), cudaMemcpyHostToDevice ) );
    const int thrPerEvt = mgOnGpu::npar * mgOnGpu::np4; // AV: transpose alg does 1 element per thread (NOT 1 event per thread)
    //const int thrPerEvt = 1; // AV: try new alg with 1 event per thread... this seems slower
    dev_transposeMomentaF2C<<<m_gpublocks*thrPerEvt, m_gputhreads>>>( m_devMomentaF.data(), m_devMomentaC.data(), m_evt );
  }
  if ( !m_goodHelsCalculated )
  {
    m_devMek.computeGoodHelicities();
    m_goodHelsCalculated = true;
  }
  if ( goodHelOnly ) return;
  m_devMek.computeMatrixElements();
  checkCuda( cudaMemcpy( mes, m_devMEsC.data(), m_devMEsC.bytes(), cudaMemcpyDeviceToHost ) );
}
#endif

#ifndef __CUDACC__
template <typename T>
void Bridge<T>::cpu_sequence( const T *momenta, double *mes, const bool goodHelOnly )
{
  hst_transposeMomentaF2C( momenta, m_hstMomentaC.data(), m_evt );
  if ( !m_goodHelsCalculated )
  {
    m_hstMek.computeGoodHelicities();
    m_goodHelsCalculated = true;
  }
  if ( goodHelOnly ) return;
  m_hstMek.computeMatrixElements();
  memcpy( mes, m_hstMEsC.data(), m_hstMEsC.bytes() );
}
#endif

// *****************************************************************************

//
// Implementations of transposition functions
//

#ifdef __CUDACC__
template <typename T>
__global__
void dev_transposeMomentaF2C( const T *in, T *out, const int evt )
{
  constexpr int part = mgOnGpu::npar;
  constexpr int mome = mgOnGpu::np4;
  constexpr int strd = MemoryAccessMomenta::neppM;
  int pos = blockDim.x * blockIdx.x + threadIdx.x;
  int arrlen = evt * part * mome;
  if (pos < arrlen)
  {
    int page_i = pos / (strd * mome * part);
    int rest_1 = pos % (strd * mome * part);
    int part_i = rest_1 / (strd * mome);
    int rest_2 = rest_1 % (strd * mome);
    int mome_i = rest_2 / strd;
    int strd_i = rest_2 % strd;
    int inpos =
      (page_i * strd + strd_i) // event number
      * (part * mome)          // event size (pos of event)
      + part_i * mome          // particle inside event
      + mome_i;                // momentum inside particle
    out[pos] = in[inpos]; // F2C (Fortran to C)
  }
  /*
  // AV attempt another implementation with 1 event per thread: this seems slower...
  // C-style: AOSOA[npagM][npar][np4][neppM]
  // F-style: AOS[nevt][npar][np4]
  constexpr int npar = mgOnGpu::npar;
  constexpr int np4 = mgOnGpu::np4;
  constexpr int neppM = MemoryAccessMomenta::neppM;
  const int npagM = evt/neppM;
  assert( evt%neppM == 0 ); // number of events is not a multiple of neppM???
  int ievt = blockDim.x * blockIdx.x + threadIdx.x;
  int ipagM = ievt/neppM;
  int ieppM = ievt%neppM;
  for ( int ip4=0; ip4<np4; ip4++ )
    for ( int ipar=0; ipar<npar; ipar++ )
    {
      int cpos = ipagM*npar*np4*neppM + ipar*np4*neppM + ip4*neppM + ieppM;
      int fpos = ievt*npar*np4 + ipar*np4 + ip4;
      out[cpos] = in[fpos]; // F2C (Fortran to C)
    }
  */
}
#endif

template <typename T>
void hst_transposeMomentaF2C( const T *in, T *out, const int evt )
{
  /*
  constexpr int part = mgOnGpu::npar;
  constexpr int mome = mgOnGpu::np4;
  constexpr int strd = MemoryAccessMomenta::neppM;
  int arrlen = evt * part * mome;
  for (int pos = 0; pos < arrlen; ++pos)
  {
    int page_i = pos / (strd * mome * part);
    int rest_1 = pos % (strd * mome * part);
    int part_i = rest_1 / (strd * mome);
    int rest_2 = rest_1 % (strd * mome);
    int mome_i = rest_2 / strd;
    int strd_i = rest_2 % strd;
    int inpos =
      (page_i * strd + strd_i) // event number
      * (part * mome)          // event size (pos of event)
      + part_i * mome          // particle inside event
      + mome_i;                // momentum inside particle
    out[pos] = in[inpos]; // F2C (Fortran to C)
  }
  */
  // AV attempt another implementation: this is slightly faster (better c++ pipelining?)
  // [NB! this is not a transposition, it is an AOS to AOSOA conversion: if neppM=1, a memcpy is enough]
  // C-style: AOSOA[npagM][npar][np4][neppM]
  // F-style: AOS[nevt][npar][np4]
  constexpr int npar = mgOnGpu::npar;
  constexpr int np4 = mgOnGpu::np4;
  constexpr int neppM = MemoryAccessMomenta::neppM;
  if ( neppM == 1 ) // eventually move to "if constexpr" (need c++17, not available in cuda 11.1)
  {
    memcpy( out, in, evt * npar * np4 * sizeof(T) );
  }
  else
  {
    const int npagM = evt/neppM;
    assert( evt%neppM == 0 ); // number of events is not a multiple of neppM???
    for ( int ipagM=0; ipagM<npagM; ipagM++ )
      for ( int ip4=0; ip4<np4; ip4++ )
        for ( int ipar=0; ipar<npar; ipar++ )
          for ( int ieppM=0; ieppM<neppM; ieppM++ )
          {
            int ievt = ipagM*neppM + ieppM;
            int cpos = ipagM*npar*np4*neppM + ipar*np4*neppM + ip4*neppM + ieppM;
            int fpos = ievt*npar*np4 + ipar*np4 + ip4;
            out[cpos] = in[fpos]; // F2C (Fortran to C)
          }
  }
}

template <typename T>
void hst_transposeMomentaC2F( const T *in, T *out, const int evt )
{
  constexpr int part = mgOnGpu::npar;
  constexpr int mome = mgOnGpu::np4;
  constexpr int strd = MemoryAccessMomenta::neppM;
  int arrlen = evt * part * mome;
  for (int pos = 0; pos < arrlen; ++pos)
  {
    int page_i = pos / (strd * mome * part);
    int rest_1 = pos % (strd * mome * part);
    int part_i = rest_1 / (strd * mome);
    int rest_2 = rest_1 % (strd * mome);
    int mome_i = rest_2 / strd;
    int strd_i = rest_2 % strd;
    int inpos =
      (page_i * strd + strd_i) // event number
      * (part * mome)          // event size (pos of event)
      + part_i * mome          // particle inside event
      + mome_i;                // momentum inside particle
    out[inpos] = in[pos]; // C2F (C to Fortran)
  }
}

// *****************************************************************************

//
// BACKUP
//

/**
   const int evnt_n = 4;  // the number of events
   const int part_n = 4;  // number of in/out particles inside an event
   const int mome_n = 3;  // number of momenta of one particle (usually 4)
   const int strd_n = 2;  // stride length for aosoa data (# adjacent events)
   const int array_bytes = evnt_n * part_n * mome_n * sizeof(T);
*/

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
