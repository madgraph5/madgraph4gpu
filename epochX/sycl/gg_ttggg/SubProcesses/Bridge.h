#ifndef BRIDGE_H
#define BRIDGE_H 1

// Includes from SYCL matrix element calculations
#include <CL/sycl.hpp>
#include "mgOnGpuConfig.h" // for mgOnGpu::npar, mgOnGpu::np4
#include "CPPProcess.h" // for CPPProcess
#include "MatrixElementKernels.h" // for MatrixElementKernelHost, MatrixElementKernelDevice
#include "MemoryBuffers.h" // for HostBufferMomenta, DeviceBufferMomenta etc

#include <cassert>
#include <cmath>
#include <cstring>
#include <iostream>
#include <memory>

namespace mg5amcGpu
{

  // Forward declare transposition kernels
  // (Inverse transposition with template bool parameter F2C and "if constexpr" would require C++17, not available in cuda 11.1)

  template <typename T>
  void dev_transposeMomentaF2C( const T *in, T *out, const int evt, const size_t pos );

  template <typename T>
  void hst_transposeMomentaF2C( const T *in, T *out, const int evt );

  template <typename T>
  void hst_transposeMomentaC2F( const T *in, T *out, const int evt );

  // *****************************************************************************

  /**
   * A templated class for calling the SYCL matrix element calculations of
   * the event generation workflow. The template parameter is used for the
   * precision of the calculations in SYCL (float or double).
   *
   * The fortran momenta passed in are in the form of
   *   DOUBLE PRECISION P_MULTI(0:3, NEXTERNAL, NB_PAGE)
   * where the dimensions are <# momenta>, <# of particles>, <# events>
   *
   * FIXME: eventually gpu sequence should take double* (not T*) momenta/MEs.
   * This would allow using double in MadEvent Fortran and float in SYCL sigmaKin.
   * This would require significant changes to the "--bridge" test in check_sa.
   * A mixing of float and double is not yet possible in the SYCL code.
   * For the moment, use T* in gpu sequence (emulate using float or double throughout).
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
    Bridge( int evt, int par, int mom, size_t d_id, int /*dummy*/str=0, int /*dummy*/ncomb=0 );

    /**
     * set the gpublocks and gputhreads for the gpusequence - throws if evnt != gpublocks*gputhreads
     * (this is needed for BridgeKernel tests rather than for actual production use in Fortran)
     *
     * @param gpublocks number of gpublocks
     * @param gputhreads number of gputhreads
     */
    void set_gpugrid( const int gpublocks, const int gputhreads );

    /**
     * sequence to be executed for the SYCL matrix element calculation
     *
     * @param momenta memory address of the input 4-momenta
     * @param mes memory address of the output matrix elements
     * @param goodHelOnly quit after computing good helicities
     */
    void gpu_sequence( const T *momenta, T *mes, const bool goodHelOnly=false );

  private:

    int m_evt;                 ///< number of events
    //int m_part;                ///< number of particles per event [NO LONGER NEEDED!]
    //int m_mome;                ///< number of momenta per particle (usually 4) [NO LONGER NEEDED!]
    //int m_strd;                ///< stride length of the AOSOA structure [NO LONGER NEEDED!]
    //int m_ncomb;               ///< number of good helicities [NO LONGER NEEDED!]
    bool m_goodHelsCalculated; ///< have the good helicities been calculated?

    int m_gputhreads;          ///< number of gpu threads (default set from number of events, can be modified)
    int m_gpublocks;           ///< number of gpu blocks (default set from number of events, can be modified)

    size_t m_d_id;
    std::vector<sycl::device> m_devices;
    sycl::queue m_q;

    mg5amcGpu::DeviceBufferMomenta m_devMomentaF;
    mg5amcGpu::DeviceBufferMomenta m_devMomentaC;
    mg5amcGpu::DeviceBufferMatrixElements m_devMEsC;
    mg5amcGpu::MatrixElementKernelDevice m_devMek;
    
    mg5amcGpu::CPPProcess m_process;

  };

  // *****************************************************************************

  //
  // Implementations of class Bridge member functions
  //

  template <typename T>
  Bridge<T>::Bridge( int evnt, int part, int mome, size_t d_id, int /*strd*/, int /*ncomb*/ )
    : m_evt( evnt )
      //, m_part( part )
      //, m_mome( mome )
      //, m_strd( strd )
      //, m_ncomb( ncomb )
    , m_d_id( d_id )
    , m_devices( sycl::device::get_devices() )
    , m_q( sycl::queue(m_devices[m_d_id]) )
    , m_goodHelsCalculated( false )
    , m_gputhreads( 256 ) // default number of gpu threads
    , m_gpublocks( ceil(double(m_evt)/m_gputhreads) ) // this ensures m_evt <= m_gpublocks*m_gputhreads
    , m_devMomentaF( evnt, m_q)
    , m_devMomentaC( evnt, m_q)
    , m_devMEsC( evnt, m_q)
    , m_devMek( m_devMomentaC, m_devMEsC, m_q, m_gpublocks, m_gputhreads )
    , m_process( 1, m_gpublocks, m_gputhreads, false )
  {
    if ( part != mgOnGpu::npar ) throw std::runtime_error( "Bridge constructor: npar mismatch" );
    if ( mome != mgOnGpu::np4 ) throw std::runtime_error( "Bridge constructor: np4 mismatch" );
    std::cout << "WARNING! Instantiate Bridge (nevt=" << m_evt << ", gpublocks=" << m_gpublocks << ", gputhreads=" << m_gputhreads
              << ", gpublocks*gputhreads=" << m_gpublocks*m_gputhreads << ")" << std::endl;
    m_process.initProc( "../../Cards/param_card.dat" );
    m_devMek.setDeviceArrays( m_process.get_tHel_ptr(), m_process.get_tIPC_ptr(), m_process.get_tIPD_ptr() );
  }

  template <typename T>
  void Bridge<T>::set_gpugrid(const int gpublocks, const int gputhreads)
  {
    if ( m_evt != gpublocks*gputhreads )
      throw std::runtime_error( "Bridge: gpublocks*gputhreads must equal m_evt in set_gpugrid" );
    m_gpublocks = gpublocks;
    m_gputhreads = gputhreads;
    std::cout << "WARNING! Set grid in Bridge (nevt=" << m_evt << ", gpublocks=" << m_gpublocks << ", gputhreads=" << m_gputhreads
              << ", gpublocks*gputhreads=" << m_gpublocks*m_gputhreads << ")" << std::endl;
    m_devMek.setGrid( m_gpublocks, m_gputhreads );
  }

  template <typename T>
  void Bridge<T>::gpu_sequence( const T *momenta, T *mes, const bool goodHelOnly )
  {
    constexpr int neppM = mgOnGpu::neppM;
    if constexpr ( neppM == 1 ) // needs c++17
    {
      m_q.memcpy( m_devMomentaC.data(), momenta, m_devMomentaC.bytes() ).wait();
    }
    else
    {
      m_q.memcpy( m_devMomentaF.data(), momenta, m_devMomentaF.bytes() ).wait();
      const int thrPerEvt = mgOnGpu::npar * mgOnGpu::np4; // AV: transpose alg does 1 element per thread (NOT 1 event per thread)
      //const int thrPerEvt = 1; // AV: try new alg with 1 event per thread... this seems slower
      m_q.submit([&](sycl::handler& cgh) {
          cgh.parallel_for_work_group(sycl::range<1>{m_gpublocks}, sycl::range<1>{m_gputhreads}, ([=](sycl::group<1> wGroup) {
              wGroup.parallel_for_work_item([&](sycl::h_item<1> index) {
                  size_t ievt = index.get_global_id(0);
                  dev_transposeMomentaF2C( m_devMomentaF.data(), m_devMomentaC.data(), m_evt, ievt );
              });
          }));
      });
      m_q.wait();

    }
    if ( !m_goodHelsCalculated )
    {
      m_devMek.computeGoodHelicities();
      m_goodHelsCalculated = true;
    }
    if ( goodHelOnly ) return;
    m_devMek.computeMatrixElements();
    m_q.memcpy( mes, m_devMEsC.data(), m_devMEsC.bytes() ).wait();
  }

  // *****************************************************************************

  //
  // Implementations of transposition functions
  //

  template <typename T>
  void dev_transposeMomentaF2C( const T *in, T *out, const int evt, const size_t pos )
  {
    constexpr bool oldImplementation = true; // default: use old implementation
    if constexpr ( oldImplementation )
    {
      // Stefan's initial implementation
      constexpr int part = mgOnGpu::npar;
      constexpr int mome = mgOnGpu::np4;
      constexpr int strd = mgOnGpu::neppM;
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
    }
    else
    {
      // AV attempt another implementation with 1 event per thread: this seems slower...
      // C-style: AOSOA[npagM][npar][np4][neppM]
      // F-style: AOS[nevt][npar][np4]
      constexpr int npar = mgOnGpu::npar;
      constexpr int np4 = mgOnGpu::np4;
      constexpr int neppM = mgOnGpu::neppM;
      assert( evt%neppM == 0 ); // number of events is not a multiple of neppM???
      int ievt = pos;
      int ipagM = ievt/neppM;
      int ieppM = ievt%neppM;
      for ( int ip4=0; ip4<np4; ip4++ )
        for ( int ipar=0; ipar<npar; ipar++ )
        {
          int cpos = ipagM*npar*np4*neppM + ipar*np4*neppM + ip4*neppM + ieppM;
          int fpos = ievt*npar*np4 + ipar*np4 + ip4;
          out[cpos] = in[fpos]; // F2C (Fortran to C)
        }
    }
  }

  template <typename T, bool F2C>
  void hst_transposeMomenta( const T *in, T *out, const int evt )
  {
    constexpr bool oldImplementation = false; // default: use new implementation
    if constexpr ( oldImplementation )
    {
      // Stefan's initial implementation
      constexpr int part = mgOnGpu::npar;
      constexpr int mome = mgOnGpu::np4;
      constexpr int strd = mgOnGpu::neppM;
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
        if constexpr ( F2C ) out[pos] = in[inpos]; // F2C (Fortran to C)
        else out[inpos] = in[pos]; // C2F (C to Fortran)
      }
    }
    else
    {
      // AV attempt another implementation: this is slightly faster (better c++ pipelining?)
      // [NB! this is not a transposition, it is an AOS to AOSOA conversion: if neppM=1, a memcpy is enough]
      // C-style: AOSOA[npagM][npar][np4][neppM]
      // F-style: AOS[nevt][npar][np4]
      constexpr int npar = mgOnGpu::npar;
      constexpr int np4 = mgOnGpu::np4;
      constexpr int neppM = mgOnGpu::neppM;
      if constexpr ( neppM == 1 ) // needs c++17 and cuda >=11.2 (#333)
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
                if constexpr ( F2C ) out[cpos] = in[fpos]; // F2C (Fortran to C)
                else out[fpos] = in[cpos]; // C2F (C to Fortran)
              }
      }
    }
  }

  template <typename T>
  void hst_transposeMomentaF2C( const T *in, T *out, const int evt )
  {
    constexpr bool F2C = true;
    hst_transposeMomenta<T, F2C>( in, out, evt );
  }

  template <typename T>
  void hst_transposeMomentaC2F( const T *in, T *out, const int evt )
  {
    constexpr bool F2C = false;
    hst_transposeMomenta<T, F2C>( in, out, evt );
  }

}

#endif // BRIDGE_H
