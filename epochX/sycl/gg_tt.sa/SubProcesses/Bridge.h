#ifndef BRIDGE_H
#define BRIDGE_H 1

#include <CL/sycl.hpp>
#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"
#include "mgOnGpuVectors.h"

#include "CPPProcess.h"           // for CPPProcess

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <iostream>
#include <memory>
#include <type_traits>

namespace mg5amcGpu
{
  //--------------------------------------------------------------------------
  /**
   * A base class for a class whose pointer is passed between Fortran and C++.
   * This is not really necessary, but it allows minimal type checks on all such pointers.
   */
  struct CppObjectInFortran
  {
    CppObjectInFortran() {}
    virtual ~CppObjectInFortran() {}
  };

  //--------------------------------------------------------------------------
  /**
   * A templated class for calling the SYCL device matrix element calculations of the event generation workflow.
   * The FORTRANFPTYPE template parameter indicates the precision of the Fortran momenta from MadEvent (float or double).
   * The precision of the matrix element calculation is hardcoded in the fptype typedef in SYCL device.
   *
   * The Fortran momenta passed in are in the form of
   *   DOUBLE PRECISION P_MULTI(0:3, NEXTERNAL, NB_PAGE_LOOP)
   * where the dimensions are <np4F(#momenta)>, <nparF(#particles)>, <nevtF(#events)>.
   * In memory, this is stored in a way that C reads as an array P_MULTI[nevtF][nparF][np4F].
   * The SYCL device momenta are stored as an array[npagM][npar][np4][neppM] with nevt=npagM*neppM.
   * The Bridge is configured to store nevt==nevtF events in SYCL device.
   * It also checks that Fortran and C++ parameters match, nparF==npar and np4F==np4.
   *
   * The cpu/gpu sequences take FORTRANFPTYPE* (not fptype*) momenta/MEs.
   * This allows mixing double in MadEvent Fortran with float in SYCL device sigmaKin.
   * In the fcheck_sa.f test, Fortran uses double while SYCL device may use double or float.
   * In the check_sa "--bridge" test, everything is implemented in fptype (double or float).
   */
  template<typename FORTRANFPTYPE>
  class Bridge final : public CppObjectInFortran
  {
  public:
    /**
     * Constructor
     *
     * @param nevtF (NB_PAGE_LOOP, vector.inc) number of events in Fortran array loops (NB_PAGE_LOOP <= NB_PAGE_MAX)
     * @param nparF (NEXTERNAL, nexternal.inc) number of external particles in Fortran arrays (KEPT FOR SANITY CHECKS ONLY)
     * @param np4F number of momenta components, usually 4, in Fortran arrays (KEPT FOR SANITY CHECKS ONLY)
     */
    Bridge( unsigned int nevtF, unsigned int nparF, unsigned int np4F );

    /**
     * Destructor
     */
    ~Bridge() {
        sycl::free( m_devMomentaF  , m_q);
        sycl::free( m_devMomentaC  , m_q);
        sycl::free( m_devMEsC      , m_q);
        sycl::free( m_hstMEsC      , m_q);
        sycl::free( m_devIsGoodHel , m_q);
        sycl::free( m_hstIsGoodHel , m_q);
        sycl::free( m_devcHel      , m_q);
        sycl::free( m_devcIPC      , m_q);
        sycl::free( m_devcIPD      , m_q);
        sycl::free( m_devcNGoodHel , m_q); 
        sycl::free( m_devcGoodHel  , m_q); 
    }

    // Delete copy/move constructors and assignment operators
    Bridge( const Bridge& ) = delete;
    Bridge( Bridge&& ) = delete;
    Bridge& operator=( const Bridge& ) = delete;
    Bridge& operator=( Bridge&& ) = delete;

    /**
     * Set the gpublocks and gputhreads for the gpusequence - throws if evnt != gpublocks*gputhreads
     * (this is needed for BridgeKernel tests rather than for actual production use in Fortran)
     *
     * @param gpublocks number of gpublocks
     * @param gputhreads number of gputhreads
     */
    void set_gpugrid( const int gpublocks, const int gputhreads );

    /**
     * Sequence to be executed for the matrix element calculation
     *
     * @param momenta the pointer to the input 4-momenta
     * @param gs the pointer to the input Gs (running QCD coupling constant alphas)
     * @param mes the pointer to the output matrix elements
     * @param channelId the Feynman diagram to enhance in multi-channel mode if 1 to n (disable multi-channel if 0)
     * @param goodHelOnly quit after computing good helicities?
     */
    void gpu_sequence( const FORTRANFPTYPE* __restrict__ momenta,
                       const FORTRANFPTYPE* __restrict__ gs,
                       FORTRANFPTYPE* __restrict__ mes,
                       const unsigned int channelId,
                       const bool goodHelOnly = false );

  private:
    unsigned int m_nevt;       // number of events
    bool m_goodHelsCalculated; // have the good helicities been calculated?

    std::vector<sycl::device> m_devices;
    sycl::queue m_q;

    int m_gputhreads; // number of gpu threads (default set from number of events, can be modified)
    int m_gpublocks;  // number of gpu blocks (default set from number of events, can be modified)
    sycl::malloc_device<fptype> m_devMomentaF;
    sycl::malloc_device<fptype> m_devMomentaC;
    //sycl::malloc_device<fptype> m_devGsC;
    //std::unique_ptr<fptype> m_hstGsC;
    sycl::malloc_device<fptype> m_devMEsC;
    sycl::malloc_host<fptype  > m_hstMEsC;
    sycl::malloc_device<bool  > m_devIsGoodHel;
    sycl::malloc_host<bool    > m_hstIsGoodHel;
    sycl::malloc_device<short > m_devcHel;
    sycl::malloc_device<fptype> m_devcIPC;
    sycl::malloc_device<fptype> m_devcIPD;
    sycl::malloc_device<int   > m_devcNGoodHel; 
    sycl::malloc_device<int   > m_devcGoodHel; 
    
    //static constexpr int s_gputhreadsmin = 16; // minimum number of gpu threads (TEST VALUE FOR MADEVENT)
    static constexpr int s_gputhreadsmin = 32; // minimum number of gpu threads (DEFAULT)
  };

  //--------------------------------------------------------------------------
  //
  // Forward declare transposition methods
  //

  template<typename Tin, typename Tout>
  void dev_transposeMomentaF2C( Tout* __restrict__ out, const Tin* __restrict__ in, size_t pos, const unsigned int nevt );

  template<typename Tin, typename Tout>
  void hst_transposeMomentaF2C( const Tin* __restrict__ in, Tout* __restrict__ out, const unsigned int nevt );

  template<typename Tin, typename Tout>
  void hst_transposeMomentaC2F( const Tin* __restrict__ in, Tout* __restrict__ out, const unsigned int nevt );

  //--------------------------------------------------------------------------
  //
  // Implementations of member functions of class Bridge
  //

  template<typename FORTRANFPTYPE>
  Bridge<FORTRANFPTYPE>::Bridge( unsigned int nevtF, unsigned int nparF, unsigned int np4F )
    : m_nevt( nevtF )
    , m_devices( sycl::device::get_devices() )
#ifdef MGONGPU_DEVICE_ID
    , m_q( sycl::queue(m_devices[MGONGPU_DEVICE_ID]) )
#else
    , m_q( sycl::queue(m_devices[0]) )
#endif
    , m_goodHelsCalculated( false )
    , m_gputhreads( 256 )                  // default number of gpu threads
    , m_gpublocks( m_nevt / m_gputhreads ) // this ensures m_nevt <= m_gpublocks*m_gputhreads
    , m_devMomentaF( m_nevt, m_q )
    , m_devMomentaC( m_nevt, m_q )
    //, m_devGsC( m_nevt )
    //, m_hstGsC( m_nevt )
    , m_devMEsC( m_nevt, m_q )
    , m_hstMEsC( m_nevt, m_q )
    , m_devIsGoodHel( mgOnGpu::ncomb, m_q )
    , m_hstIsGoodHel( mgOnGpu::ncomb, m_q )
    , m_devcHel( mgOnGpu::ncomb*mgOnGpu::npar, m_q )
    , m_devcIPC( mgOnGpu::ncouplingstimes2, m_q )
    , m_devcIPD( mgOnGpu::nparams, m_q )
    , m_devcNGoodHel( 1, m_q ) 
    , m_devcGoodHel( mgOnGpu::ncomb, m_q ) 
  {
    if( nparF != mgOnGpu::npar ) throw std::runtime_error( "Bridge constructor: npar mismatch" );
    if( np4F != mgOnGpu::np4 ) throw std::runtime_error( "Bridge constructor: np4 mismatch" );
    if( ( m_nevt < s_gputhreadsmin ) || ( m_nevt % s_gputhreadsmin != 0 ) )
      throw std::runtime_error( "Bridge constructor: nevt should be a multiple of " + std::to_string( s_gputhreadsmin ) );
    while( m_nevt != m_gpublocks * m_gputhreads )
    {
      m_gputhreads /= 2;
      if( m_gputhreads < s_gputhreadsmin )
        throw std::logic_error( "Bridge constructor: FIXME! cannot choose gputhreads" ); // this should never happen!
      m_gpublocks = m_nevt / m_gputhreads;
    }
    std::cout << "WARNING! Instantiate device Bridge (nevt=" << m_nevt << ", gpublocks=" << m_gpublocks << ", gputhreads=" << m_gputhreads
              << ", gpublocks*gputhreads=" << m_gpublocks * m_gputhreads << ")" << std::endl;

    Proc::CPPProcess process( 1, m_gpublocks, m_gputhreads, false, false );
    process.initProc( "../../Cards/param_card.dat" );

    m_q.memcpy( m_devcHel, process.get_tHel_ptr(), mgOnGpu::ncomb*mgOnGpu::npar*sizeof(short) );
    m_q.memcpy( m_devcIPC, process.get_tIPC_ptr(), mgOnGpu::ncouplingstimes2*sizeof(fptype) );
    m_q.memcpy( m_devcIPD, process.get_tIPD_ptr(), mgOnGpu::nparams*sizeof(fptype) ).wait();
  }

  template<typename FORTRANFPTYPE>
  void Bridge<FORTRANFPTYPE>::set_gpugrid( const int gpublocks, const int gputhreads )
  {
    if( m_nevt != gpublocks * gputhreads )
      throw std::runtime_error( "Bridge: gpublocks*gputhreads must equal m_nevt in set_gpugrid" );
    m_gpublocks = gpublocks;
    m_gputhreads = gputhreads;
    std::cout << "WARNING! Set grid in Bridge (nevt=" << m_nevt << ", gpublocks=" << m_gpublocks << ", gputhreads=" << m_gputhreads
              << ", gpublocks*gputhreads=" << m_gpublocks * m_gputhreads << ")" << std::endl;
  }

  template<typename FORTRANFPTYPE>
  void Bridge<FORTRANFPTYPE>::gpu_sequence( const FORTRANFPTYPE* __restrict__ momenta,
                                            const FORTRANFPTYPE* __restrict__ gs,
                                            FORTRANFPTYPE* __restrict__ mes,
                                            const unsigned int channelId,
                                            const bool goodHelOnly )
  {
    static constexpr int neppM = mgOnGpu::neppM;
    static constexpr int np4 =  mgOnGpu::np4;
    static constexpr int nparf = mgOnGpu::nparf;
    static constexpr int npar = mgOnGpu::npar;
    static constexpr int ncomb = mgOnGpu::ncomb;

    if constexpr( neppM == 1 && std::is_same_v<FORTRANFPTYPE, fptype> )
    {
      m_q.memcpy( m_devMomentaC, momenta, np4*npar*m_nevt*sizeof(fptype) ).wait();
    }
    else
    {
      m_q.memcpy( m_devMomentaF, momenta, np4*npar*m_nevt*sizeof(fptype) ).wait();
      static constexpr int thrPerEvt = npar * np4; // AV: transpose alg does 1 element per thread (NOT 1 event per thread)

      //const int thrPerEvt = 1; // AV: try new alg with 1 event per thread... this seems slower
      m_q.submit([&](sycl::handler& cgh) {
          cgh.parallel_for_work_group(sycl::range<1>{m_gpublocks * thrPerEvt}, sycl::range<1>{m_gputhreads}, ([=](sycl::group<1> wGroup) {
              wGroup.parallel_for_work_item([&](sycl::h_item<1> index) {
                  size_t ievt = index.get_global_id(0);
                  dev_transposeMomentaF2C( m_devMomentaC, m_devMomentaF, ievt, m_nevt );
              });
          }));
      });
      m_q.wait();
    }
    //std::copy( gs, gs + m_nevt, m_hstGsC );
    //m_q.memcpy( m_devGsC, m_hstGsC, m_nevt*sizeof(fptype) ).wait();
    if( !m_goodHelsCalculated )
    {
      m_q.submit([&](sycl::handler& cgh) {
          cgh.parallel_for_work_group(sycl::range<1>{m_gpublocks}, sycl::range<1>{m_gputhreads}, ([=](sycl::group<1> wGroup) {
              wGroup.parallel_for_work_item([&](sycl::h_item<1> index) {
                  size_t ievt = index.get_global_id(0);
                  const int ipagM = ievt/neppM; // #eventpage in this iteration
                  const int ieppM = ievt%neppM; // #event in the current eventpage in this iteration
                  Proc::sigmaKin_getGoodHel( m_devMomentaC + ipagM * npar * np4 * neppM + ieppM, m_devIsGoodHel, m_devcHel, m_devcIPC, m_devcIPD );
              });
          }));
      });
      m_q.wait();

      m_q.memcpy(m_hstIsGoodHel, m_devIsGoodHel, ncomb*sizeof(bool)).wait();

      int goodHel[mgOnGpu::ncomb] = {0};
      int nGoodHel = Proc::sigmaKin_setGoodHel( m_hstIsGoodHel, goodHel );

      m_q.memcpy( m_devcNGoodHel, &nGoodHel, sizeof(int) ).wait();
      m_q.memcpy( m_devcGoodHel, goodHel, ncomb*sizeof(int) ).wait();
      m_goodHelsCalculated = true;
    }
    if( goodHelOnly ) return;

    m_q.submit([&](sycl::handler& cgh) {
        cgh.parallel_for_work_group(sycl::range<1>{m_gpublocks}, sycl::range<1>{m_gputhreads}, ([=](sycl::group<1> wGroup) {
            wGroup.parallel_for_work_item([&](sycl::h_item<1> index) {
                size_t ievt = index.get_global_id(0);
                const int ipagM = ievt/neppM;
                const int ieppM = ievt%neppM;
                m_devMEsC[ievt] = Proc::sigmaKin( m_devMomenta + ipagM * npar * np4 * neppM + ieppM, m_devcHel, m_devcIPC, m_devcIPD, m_devcNGoodHel, m_devcGoodHel );
            });
        }));
    });
    m_q.wait();

    if constexpr( std::is_same_v<FORTRANFPTYPE, fptype> )
    {
      m_q.memcpy( mes, m_devMEsC, m_nevt*sizeof(fptype) ).wait();
      //flagAbnormalMEs( mes, m_nevt );
    }
    else
    {
      m_q.memcpy( m_hstMEsC, m_devMEsC, m_nevt*sizeof(fptype) ).wait();
      //flagAbnormalMEs( m_hstMEsC, m_nevt );
      std::copy( m_hstMEsC, m_hstMEsC + m_nevt, mes );
    }

  }

  //--------------------------------------------------------------------------
  //
  // Implementations of transposition methods
  // - FORTRAN arrays: P_MULTI(0:3, NEXTERNAL, NB_PAGE_LOOP) ==> p_multi[nevtF][nparF][np4F] in C++ (AOS)
  // - C++ array: momenta[npagM][npar][np4][neppM] with nevt=npagM*neppM (AOSOA)
  //

  template<typename Tin, typename Tout>
  void dev_transposeMomentaF2C( Tout* __restrict__ out, const Tin* __restrict__ in, size_t pos, const unsigned int nevt )
  {
    constexpr bool oldImplementation = true; // default: use old implementation
    if constexpr( oldImplementation )
    {
      // SR initial implementation
      constexpr int part = mgOnGpu::npar;
      constexpr int mome = mgOnGpu::np4;
      constexpr int strd = MemoryAccessMomenta::neppM;
      int arrlen = nevt * part * mome;
      if( pos < arrlen )
      {
        int page_i = pos / ( strd * mome * part );
        int rest_1 = pos % ( strd * mome * part );
        int part_i = rest_1 / ( strd * mome );
        int rest_2 = rest_1 % ( strd * mome );
        int mome_i = rest_2 / strd;
        int strd_i = rest_2 % strd;
        int inpos =
          ( page_i * strd + strd_i ) // event number
            * ( part * mome )        // event size (pos of event)
          + part_i * mome            // particle inside event
          + mome_i;                  // momentum inside particle
        out[pos] = in[inpos];        // F2C (Fortran to C)
      }
    }
    else
    {
      // AV attempt another implementation with 1 event per thread: this seems slower...
      // F-style: AOS[nevtF][nparF][np4F]
      // C-style: AOSOA[npagM][npar][np4][neppM] with nevt=npagM*neppM
      constexpr int npar = mgOnGpu::npar;
      constexpr int np4 = mgOnGpu::np4;
      constexpr int neppM = MemoryAccessMomenta::neppM;
      assert( nevt % neppM == 0 ); // number of events is not a multiple of neppM???
      size_t ievt = pos;
      int ipagM = ievt / neppM;
      int ieppM = ievt % neppM;
      for( int ip4 = 0; ip4 < np4; ip4++ )
        for( int ipar = 0; ipar < npar; ipar++ )
        {
          int cpos = ipagM * npar * np4 * neppM + ipar * np4 * neppM + ip4 * neppM + ieppM;
          int fpos = ievt * npar * np4 + ipar * np4 + ip4;
          out[cpos] = in[fpos]; // F2C (Fortran to C)
        }
    }
  }

  template<typename Tin, typename Tout, bool F2C>
  void hst_transposeMomenta( const Tin* __restrict__ in, Tout* __restrict__ out, const unsigned int nevt )
  {
    constexpr bool oldImplementation = false; // default: use new implementation
    if constexpr( oldImplementation )
    {
      // SR initial implementation
      constexpr unsigned int part = mgOnGpu::npar;
      constexpr unsigned int mome = mgOnGpu::np4;
      constexpr unsigned int strd = MemoryAccessMomenta::neppM;
      unsigned int arrlen = nevt * part * mome;
      for( unsigned int pos = 0; pos < arrlen; ++pos )
      {
        unsigned int page_i = pos / ( strd * mome * part );
        unsigned int rest_1 = pos % ( strd * mome * part );
        unsigned int part_i = rest_1 / ( strd * mome );
        unsigned int rest_2 = rest_1 % ( strd * mome );
        unsigned int mome_i = rest_2 / strd;
        unsigned int strd_i = rest_2 % strd;
        unsigned int inpos =
          ( page_i * strd + strd_i ) // event number
            * ( part * mome )        // event size (pos of event)
          + part_i * mome            // particle inside event
          + mome_i;                  // momentum inside particle
        if constexpr( F2C )          // needs c++17
          out[pos] = in[inpos];      // F2C (Fortran to C)
        else
          out[inpos] = in[pos]; // C2F (C to Fortran)
      }
    }
    else
    {
      // AV attempt another implementation: this is slightly faster (better c++ pipelining?)
      // [NB! this is not a transposition, it is an AOS to AOSOA conversion: if neppM=1, a memcpy is enough]
      // F-style: AOS[nevtF][nparF][np4F]
      // C-style: AOSOA[npagM][npar][np4][neppM] with nevt=npagM*neppM
      constexpr unsigned int npar = mgOnGpu::npar;
      constexpr unsigned int np4 = mgOnGpu::np4;
      constexpr unsigned int neppM = MemoryAccessMomenta::neppM;
      if constexpr( neppM == 1 && std::is_same_v<Tin, Tout> )
      {
        memcpy( out, in, nevt * npar * np4 * sizeof( Tin ) );
      }
      else
      {
        const unsigned int npagM = nevt / neppM;
        assert( nevt % neppM == 0 ); // number of events is not a multiple of neppM???
        for( unsigned int ipagM = 0; ipagM < npagM; ipagM++ )
          for( unsigned int ip4 = 0; ip4 < np4; ip4++ )
            for( unsigned int ipar = 0; ipar < npar; ipar++ )
              for( unsigned int ieppM = 0; ieppM < neppM; ieppM++ )
              {
                unsigned int ievt = ipagM * neppM + ieppM;
                unsigned int cpos = ipagM * npar * np4 * neppM + ipar * np4 * neppM + ip4 * neppM + ieppM;
                unsigned int fpos = ievt * npar * np4 + ipar * np4 + ip4;
                if constexpr( F2C )
                  out[cpos] = in[fpos]; // F2C (Fortran to C)
                else
                  out[fpos] = in[cpos]; // C2F (C to Fortran)
              }
      }
    }
  }

  template<typename Tin, typename Tout>
  void hst_transposeMomentaF2C( const Tin* __restrict__ in, Tout* __restrict__ out, const unsigned int nevt )
  {
    constexpr bool F2C = true;
    hst_transposeMomenta<Tin, Tout, F2C>( in, out, nevt );
  }

  template<typename Tin, typename Tout>
  void hst_transposeMomentaC2F( const Tin* __restrict__ in, Tout* __restrict__ out, const unsigned int nevt )
  {
    constexpr bool F2C = false;
    hst_transposeMomenta<Tin, Tout, F2C>( in, out, nevt );
  }

  //--------------------------------------------------------------------------
}
#endif // BRIDGE_H
