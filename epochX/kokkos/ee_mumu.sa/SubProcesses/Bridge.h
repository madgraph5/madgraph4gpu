#ifndef BRIDGE_H
#define BRIDGE_H 1

#include "Kokkos_Core.hpp"
#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"

#include "CPPProcess.h"           // for CPPProcess

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <iostream>
#include <memory>
#include <type_traits>

namespace mg5amcKokkos
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
   * A templated class for calling the KOKKOS device matrix element calculations of the event generation workflow.
   * The FORTRANFPTYPE template parameter indicates the precision of the Fortran momenta from MadEvent (float or double).
   * The precision of the matrix element calculation is hardcoded in the fptype typedef in KOKKOS device.
   *
   * The Fortran momenta passed in are in the form of
   *   DOUBLE PRECISION P_MULTI(0:3, NEXTERNAL, NB_PAGE_LOOP)
   * where the dimensions are <np4F(#momenta)>, <nparF(#particles)>, <nevtF(#events)>.
   * In memory, this is stored in a way that C reads as an array P_MULTI[nevtF][nparF][np4F].
   * The KOKKOS device momenta are stored as an array[npagM][npar][np4][neppM] with nevt=npagM*neppM.
   * The Bridge is configured to store nevt==nevtF events in KOKKOS device.
   * It also checks that Fortran and C++ parameters match, nparF==npar and np4F==np4.
   *
   * The cpu/gpu sequences take FORTRANFPTYPE* (not fptype*) momenta/MEs.
   * This allows mixing double in MadEvent Fortran with float in KOKKOS device sigmaKin.
   * In the fcheck_sa.f test, Fortran uses double while KOKKOS device may use double or float.
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
    virtual ~Bridge() {}

    // Delete copy/move constructors and assignment operators
    Bridge( const Bridge& ) = delete;
    Bridge( Bridge&& ) = delete;
    Bridge& operator=( const Bridge& ) = delete;
    Bridge& operator=( Bridge&& ) = delete;

    /**
     * Set the league_size and team_size for the gpusequence - throws if evnt != league_size*team_size
     * (this is needed for BridgeKernel tests rather than for actual production use in Fortran)
     *
     * @param league_size number of league_size
     * @param team_size number of team_size
     */
    void set_gpugrid( const int league_size, const int team_size );

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

    int m_team_size; // number of gpu threads (default set from number of events, can be modified)
    int m_league_size;  // number of gpu blocks (default set from number of events, can be modified)
    Kokkos::View<fptype***> m_devMomenta;
    Kokkos::View<fptype*> m_devGs;
    Kokkos::View<fptype*> m_devMEs;
    Kokkos::View<fptype*,Kokkos::HostSpace> m_hstMEs;
    Kokkos::View<bool*> m_devIsGoodHel;
    Kokkos::View<short** > m_devcHel;
    Kokkos::View<fptype*> m_devcIPC;
    Kokkos::View<fptype*> m_devcIPD;
    Kokkos::View<int*> m_devNGoodHel; 
    Kokkos::View<cxtype*> m_dev_independent_couplings;
    Kokkos::View<fptype*> m_dev_independent_parameters;
    
    //static constexpr int s_team_sizemin = 16; // minimum number of gpu threads (TEST VALUE FOR MADEVENT)
    static constexpr int s_team_sizemin = 32; // minimum number of gpu threads (DEFAULT)
  };

  //--------------------------------------------------------------------------
  //
  // Forward declare transposition methods
  //

  // template<typename Tin, typename Tout>
  // void dev_transposeMomentaF2C( Tout* __restrict__ out, const Tin* __restrict__ in, size_t pos, const unsigned int nevt );

  // template<typename Tin, typename Tout>
  // void hst_transposeMomentaF2C( const Tin* __restrict__ in, Tout* __restrict__ out, const unsigned int nevt );

  // template<typename Tin, typename Tout>
  // void hst_transposeMomentaC2F( const Tin* __restrict__ in, Tout* __restrict__ out, const unsigned int nevt );

  //--------------------------------------------------------------------------
  //
  // Implementations of member functions of class Bridge
  //

  template<typename FORTRANFPTYPE>
  Bridge<FORTRANFPTYPE>::Bridge( unsigned int nevtF, unsigned int nparF, unsigned int np4F )
    : m_nevt( nevtF )
    , m_goodHelsCalculated( false )
    , m_team_size( 256 )                  // default number of gpu threads
    , m_league_size( m_nevt / m_team_size ) // this ensures m_nevt <= m_league_size*m_team_size
    , m_devMomenta( Kokkos::ViewAllocateWithoutInitializing("m_devMomenta"), m_nevt, mgOnGpu::npar, mgOnGpu::np4 )
    , m_devGs( Kokkos::ViewAllocateWithoutInitializing("m_devGs"), m_nevt )
    , m_devMEs( Kokkos::ViewAllocateWithoutInitializing("m_devMEs"), m_nevt )
    , m_hstMEs( Kokkos::ViewAllocateWithoutInitializing("m_hstMEs"), m_nevt )
    , m_devIsGoodHel( Kokkos::ViewAllocateWithoutInitializing("m_devIsGoodHel"), mgOnGpu::ncomb )
    , m_devNGoodHel( Kokkos::ViewAllocateWithoutInitializing("m_devNGoodHel"),1) 
    , m_dev_independent_couplings(Kokkos::ViewAllocateWithoutInitializing("m_dev_independent_couplings"), independentCouplings::nicoup )
    , m_dev_independent_parameters(Kokkos::ViewAllocateWithoutInitializing("m_dev_independent_parameters"), mgOnGpu::nparams )
  {
    if( nparF != mgOnGpu::npar ) throw std::runtime_error( "Bridge constructor: npar mismatch" );
    if( np4F != mgOnGpu::np4 ) throw std::runtime_error( "Bridge constructor: np4 mismatch" );
    if( ( m_nevt < s_team_sizemin ) || ( m_nevt % s_team_sizemin != 0 ) )
      throw std::runtime_error( "Bridge constructor: nevt should be a multiple of " + std::to_string( s_team_sizemin ) );
    while( m_nevt != m_league_size * m_team_size )
    {
      m_team_size /= 2;
      if( m_team_size < s_team_sizemin )
        throw std::logic_error( "Bridge constructor: FIXME! cannot choose team_size" ); // this should never happen!
      m_league_size = m_nevt / m_team_size;
    }
    std::cout << "WARNING! Instantiate device Bridge (nevt=" << m_nevt << ", league_size=" << m_league_size << ", team_size=" << m_team_size
              << ", league_size*team_size=" << m_league_size * m_team_size << ")" << std::endl;

    CPPProcess process( false, false );
    process.initProc( "../../Cards/param_card.dat" );
  }

  template<typename FORTRANFPTYPE>
  void Bridge<FORTRANFPTYPE>::set_gpugrid( const int league_size, const int team_size )
  {
    if( m_nevt != league_size * team_size )
      throw std::runtime_error( "Bridge: league_size*team_size must equal m_nevt in set_gpugrid" );
    m_league_size = league_size;
    m_team_size = team_size;
    std::cout << "WARNING! Set grid in Bridge (nevt=" << m_nevt << ", league_size=" << m_league_size << ", team_size=" << m_team_size
              << ", league_size*team_size=" << m_league_size * m_team_size << ")" << std::endl;
  }

  template<typename FORTRANFPTYPE>
  void Bridge<FORTRANFPTYPE>::gpu_sequence( const FORTRANFPTYPE* __restrict__ momenta,
                                            const FORTRANFPTYPE* __restrict__ gs,
                                            FORTRANFPTYPE* __restrict__ mes,
                                            const unsigned int channelId,
                                            const bool goodHelOnly )
  {
    static constexpr int np4 =  mgOnGpu::np4;
    static constexpr int nparf = mgOnGpu::nparf;
    static constexpr int npar = mgOnGpu::npar;
    static constexpr int ncomb = mgOnGpu::ncomb;

    // TC: warning that I am assuming incoming momenta from FORTRAN
    // TC: are already in the format [ievt][ipar][i4]

    // TC: Might need this -> if constexpr( neppM == 1 && std::is_same_v<FORTRANFPTYPE, fptype> )
    Kokkos::View<FORTRANFPTYPE const***,Kokkos::HostSpace> hst_momenta(momenta,m_nevt,npar,np4);
    // Kokkos::resize(hst_momenta,m_nevt,npar,np4);
    Kokkos::deep_copy(m_devMomenta,hst_momenta);

    Kokkos::View<FORTRANFPTYPE const *,Kokkos::HostSpace> hstGs(gs,m_nevt);
    Kokkos::deep_copy(m_devGs,hstGs);

    if( !m_goodHelsCalculated )
    {
      sigmaKin_setup(m_devMomenta,m_devIsGoodHel, m_devNGoodHel,
                     m_devGs, m_dev_independent_couplings,
                     m_dev_independent_parameters, 
                     m_league_size, m_team_size);

      m_goodHelsCalculated = true;
    }
    if( goodHelOnly ) return;

    sigmaKin(m_devMomenta,m_devIsGoodHel,m_devNGoodHel,
            m_devGs, m_dev_independent_couplings,
            m_dev_independent_parameters, 
            m_league_size, m_team_size,
            m_devMEs);

    Kokkos::deep_copy(m_hstMEs,m_devMEs);

  }

  //--------------------------------------------------------------------------
  //
  // Implementations of transposition methods
  // - FORTRAN arrays: P_MULTI(0:3, NEXTERNAL, NB_PAGE_LOOP) ==> p_multi[nevtF][nparF][np4F] in C++ (AOS)
  // - C++ array: momenta[npagM][npar][np4][neppM] with nevt=npagM*neppM (AOSOA)
  //

  // template<typename Tin, typename Tout>
  // void dev_transposeMomentaF2C( Tout* __restrict__ out, const Tin* __restrict__ in, size_t pos, const unsigned int nevt )
  // {
  //   constexpr bool oldImplementation = true; // default: use old implementation
  //   if constexpr( oldImplementation )
  //   {
  //     // SR initial implementation
  //     constexpr int part = mgOnGpu::npar;
  //     constexpr int mome = mgOnGpu::np4;
  //     constexpr int strd = mgOnGpu::neppM;
  //     unsigned int arrlen = nevt * part * mome;
  //     if( pos < arrlen )
  //     {
  //       int page_i = pos / ( strd * mome * part );
  //       int rest_1 = pos % ( strd * mome * part );
  //       int part_i = rest_1 / ( strd * mome );
  //       int rest_2 = rest_1 % ( strd * mome );
  //       int mome_i = rest_2 / strd;
  //       int strd_i = rest_2 % strd;
  //       int inpos =
  //         ( page_i * strd + strd_i ) // event number
  //           * ( part * mome )        // event size (pos of event)
  //         + part_i * mome            // particle inside event
  //         + mome_i;                  // momentum inside particle
  //       out[pos] = in[inpos];        // F2C (Fortran to C)
  //     }
  //   }
  //   else
  //   {
  //     // AV attempt another implementation with 1 event per thread: this seems slower...
  //     // F-style: AOS[nevtF][nparF][np4F]
  //     // C-style: AOSOA[npagM][npar][np4][neppM] with nevt=npagM*neppM
  //     constexpr int npar = mgOnGpu::npar;
  //     constexpr int np4 = mgOnGpu::np4;
  //     constexpr int neppM = mgOnGpu::neppM;
  //     assert( nevt % neppM == 0 ); // number of events is not a multiple of neppM???
  //     size_t ievt = pos;
  //     int ipagM = ievt / neppM;
  //     int ieppM = ievt % neppM;
  //     for( int ip4 = 0; ip4 < np4; ip4++ )
  //       for( int ipar = 0; ipar < npar; ipar++ )
  //       {
  //         int cpos = ipagM * npar * np4 * neppM + ipar * np4 * neppM + ip4 * neppM + ieppM;
  //         int fpos = ievt * npar * np4 + ipar * np4 + ip4;
  //         out[cpos] = in[fpos]; // F2C (Fortran to C)
  //       }
  //   }
  // }

  // template<typename Tin, typename Tout, bool F2C>
  // void hst_transposeMomenta( const Tin* __restrict__ in, Tout* __restrict__ out, const unsigned int nevt )
  // {
  //   constexpr bool oldImplementation = false; // default: use new implementation
  //   if constexpr( oldImplementation )
  //   {
  //     // SR initial implementation
  //     constexpr unsigned int part = mgOnGpu::npar;
  //     constexpr unsigned int mome = mgOnGpu::np4;
  //     constexpr unsigned int strd = mgOnGpu::neppM;
  //     unsigned int arrlen = nevt * part * mome;
  //     for( unsigned int pos = 0; pos < arrlen; ++pos )
  //     {
  //       unsigned int page_i = pos / ( strd * mome * part );
  //       unsigned int rest_1 = pos % ( strd * mome * part );
  //       unsigned int part_i = rest_1 / ( strd * mome );
  //       unsigned int rest_2 = rest_1 % ( strd * mome );
  //       unsigned int mome_i = rest_2 / strd;
  //       unsigned int strd_i = rest_2 % strd;
  //       unsigned int inpos =
  //         ( page_i * strd + strd_i ) // event number
  //           * ( part * mome )        // event size (pos of event)
  //         + part_i * mome            // particle inside event
  //         + mome_i;                  // momentum inside particle
  //       if constexpr( F2C )          // needs c++17
  //         out[pos] = in[inpos];      // F2C (Fortran to C)
  //       else
  //         out[inpos] = in[pos]; // C2F (C to Fortran)
  //     }
  //   }
  //   else
  //   {
  //     // AV attempt another implementation: this is slightly faster (better c++ pipelining?)
  //     // [NB! this is not a transposition, it is an AOS to AOSOA conversion: if neppM=1, a memcpy is enough]
  //     // F-style: AOS[nevtF][nparF][np4F]
  //     // C-style: AOSOA[npagM][npar][np4][neppM] with nevt=npagM*neppM
  //     constexpr unsigned int npar = mgOnGpu::npar;
  //     constexpr unsigned int np4 = mgOnGpu::np4;
  //     constexpr unsigned int neppM = mgOnGpu::neppM;
  //     if constexpr( neppM == 1 && std::is_same_v<Tin, Tout> )
  //     {
  //       memcpy( out, in, nevt * npar * np4 * sizeof( Tin ) );
  //     }
  //     else
  //     {
  //       const unsigned int npagM = nevt / neppM;
  //       assert( nevt % neppM == 0 ); // number of events is not a multiple of neppM???
  //       for( unsigned int ipagM = 0; ipagM < npagM; ipagM++ )
  //         for( unsigned int ip4 = 0; ip4 < np4; ip4++ )
  //           for( unsigned int ipar = 0; ipar < npar; ipar++ )
  //             for( unsigned int ieppM = 0; ieppM < neppM; ieppM++ )
  //             {
  //               unsigned int ievt = ipagM * neppM + ieppM;
  //               unsigned int cpos = ipagM * npar * np4 * neppM + ipar * np4 * neppM + ip4 * neppM + ieppM;
  //               unsigned int fpos = ievt * npar * np4 + ipar * np4 + ip4;
  //               if constexpr( F2C )
  //                 out[cpos] = in[fpos]; // F2C (Fortran to C)
  //               else
  //                 out[fpos] = in[cpos]; // C2F (C to Fortran)
  //             }
  //     }
  //   }
  // }

  // template<typename Tin, typename Tout>
  // void hst_transposeMomentaF2C( const Tin* __restrict__ in, Tout* __restrict__ out, const unsigned int nevt )
  // {
  //   constexpr bool F2C = true;
  //   hst_transposeMomenta<Tin, Tout, F2C>( in, out, nevt );
  // }

  // template<typename Tin, typename Tout>
  // void hst_transposeMomentaC2F( const Tin* __restrict__ in, Tout* __restrict__ out, const unsigned int nevt )
  // {
  //   constexpr bool F2C = false;
  //   hst_transposeMomenta<Tin, Tout, F2C>( in, out, nevt );
  // }

  //--------------------------------------------------------------------------
}
#endif // BRIDGE_H
