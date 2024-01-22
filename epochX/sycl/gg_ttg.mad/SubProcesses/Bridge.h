// Copyright (C) 2020-2023 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: S. Roiser (Nov 2021) for the MG5aMC CUDACPP plugin.
// Further modified by: S. Roiser, A. Valassi (2021-2023) for the MG5aMC CUDACPP plugin.
//==========================================================================
// Copyright (C) 2021-2023 Argonne National Laboratory.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Modified by: N. Nichols (2021-2023) for the MG5aMC SYCL plugin.
//==========================================================================

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

  
  // A base class encapsulating a memory buffer (not necessarily an event buffer)
  template<typename T>
  class buffer_base {
  protected:
    buffer_base( const size_t size, const bool on_device, sycl::queue q) : m_size( size ), m_data( nullptr ), m_on_device( on_device ), m_q( q ){}
    virtual ~buffer_base(){}
  public:
    T* data(){ return m_data; }
    const T* data() const{ return m_data; }
    sycl::queue get_queue(){ return m_q; }
    const sycl::queue get_queue() const{ return m_q; }
    T& operator[]( const size_t index ){ return m_data[index]; }
    const T& operator[]( const size_t index ) const { return m_data[index]; }
    size_t size() const{ return m_size; }
    size_t bytes() const{ return m_size * sizeof(T); }
    bool isOnDevice() const { return m_on_device; }
  protected:
    const size_t m_size;
    T* m_data;
    const bool m_on_device;
    sycl::queue m_q;
 };

  // A class encapsulating a SYCL device buffer
  template<typename T>
  class device_buffer : public buffer_base<T> {
  public:
    device_buffer( const size_t size, sycl::queue q ) : buffer_base<T>( size, true, q ) {
        this->m_data = sycl::malloc_device<T>(this->size(), this->m_q);
    }
    virtual ~device_buffer() {
        sycl::free( this->m_data, this->m_q );
    }
  };

  // A class encapsulating a SYCL host buffer
  template<typename T>
  class host_buffer : public buffer_base<T> {
  public:
    host_buffer( const size_t size, sycl::queue q ) : buffer_base<T>( size, false, q ) {
        this->m_data = sycl::malloc_host<T>(this->size(), this->m_q);
    }
    virtual ~host_buffer() {
        sycl::free( this->m_data, this->m_q );
    }
  };

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
    Bridge( size_t nevtF, size_t nparF, size_t np4F );

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
    void set_gpugrid( const size_t gpublocks, const size_t gputhreads );

    /**
     * Sequence to be executed for the matrix element calculation
     *
     * @param momenta the pointer to the input 4-momenta
     * @param gs the pointer to the input Gs (running QCD coupling constant alphas)
     * @param rndhel the pointer to the input random numbers for helicity selection
     * @param rndcol the pointer to the input random numbers for color selection
     * @param channelId the Feynman diagram to enhance in multi-channel mode if 1 to n (disable multi-channel if 0)
     * @param mes the pointer to the output matrix elements
     * @param goodHelOnly quit after computing good helicities?
     * @param selhel the pointer to the output selected helicities
     * @param selcol the pointer to the output selected colors
     */
    void gpu_sequence( const FORTRANFPTYPE* __restrict__ momenta,
                       const FORTRANFPTYPE* __restrict__ gs,
                       const FORTRANFPTYPE* __restrict__ rndhel,
                       const FORTRANFPTYPE* __restrict__ rndcol,
                       const size_t channelId,
                       FORTRANFPTYPE* __restrict__ mes,
                       int* __restrict__ selhel,
                       int* __restrict__ selcol,
                       const bool goodHelOnly = false );

  // Return the number of good helicities (-1 initially when they have not yet been calculated)
  int nGoodHel() const { return m_nGoodHel; }

  // Return the total number of helicities (expose cudacpp ncomb in the Bridge interface to Fortran)
  constexpr int nTotHel() const { return int(CPPPROCESS_NCOMB); }

  private:
    size_t m_nevt;       // number of events
    int m_nGoodHel;      // the number of good helicities (-1 initially when they have not yet been calculated)
    bool m_goodHelsCalculated; // have the good helicities been calculated?

    std::vector<sycl::device> m_devices;
    sycl::queue m_q;

    size_t m_gputhreads; // number of gpu threads (default set from number of events, can be modified)
    size_t m_gpublocks;  // number of gpu blocks (default set from number of events, can be modified)
    device_buffer<vector4  > m_devMomentaC;
    device_buffer<fptype_sv> m_devGsC;
    device_buffer<fptype_sv> m_devRndHel;
    device_buffer<fptype_sv> m_devRndCol;
    device_buffer<fptype_sv> m_devMEsC;
    device_buffer<int_sv   > m_devSelHel;
    device_buffer<int_sv   > m_devSelCol;
    host_buffer<fptype     > m_hstMEsC;
    host_buffer<long       > m_hstSelHel;
    host_buffer<long       > m_hstSelCol;
    device_buffer<bool     > m_devIsGoodHel;
    host_buffer<bool       > m_hstIsGoodHel;
    #ifndef MGONGPU_HARDCODE_PARAM
        device_buffer<cxtype> m_dev_independent_couplings;
        device_buffer<fptype> m_dev_independent_parameters;
    #endif
    device_buffer<size_t   > m_devcNGoodHel; 
    device_buffer<size_t   > m_devcGoodHel; 
    
    static constexpr size_t s_gputhreadsmin = 32; // minimum number of gpu threads (DEFAULT)
  };

  //--------------------------------------------------------------------------
  //
  // Forward declare transposition methods
  //

  template<typename FPTypeDst, typename FPTypeSrc>
  void hst_transposeMomenta(FPTypeDst* __restrict__ dst, const FPTypeSrc* __restrict__ src, const size_t N );

  //--------------------------------------------------------------------------
  //
  // Implementations of member functions of class Bridge
  //

  template<typename FORTRANFPTYPE>
  Bridge<FORTRANFPTYPE>::Bridge( size_t nevtF, size_t nparF, size_t np4F )
    : m_nevt( nevtF )
    , m_nGoodHel( -1 )
    , m_goodHelsCalculated( false )
    , m_devices( sycl::device::get_devices() )
    #ifdef MGONGPU_DEVICE_ID
        , m_q( sycl::queue(m_devices[MGONGPU_DEVICE_ID]) )
    #else
        , m_q( sycl::queue(m_devices[0]) )
    #endif
    , m_gputhreads( 256 )                  // default number of gpu threads
    #if MGONGPU_VEC_DIM > 1
        , m_gpublocks( m_nevt / m_gputhreads / MGONGPU_VEC_DIM) // this ensures m_nevt <= m_gpublocks*m_gputhreads
    #else
        , m_gpublocks( m_nevt / m_gputhreads ) // this ensures m_nevt <= m_gpublocks*m_gputhreads
    #endif
    , m_devMomentaC( m_nevt*CPPPROCESS_NPAR/MGONGPU_VEC_DIM, m_q )
    , m_devGsC( m_nevt/MGONGPU_VEC_DIM, m_q )
    , m_devRndHel( m_nevt/MGONGPU_VEC_DIM, m_q )
    , m_devRndCol( m_nevt/MGONGPU_VEC_DIM, m_q )
    , m_devMEsC( m_nevt/MGONGPU_VEC_DIM, m_q )
    , m_devSelHel( m_nevt/MGONGPU_VEC_DIM, m_q )
    , m_devSelCol( m_nevt/MGONGPU_VEC_DIM, m_q )
    , m_hstMEsC( m_nevt, m_q )
    , m_hstSelHel( m_nevt, m_q )
    , m_hstSelCol( m_nevt, m_q )
    , m_devIsGoodHel( CPPPROCESS_NCOMB, m_q )
    , m_hstIsGoodHel( CPPPROCESS_NCOMB, m_q )
    #ifndef MGONGPU_HARDCODE_PARAM
    , m_dev_independent_couplings( Proc::independentCouplings::nicoup, m_q )
    , m_dev_independent_parameters( mgOnGpu::nparams, m_q )
    #endif
    , m_devcNGoodHel( 1, m_q ) 
    , m_devcGoodHel( CPPPROCESS_NCOMB, m_q ) 
  {
    if( nparF != CPPPROCESS_NPAR ) throw std::runtime_error( "Bridge constructor: npar mismatch" );
    if( np4F != CPPPROCESS_NP4 ) throw std::runtime_error( "Bridge constructor: np4 mismatch" );
    if( ( m_nevt < s_gputhreadsmin ) || ( m_nevt % s_gputhreadsmin != 0 ) )
      throw std::runtime_error( "Bridge constructor: nevt should be a multiple of " + std::to_string( s_gputhreadsmin ) );
    #if MGONGPU_VEC_DIM > 1
        while( m_nevt != m_gpublocks * m_gputhreads * MGONGPU_VEC_DIM) {
    #else
        while( m_nevt != m_gpublocks * m_gputhreads ) {
    #endif
      m_gputhreads /= 2;
      if( m_gputhreads < s_gputhreadsmin )
        throw std::logic_error( "Bridge constructor: FIXME! cannot choose gputhreads" ); // this should never happen!
      #if MGONGPU_VEC_DIM > 1
          m_gpublocks = m_nevt / m_gputhreads / MGONGPU_VEC_DIM;
      #else
          m_gpublocks = m_nevt / m_gputhreads;
      #endif
    }
    std::cout << "WARNING! Instantiate device Bridge (nevt=" << m_nevt << ", gpublocks=" << m_gpublocks << ", gputhreads=" << m_gputhreads
              << ", gpublocks*gputhreads=" << m_gpublocks * m_gputhreads << ")" << std::endl;

    #ifndef MGONGPU_HARDCODE_PARAM
    Proc::CPPProcess process( 1, m_gpublocks, m_gputhreads, false, false );
    process.initProc( "../../Cards/param_card.dat" );

    m_q.memcpy( m_dev_independent_couplings.data(), process.get_tIPC_ptr(), 2*Proc::independentCouplings::nicoup*sizeof(fptype) );
    m_q.memcpy( m_dev_independent_parameters.data(), process.get_tIPD_ptr(), mgOnGpu::nparams*sizeof(fptype) ).wait();
    #endif
  }

  template<typename FORTRANFPTYPE>
  void Bridge<FORTRANFPTYPE>::set_gpugrid( const size_t gpublocks, const size_t gputhreads )
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
                                            const FORTRANFPTYPE* __restrict__ rndhel,
                                            const FORTRANFPTYPE* __restrict__ rndcol,
                                            const size_t channelId,
                                            FORTRANFPTYPE* __restrict__ mes,
                                            int* __restrict__ selhel,
                                            int* __restrict__ selcol,
                                            const bool goodHelOnly ) {
    #if MGONGPU_VEC_DIM > 1
        host_buffer<fptype> hstMomentaC(CPPPROCESS_NPAR*MGONGPU_FOURVECTOR_DIM*m_nevt, m_q);
        hst_transposeMomenta(hstMomentaC.data(), momenta, m_nevt);
        m_q.memcpy(m_devMomentaC.data(), hstMomentaC.data(), CPPPROCESS_NPAR*MGONGPU_FOURVECTOR_DIM*m_nevt*sizeof(fptype));
    #else
        if constexpr (std::is_same_v<FORTRANFPTYPE, fptype>) {
            m_q.memcpy(m_devMomentaC.data(), momenta, CPPPROCESS_NPAR*MGONGPU_FOURVECTOR_DIM*m_nevt*sizeof(fptype));
        }
        else {
            host_buffer<fptype> hstMomentaC(CPPPROCESS_NPAR*MGONGPU_FOURVECTOR_DIM*m_nevt, m_q);
            std::copy(momenta, momenta + CPPPROCESS_NPAR*MGONGPU_FOURVECTOR_DIM*m_nevt, hstMomentaC.data());
            m_q.memcpy(m_devMomentaC.data(), hstMomentaC.data(), CPPPROCESS_NPAR*MGONGPU_FOURVECTOR_DIM*m_nevt*sizeof(fptype));
        }
    #endif

    if constexpr (std::is_same_v<FORTRANFPTYPE, fptype>) {
        m_q.memcpy(m_devGsC.data(), gs, m_nevt*sizeof(fptype));
        m_q.memcpy(m_devRndHel.data(), rndhel, m_nevt*sizeof(fptype));
        m_q.memcpy(m_devRndCol.data(), rndcol, m_nevt*sizeof(fptype)).wait();
    }
    else {
        host_buffer<fptype> hstGsC(m_nevt, m_q);
        host_buffer<fptype> hst_rndhel(m_nevt, m_q);
        host_buffer<fptype> hst_rndcol(m_nevt, m_q);
        std::copy(gs, gs + m_nevt, hstGsC.data());
        std::copy(rndhel, rndhel + m_nevt, hst_rndhel.data());
        std::copy(rndcol, rndcol + m_nevt, hst_rndcol.data());
        m_q.memcpy(m_devGsC.data(), hstGsC.data(), m_nevt*sizeof(fptype));
        m_q.memcpy(m_devRndHel.data(), hst_rndhel.data(), m_nevt*sizeof(fptype));
        m_q.memcpy(m_devRndCol.data(), hst_rndcol.data(), m_nevt*sizeof(fptype)).wait();
    }

    if (!m_goodHelsCalculated) {
        m_q.submit([&](sycl::handler& cgh) {
            auto devMomentaC = m_devMomentaC.data();
            auto devGsC = m_devGsC.data();
            auto devIsGoodHel = m_devIsGoodHel.data();
            #ifndef MGONGPU_HARDCODE_PARAM
                //Get pointer to independent couplings and parameters into shared memory if not hardcoded
                auto m_dev_independent_couplings_ptr = m_dev_independent_couplings.data();
                auto m_dev_parameters_ptr = m_dev_independent_parameters.data();
            #endif
            cgh.parallel_for_work_group(sycl::range<1>{m_gpublocks}, sycl::range<1>{m_gputhreads}, ([=](sycl::group<1> wGroup) {
                #ifndef MGONGPU_HARDCODE_PARAM
                    //Load independent couplings and parameters into shared memory if not hardcoded
                    auto dev_independent_couplings = m_dev_independent_couplings_ptr;
                #endif
                wGroup.parallel_for_work_item([&](sycl::h_item<1> index) {
                    size_t ievt = index.get_global_id(0);

                    #ifdef MGONGPU_HARDCODE_PARAM
                        //Load parameters into local (private) memory if hardcoded
                        auto dev_parameters = Proc::independent_parameters<fptype>;
                    #else
                        fptype dev_parameters[mgOnGpu::nparams];
                        for (size_t i = 0; i < mgOnGpu::nparams; i++) {
                            dev_parameters[i] = m_dev_parameters_ptr[i];
                        }
                    #endif

                    //Load helicities into local (private) memory
                    auto dev_helicities = Proc::helicities<signed char>;
                    cxtype_sv dev_couplings[Proc::dependentCouplings::ndcoup + Proc::independentCouplings::nicoup];

                    #if MGONGPU_NDCOUP > 0
                        Proc::dependentCouplings::set_couplings_from_G(dev_couplings, devGsC[ievt]); 
                    #endif

                    #if MGONGPU_NICOUP > 0
                        #ifdef MGONGPU_HARDCODE_PARAM
                        //Load independent couplings into local (private) memory if hardcoded
                        auto dev_independent_couplings = Proc::independentCouplings::independent_couplings<cxtype, fptype>;
                        #endif
                        for (size_t i = 0; i < Proc::independentCouplings::nicoup; i++) {
                            dev_couplings[Proc::dependentCouplings::ndcoup + i] = dev_independent_couplings[i];
                        }
                    #endif

                    //Load momenta into local (private) memory
                    vector4 p_momenta[CPPPROCESS_NPAR];
                    for (size_t i = 0; i < CPPPROCESS_NPAR; i++) {
                        p_momenta[i] = devMomentaC[CPPPROCESS_NPAR*ievt + i];
                    }

                    Proc::sigmaKin_getGoodHel( p_momenta, devIsGoodHel, dev_helicities, dev_couplings, dev_parameters );
                });
            }));
        });
        m_q.wait();

        m_q.memcpy(m_hstIsGoodHel.data(), m_devIsGoodHel.data(), CPPPROCESS_NCOMB*sizeof(bool)).wait();

        size_t goodHel[CPPPROCESS_NCOMB] = {0};
        size_t nGoodHel = Proc::sigmaKin_setGoodHel( m_hstIsGoodHel.data(), goodHel );
        m_nGoodHel = int(nGoodHel);

        m_q.memcpy( m_devcNGoodHel.data(), &nGoodHel, sizeof(size_t) ).wait();
        m_q.memcpy( m_devcGoodHel.data(), goodHel, CPPPROCESS_NCOMB*sizeof(size_t) ).wait();
        m_goodHelsCalculated = true;
    }
    if( goodHelOnly ) return;

    m_q.submit([&](sycl::handler& cgh) {
        auto devMomentaC = m_devMomentaC.data();
        auto devGsC = m_devGsC.data();
        auto devRndHel = m_devRndHel.data();
        auto devRndCol = m_devRndCol.data();
        auto devcNGoodHel = m_devcNGoodHel.data();
        auto devcGoodHel = m_devcGoodHel.data();
        auto devMEsC = m_devMEsC.data();
        auto devSelHel = m_devSelHel.data();
        auto devSelCol = m_devSelCol.data();
        #ifndef MGONGPU_HARDCODE_PARAM
            //Get pointer to independent couplings and parameters into shared memory if not hardcoded
            auto m_dev_independent_couplings_ptr = m_dev_independent_couplings.data();
            auto m_dev_parameters_ptr = m_dev_independent_parameters.data();
        #endif
        cgh.parallel_for_work_group(sycl::range<1>{m_gpublocks}, sycl::range<1>{m_gputhreads}, ([=](sycl::group<1> wGroup) {
            #ifndef MGONGPU_HARDCODE_PARAM
                //Load independent couplings and parameters into shared memory if not hardcoded
                auto dev_independent_couplings = m_dev_independent_couplings_ptr;
            #endif
            auto l_channelId = channelId;
            wGroup.parallel_for_work_item([&](sycl::h_item<1> index) {
                size_t ievt = index.get_global_id(0);

                #ifdef MGONGPU_HARDCODE_PARAM
                    //Load parameters into local (private) memory if hardcoded
                    auto dev_parameters = Proc::independent_parameters<fptype>;
                #else
                    fptype dev_parameters[mgOnGpu::nparams];
                    for (size_t i = 0; i < mgOnGpu::nparams; i++) {
                        dev_parameters[i] = m_dev_parameters_ptr[i];
                    }
                #endif

                //Load helicities and couplings into local (private) memory
                auto dev_helicities = Proc::helicities<signed char>;
                cxtype_sv dev_couplings[Proc::dependentCouplings::ndcoup + Proc::independentCouplings::nicoup];

                #if MGONGPU_NDCOUP > 0
                    Proc::dependentCouplings::set_couplings_from_G(dev_couplings, devGsC[ievt]); 
                #endif

                #if MGONGPU_NICOUP > 0
                    #ifdef MGONGPU_HARDCODE_PARAM
                        //Load independent couplings into local (private) memory if hardcoded
                        auto dev_independent_couplings = Proc::independentCouplings::independent_couplings<cxtype, fptype>;
                    #endif
                    for (size_t i = 0; i < Proc::independentCouplings::nicoup; i++) {
                        dev_couplings[Proc::dependentCouplings::ndcoup + i] = dev_independent_couplings[i];
                    }
                #endif

                //Load momenta into local (private) memory
                vector4 p_momenta[CPPPROCESS_NPAR];
                for (size_t i = 0; i < CPPPROCESS_NPAR; i++) {
                    p_momenta[i] = devMomentaC[CPPPROCESS_NPAR*ievt + i];
                }

                #ifdef MGONGPU_SUPPORTS_MULTICHANNEL
                    devMEsC[ievt] = Proc::sigmaKin( p_momenta, devRndHel + ievt, devRndCol + ievt, devSelHel + ievt, devSelCol + ievt, l_channelId, dev_helicities, dev_couplings, dev_parameters, devcNGoodHel, devcGoodHel );
                #else
                    devMEsC[ievt] = Proc::sigmaKin( p_momenta, devRndHel + ievt, devRndCol + ievt, devSelHel + ievt, devSelCol + ievt, dev_helicities, dev_couplings, dev_parameters, devcNGoodHel, devcGoodHel );
                #endif
            });
        }));
    });
    m_q.wait();

    m_q.memcpy(m_hstSelHel.data(), m_devSelHel.data(), m_nevt*sizeof(long));
    m_q.memcpy(m_hstSelCol.data(), m_devSelCol.data(), m_nevt*sizeof(long));
    if constexpr(std::is_same_v<FORTRANFPTYPE, fptype>) {
        m_q.memcpy(mes, m_devMEsC.data(), m_nevt*sizeof(fptype)).wait();
    }
    else {
        m_q.memcpy( m_hstMEsC.data(), m_devMEsC.data(), m_nevt*sizeof(fptype)).wait();
        std::copy( m_hstMEsC.data(), m_hstMEsC.data() + m_nevt, mes);
    }
    std::copy( m_hstSelHel.data(), m_hstSelHel.data() + m_nevt, selhel);
    std::copy( m_hstSelCol.data(), m_hstSelCol.data() + m_nevt, selcol);

  }

  template<typename FPTypeDst, typename FPTypeSrc>
  void hst_transposeMomenta(FPTypeDst* __restrict__ dst, const FPTypeSrc* __restrict__ src, const size_t N ) {
  /* Transpose from [ evt0par0wxyz, evt0par1wxyz, ... , evt0parNPARwxyz, evt1par0wxyz, ... ] to
     [ evt0par0w, evt1par0w, ... , evt(MGONGPU_VEC_DIM - 1)par0w,
       evt0par0x, evt1par0x, ... , evt(MGONGPU_VEC_DIM - 1)par0x,
       evt0par0y, evt1par0y, ... , evt(MGONGPU_VEC_DIM - 1)par0y,
       evt0par0z, evt1par0z, ... , evt(MGONGPU_VEC_DIM - 1)par0z,
       evt0par1w, evt1par1w, ... , evt(MGONGPU_VEC_DIM - 1)par1w,
       ... ,
       evt0parNPARz, evt1parNPARz, ... , evtMGONGPU_VEC_DIMparNPARz,
       evt(MGONGPU_VEC_DIM)par0w, evt(MGONGPU_VEC_DIM + 1)par0w, ... , evt(2*MGONGPU_VEC_DIM - 1)par0w,
       ... ,
       evt(N - MGONGPU_VEC_DIM)parNPARz, evt(N - MGONGPU_VEC_DIM + 1)parNPARz, ... , evtNparNPARz
     ]
  */
      for ( size_t h = 0; h < N/MGONGPU_VEC_DIM; h++ ) {
          for ( size_t i = 0; i < MGONGPU_VEC_DIM; i++ ) {
              size_t idx_evt = h*MGONGPU_VEC_DIM + i;
              for ( size_t j = 0; j < CPPPROCESS_NPAR; j++ ) { // size_t idx_par = j
                  for ( size_t k = 0; k < MGONGPU_FOURVECTOR_DIM; k++ ) { // size_t idx_wxyz = k
                      size_t idx_src = CPPPROCESS_NPAR*MGONGPU_FOURVECTOR_DIM*idx_evt + MGONGPU_FOURVECTOR_DIM*j + k;
                      size_t idx_dst = CPPPROCESS_NPAR*MGONGPU_FOURVECTOR_DIM*MGONGPU_VEC_DIM*h + MGONGPU_FOURVECTOR_DIM*MGONGPU_VEC_DIM*j + MGONGPU_VEC_DIM*k + i;
                      dst[idx_dst] = src[idx_src];
                  }
              }
          }
      }
  }

  //--------------------------------------------------------------------------
}
#endif // BRIDGE_H
