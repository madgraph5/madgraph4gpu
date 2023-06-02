// Copyright (C) 2020-2023 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Feb 2022) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2022-2023) for the MG5aMC CUDACPP plugin.
//==========================================================================
// Copyright (C) 2021-2023 Argonne National Laboratory.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Modified by: N. Nichols (2021-2023) for the MG5aMC SYCL plugin.
//==========================================================================

#include "Bridge.h"
#include "Memory.h"
#include "rambo.h"
#include "CommonRandomNumbers.h"

//--------------------------------------------------------------------------

namespace mg5amcGpu
{
  template<typename FORTRANFPTYPE>
  class Sampler final : public CppObjectInFortran
  {
  public:
    // Constructor
    // @param nevtF (NB_PAGE_LOOP, vector.inc) number of events in Fortran arrays
    // @param nparF (NEXTERNAL, nexternal.inc) number of external particles in Fortran arrays (KEPT FOR SANITY CHECKS ONLY: remove it?)
    // @param np4F number of momenta components, usually 4, in Fortran arrays (KEPT FOR SANITY CHECKS ONLY: remove it?)
    Sampler( int nevtF, int nparF, int np4F );
    // Destructor
    virtual ~Sampler() {}
    // Delete copy/move constructors and assignment operators
    Sampler( const Sampler& ) = delete;
    Sampler( Sampler&& ) = delete;
    Sampler& operator=( const Sampler& ) = delete;
    Sampler& operator=( Sampler&& ) = delete;
    // Draw random numbers and convert them to momenta in C++, then transpose them to Fortran momenta
    void samplerHostSequence( FORTRANFPTYPE* fortranMomenta );
  private:
    const int m_nevt; // The number of events in each iteration
    int m_iiter;      // The iteration counter (for random number seeding)

    host_buffer_unique<fptype> m_hstRnarray; // Memory buffers for random numbers
    host_buffer_unique<fptype> m_hstMomenta; // Memory buffers for momenta
    host_buffer_unique<fptype> m_hstWeights; // Memory buffers for sampling weights

    // HARDCODED DEFAULTS
    static constexpr fptype energy = 1500; // historical default, Ecms = 1500 GeV = 1.5 TeV (above the Z peak)
  };

  template<typename FORTRANFPTYPE>
  Sampler<FORTRANFPTYPE>::Sampler( int nevtF, int nparF, int np4F )
    : m_nevt( nevtF )
    , m_iiter( 0 )
    , m_hstRnarray( CPPPROCESS_NP4*CPPPROCESS_NPARF*m_nevt )
    , m_hstMomenta( CPPPROCESS_NP4*CPPPROCESS_NPAR*m_nevt )
    , m_hstWeights( m_nevt )
  {
    if( nparF != CPPPROCESS_NPAR ) throw std::runtime_error( "Sampler constructor: npar mismatch" );
    if( np4F != CPPPROCESS_NP4 ) throw std::runtime_error( "Sampler constructor: np4 mismatch" );
    std::cout << "WARNING! Instantiate host Sampler (nevt=" << m_nevt << ")" << std::endl;
  }

  // Draw random numbers and convert them to momenta in C++, then transpose them to Fortran momenta
  template<typename FORTRANFPTYPE>
  void Sampler<FORTRANFPTYPE>::samplerHostSequence( FORTRANFPTYPE* fortranMomenta )
  {
    std::cout << "Iteration #" << m_iiter + 1 << std::endl;
    // === STEP 1 OF 3
    // --- 1a. Seed rnd generator (to get same results on host and device in curand)
    // [NB This should not be necessary using the host API: "Generation functions
    // can be called multiple times on the same generator to generate successive
    // blocks of results. For pseudorandom generators, multiple calls to generation
    // functions will yield the same result as a single call with a large size."]
    // *** NB! REMEMBER THAT THE FORTRAN SAMPLER ALWAYS USES COMMON RANDOM NUMBERS! ***
    constexpr unsigned long long seed = 20200805;
    std::vector<double> rnd = CommonRandomNumbers::generate<double>( m_hstRnarray.size(), seed + m_iiter );
    std::copy( rnd.begin(), rnd.end(), m_hstRnarray.data() );
    m_iiter++;

    // ** START LOOP ON IEVT **
    for( size_t ievt = 0; ievt < m_nevt; ++ievt )
    {
      rambo2toNm0::ramboGetMomentaInitial( energy, m_hstMomenta.data(), ievt );
      rambo2toNm0::ramboGetMomentaFinal( energy, m_hstRnarray.data(), m_hstMomenta.data(), m_hstWeights.data(), ievt );
    }
    hst_transposeMomentaC2F( m_hstMomenta.data(), fortranMomenta, m_nevt );
  }
}

//--------------------------------------------------------------------------

extern "C"
{
  using namespace mg5amcGpu;

  /**
   * The floating point precision used in Fortran arrays.
   * This is presently hardcoded to double precision (REAL*8).
   */
  using FORTRANFPTYPE = double; // for Fortran double precision (REAL*8) arrays
  //using FORTRANFPTYPE = float; // for Fortran single precision (REAL*4) arrays

  /**
   * Create a Sampler and return its pointer.
   * This is a C symbol that should be called from the Fortran code (in auto_dsig1.f).
   *
   * @param ppsampler the pointer to the Sampler pointer (the Sampler pointer is handled in Fortran as an INTEGER*8 variable)
   * @param nevtF the pointer to the number of events in the Fortran arrays
   * @param nparF the pointer to the number of external particles in the Fortran arrays (KEPT FOR SANITY CHECKS ONLY)
   * @param np4F the pointer to the number of momenta components, usually 4, in the Fortran arrays (KEPT FOR SANITY CHECKS ONLY)
   */
  void fsamplercreate_( CppObjectInFortran** ppsampler, const int* pnevtF, const int* pnparF, const int* pnp4F )
  {
    *ppsampler = new Sampler<FORTRANFPTYPE>( *pnevtF, *pnparF, *pnp4F );
  }

  /**
   * Delete a Sampler.
   * This is a C symbol that should be called from the Fortran code (in auto_dsig1.f).
   *
   * @param ppsampler the pointer to the Sampler pointer (the Sampler pointer is handled in Fortran as an INTEGER*8 variable)
   */
  void fsamplerdelete_( CppObjectInFortran** ppsampler )
  {
    Sampler<FORTRANFPTYPE>* psampler = dynamic_cast<Sampler<FORTRANFPTYPE>*>( *ppsampler );
    if( psampler == 0 ) throw std::runtime_error( "fsamplerdelete_: invalid Sampler address" );
    delete psampler;
  }

  /**
   * Execute the matrix-element calculation "sequence" via a Sampler on SYCL devices.
   * This is a C symbol that should be called from the Fortran code (in auto_dsig1.f).
   *
   * @param ppsampler the pointer to the Sampler pointer (the Sampler pointer is handled in Fortran as an INTEGER*8 variable)
   * @param momenta the pointer to the input 4-momenta
   * @param mes the pointer to the output matrix elements
   */
  void fsamplersequence_( CppObjectInFortran** ppsampler, FORTRANFPTYPE* momenta )
  {
    Sampler<FORTRANFPTYPE>* psampler = dynamic_cast<Sampler<FORTRANFPTYPE>*>( *ppsampler );
    if( psampler == 0 ) throw std::runtime_error( "fsamplersequence_: invalid Sampler address" );
    // Use the host/CPU implementation (there is no device implementation)
    psampler->samplerHostSequence( momenta );
  }
}

//--------------------------------------------------------------------------
