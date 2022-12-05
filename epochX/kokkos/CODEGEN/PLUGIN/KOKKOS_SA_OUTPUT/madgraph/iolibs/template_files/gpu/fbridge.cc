#include "Bridge.h"

extern "C"
{
  /**
   * The namespace where the Bridge class is taken from.
   *
   * In the current implementation, shared libraries are created for the SYCL device implementation.
   * A single fcreatebridge_ symbol is created in each library with the same name, connected to the appropriate Bridge on SYCL devices.
   * The Fortran MadEvent code is always the same: the choice whether to use a specific SYCL device implementation is done by linking the appropriate library and setting an evironment variable.
   */
  using namespace mg5amcKokkos;

  /**
   * The floating point precision used in Fortran arrays.
   * This is presently hardcoded to double precision (REAL*8).
   */
  using FORTRANFPTYPE = double; // for Fortran double precision (REAL*8) arrays
  //using FORTRANFPTYPE = float; // for Fortran single precision (REAL*4) arrays

  /**
   * Create a Bridge and return its pointer.
   * This is a C symbol that should be called from the Fortran code (in auto_dsig1.f).
   *
   * @param ppbridge the pointer to the Bridge pointer (the Bridge pointer is handled in Fortran as an INTEGER*8 variable)
   * @param nevtF the pointer to the number of events in the Fortran arrays
   * @param nparF the pointer to the number of external particles in the Fortran arrays (KEPT FOR SANITY CHECKS ONLY)
   * @param np4F the pointer to the number of momenta components, usually 4, in the Fortran arrays (KEPT FOR SANITY CHECKS ONLY)
   */
  void fbridgecreate_( CppObjectInFortran** ppbridge, const int* pnevtF, const int* pnparF, const int* pnp4F )
  {
    bool _kokkos_initialized = kokkos_initialize();
    *ppbridge = new Bridge<FORTRANFPTYPE>( *pnevtF, *pnparF, *pnp4F );
  }

  /**
   * Delete a Bridge.
   * This is a C symbol that should be called from the Fortran code (in auto_dsig1.f).
   *
   * @param ppbridge the pointer to the Bridge pointer (the Bridge pointer is handled in Fortran as an INTEGER*8 variable)
   */
  void fbridgedelete_( CppObjectInFortran** ppbridge )
  {
    Bridge<FORTRANFPTYPE>* pbridge = dynamic_cast<Bridge<FORTRANFPTYPE>*>( *ppbridge );
    if( pbridge == 0 ) throw std::runtime_error( "fbridgedelete_: invalid Bridge address" );
    delete pbridge;
    Kokkos::finalize();
  }

  /**
   * Execute the matrix-element calculation "sequence" via a Bridge on GPU/CUDA or CUDA/C++.
   * This is a C symbol that should be called from the Fortran code (in auto_dsig1.f).
   *
   * @param ppbridge the pointer to the Bridge pointer (the Bridge pointer is handled in Fortran as an INTEGER*8 variable)
   * @param momenta the pointer to the input 4-momenta
   * @param gs the pointer to the input Gs (running QCD coupling constant alphas)
   * @param mes the pointer to the output matrix elements
   * @param channelId the pointer to the Feynman diagram to enhance in multi-channel mode if 1 to n (disable multi-channel if 0)
   */
  void fbridgesequence_( CppObjectInFortran** ppbridge,
                         const FORTRANFPTYPE* momenta,
                         const FORTRANFPTYPE* gs,
                         FORTRANFPTYPE* mes,
                         const unsigned int* pchannelId )
  {
    Bridge<FORTRANFPTYPE>* pbridge = dynamic_cast<Bridge<FORTRANFPTYPE>*>( *ppbridge );
    if( pbridge == 0 ) throw std::runtime_error( "fbridgesequence_: invalid Bridge address" );
    pbridge->gpu_sequence( momenta, gs, mes, *pchannelId );
  }
}
