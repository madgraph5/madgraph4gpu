// Copyright (C) 2020-2023 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: S. Roiser (Oct 2021) for the MG5aMC CUDACPP plugin.
// Further modified by: S. Roiser, A. Valassi (2021-2023) for the MG5aMC CUDACPP plugin.
//
// Copyright (C) 2021-2023 Argonne National Laboratory.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Modified by: N. Nichols (2021-2023) for the MG5aMC SYCL plugin.

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
  using namespace mg5amcGpu;

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
  }

  /**
   * Execute the matrix-element calculation "sequence" via a Bridge on SYCL device.
   * This is a C symbol that should be called from the Fortran code (in auto_dsig1.f).
   *
   * @param ppbridge the pointer to the Bridge pointer (the Bridge pointer is handled in Fortran as an INTEGER*8 variable)
   * @param momenta the pointer to the input 4-momenta
   * @param gs the pointer to the input Gs (running QCD coupling constant alphas)
   * @param rndhel the pointer to the input random numbers for helicity selection
   * @param rndcol the pointer to the input random numbers for color selection
   * @param channelId the pointer to the input Feynman diagram to enhance in multi-channel mode if 1 to n (disable multi-channel if 0)
   * @param mes the pointer to the output matrix elements
   * @param selhel the pointer to the output selected helicities
   * @param selcol the pointer to the output selected colors
   */
  void fbridgesequence_( CppObjectInFortran** ppbridge,
                         const FORTRANFPTYPE* momenta,
                         const FORTRANFPTYPE* gs,
                         const FORTRANFPTYPE* rndhel,
                         const FORTRANFPTYPE* rndcol,
                         const unsigned int* pchannelId,
                         FORTRANFPTYPE* mes,
                         int* selhel,
                         int* selcol )
  {
    Bridge<FORTRANFPTYPE>* pbridge = dynamic_cast<Bridge<FORTRANFPTYPE>*>( *ppbridge );
    if( pbridge == 0 ) throw std::runtime_error( "fbridgesequence_: invalid Bridge address" );
    pbridge->gpu_sequence( momenta, gs, rndhel, rndcol, *pchannelId, mes, selhel, selcol );
  }

  /**
   * Retrieve the number of good helicities for helicity filtering in the Bridge.
   * This is a C symbol that should be called from the Fortran code (in auto_dsig1.f).
   *
   * @param ppbridge the pointer to the Bridge pointer (the Bridge pointer is handled in Fortran as an INTEGER*8 variable)
   * @param pngoodhel the pointer to the output number of good helicities
   * @param pntothel the pointer to the output total number of helicities
   */
  void fbridgegetngoodhel_( CppObjectInFortran** ppbridge,
                            unsigned int* pngoodhel,
                            unsigned int* pntothel )
  {
    Bridge<FORTRANFPTYPE>* pbridge = dynamic_cast<Bridge<FORTRANFPTYPE>*>( *ppbridge );
    if( pbridge == 0 ) throw std::runtime_error( "fbridgegetngoodhel_: invalid Bridge address" );
    *pngoodhel = pbridge->nGoodHel();
    *pntothel = pbridge->nTotHel();
  }
}
