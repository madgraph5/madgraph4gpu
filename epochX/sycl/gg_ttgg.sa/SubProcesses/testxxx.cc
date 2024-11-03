// Copyright (C) 2020-2023 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Apr 2021) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2021-2023) for the MG5aMC CUDACPP plugin.
//
// Copyright (C) 2021-2023 Argonne National Laboratory.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Modified by: N. Nichols (2021-2023) for the MG5aMC SYCL plugin.
#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"
#include "mgOnGpuVectors.h"
#include "HelAmps_sm.h"

#include "CPPProcess.h"
#include "Memory.h"

#include <array>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include <gtest/gtest.h>
#include "epoch_process_id.h"
#define TESTID_CPU( s ) s##_CPU
#define XTESTID_CPU( s ) TESTID_CPU( s )

TEST( XTESTID_CPU( MG_EPOCH_PROCESS_ID ), testxxx )
{
  constexpr bool dumpEvents = false;  // dump the expected output of the test?
  constexpr bool testEvents = !dumpEvents; // run the test?
  constexpr fptype toleranceXXXs = std::is_same<fptype, double>::value ? 1.E-15 : 1.E-5;
  // Constant parameters
  using mgOnGpu::neppM;
  const int nevt = 16; // 12 independent tests plus 4 duplicates (need a multiple of 8 for floats or for '512z')
  assert( nevt % neppM == 0 ); // nevt must be a multiple of neppM
  // Fill in the input momenta
  const int nMomenta = CPPPROCESS_NP4*CPPPROCESS_NPAR*nevt;
  auto hstMomenta = hstMakeUnique<fptype_sv>( nMomenta ); // AOSOA[npagM][CPPPROCESS_NPAR=4][CPPPROCESS_NP4=4][neppM]
  const fptype par0[CPPPROCESS_NP4* nevt] =                         // AOS[nevt][CPPPROCESS_NP4]
    { 500, 0,    0,    500,  // #0 (m=0 pT=0 E=pz>0)
      500, 0,    0,    -500, // #1 (m=0 pT=0 -E=pz<0)
      500, 300,  400,  0,    // #2 (m=0 pT>0 pz=0)
      500, 180,  240,  400,  // #3 (m=0 pT>0 pz>0)
      500, 180,  240,  -400, // #4 (m=0 pT>0 pz<0)
      500, 0,    0,    0,    // #5 (m=50>0 pT=0 pz=0)
      500, 0,    0,    300,  // #6 (m=40>0 pT=0 pz>0)
      500, 0,    0,    -300, // #7 (m=40>0 pT=0 pz<0)
      500, 180,  240,  0,    // #8 (m=40>0 pT>0 pz=0)
      500, -240, -180, 0,    // #9 (m=40>0 pT>0 pz=0)
      500, 180,  192,  144,  // #10 (m=40>0 pT>0 pz>0)
      500, 180,  192,  -144, // #11 (m=40>0 pT>0 pz<0)
      500, 0,    0,    500,  // DUPLICATE #12 == #0 (m=0 pT=0 E=pz>0)
      500, 0,    0,    -500, // DUPLICATE #13 == #1 (m=0 pT=0 -E=pz<0)
      500, 300,  400,  0,    // DUPLICATE #14 == #2 (m=0 pT>0 pz=0)
      500, 180,  240,  400   // DUPLICATE #15 == #3 (m=0 pT>0 pz>0)
    };
  fptype mass0[nevt]{};
  bool ispzgt0[nevt]{};
  bool ispzlt0[nevt]{};
  bool isptgt0[nevt]{};
  for ( int ievt = 0; ievt < nevt; ievt++ )
  {
    const fptype p0 = par0[ievt * CPPPROCESS_NP4 + 0];
    const fptype p1 = par0[ievt * CPPPROCESS_NP4 + 1];
    const fptype p2 = par0[ievt * CPPPROCESS_NP4 + 2];
    const fptype p3 = par0[ievt * CPPPROCESS_NP4 + 3];
    mass0[ievt] = sqrt( p0 * p0 - p1 * p1 - p2 * p2 - p3 * p3 );
    ispzgt0[ievt] = ( p3 > 0 );
    ispzlt0[ievt] = ( p3 < 0 );
    isptgt0[ievt] = ( p1 != 0 ) || ( p2 != 0 );
  }
  const int ipar = 0; // use only particle0 for this test
  for ( int ievt = 0; ievt < nevt; ievt++ )
  {
    for ( int ip4 = 0; ip4 < CPPPROCESS_NP4; ip4++ )
    {
      const int ipagM = ievt/neppM; // #eventpage in this iteration
      const int ieppM = ievt%neppM; // #event in the current eventpage in this iteration
#ifdef MGONGPU_CPPSIMD
      hstMomenta[ipagM*CPPPROCESS_NPAR*CPPPROCESS_NP4 + ipar*CPPPROCESS_NP4 + ip4][ieppM] = par0[ievt*CPPPROCESS_NP4 + ip4]; // AOS to AOSOA
#else
      hstMomenta[ipagM*CPPPROCESS_NPAR*CPPPROCESS_NP4*neppM + ipar*CPPPROCESS_NP4*neppM + ip4*neppM + ieppM] = par0[ievt*CPPPROCESS_NP4 + ip4]; // AOS to AOSOA
#endif
    }
  }
  // Expected output wavefunctions
  std::vector<std::array<fptype, 12>> expwfs;
#include "testxxx_cc_ref.txt" // expwfs.push_back( {...} );
  std::string dumpFileName = "testxxx_cc_ref.txt.new";
  // Compute the output wavefunctions
  // Dump new reference file if requested
  using namespace MG5_sm;
  int itest = 0; // index on the expected output vector
  std::ofstream dumpFile;
  if ( dumpEvents ) dumpFile.open( dumpFileName, std::ios::trunc );
  auto dumpwf6 = [&]( std::ostream& out, const cxtype_sv wf[6], const char* xxx, int ievt, int nsp, fptype mass ) {
    out << std::setprecision(15) << std::scientific;
    out << "  expwfs.push_back( {";
    out << "                                   // ---------" << std::endl;
    for ( int iw6 = 0; iw6<CPPPROCESS_NW6; iw6++ )
    {
#ifdef MGONGPU_CPPSIMD
      const int ieppM = ievt%neppM; // #event in the current eventpage in this iteration
#ifdef MGONGPU_HAS_CXTYPE_REF
      out << std::setw(26) << cxreal( wf[iw6][ieppM] ) << ", ";
      out << std::setw(22) << cximag( wf[iw6][ieppM] );
#else
      out << std::setw(26) << wf[iw6].real()[ieppM] << ", ";
      out << std::setw(22) << wf[iw6].imag()[ieppM];
#endif
#else
      out << std::setw(26) << wf[iw6].real();
      out << ", " << std::setw(22) << wf[iw6].imag();
#endif
      if ( iw6 < CPPPROCESS_NW6 - 1 ) out << ",    ";
      else out << " } );";
      out << " // itest=" << itest << ": " << xxx << "#" << ievt;
      out << " nsp=" << nsp << " mass=" << (int)mass << std::endl;
    }
    out << std::defaultfloat;
  };
  auto testwf6 = [&]( const cxtype_sv wf[6], const char* xxx, int ievt, int nsp, fptype mass ) {
    if ( dumpEvents ) dumpwf6( dumpFile, wf, xxx, ievt, nsp, mass );
    if ( testEvents )
    {
      //std::cout << "Testing " << std::setw(3) << itest << ": " << xxx << " #" << ievt << std::endl;
      std::array<fptype, 12>& expwf = expwfs[itest];
      for ( int iw6 = 0; iw6 < CPPPROCESS_NW6; iw6++ )
      {
        const fptype expReal = expwf[iw6*2];
        const fptype expImag = expwf[iw6*2+1];
        if ( true )
        {
#ifdef MGONGPU_CPPSIMD
          const int ieppM = ievt%neppM; // #event in the current eventpage in this iteration
#ifdef MGONGPU_HAS_CXTYPE_REF
          EXPECT_NEAR( cxreal( wf[iw6][ieppM] ), expReal, std::abs( expReal * toleranceXXXs ) )
            << " itest=" << itest << ": " << xxx << "#" << ievt;
          EXPECT_NEAR( cximag( wf[iw6][ieppM] ), expImag, std::abs( expImag * toleranceXXXs ) )
            << " itest=" << itest << ": " << xxx << "#" << ievt;
#else
          EXPECT_NEAR( wf[iw6].real()[ieppM], expReal, std::abs( expReal * toleranceXXXs ) )
            << " itest=" << itest << ": " << xxx << "#" << ievt;
          EXPECT_NEAR( wf[iw6].imag()[ieppM], expImag, std::abs( expImag * toleranceXXXs ) )
            << " itest=" << itest << ": " << xxx << "#" << ievt;
#endif
#else
          EXPECT_NEAR( cxreal( wf[iw6] ), expReal, std::abs( expReal * toleranceXXXs ) )
            << " itest=" << itest << ": " << xxx << "#" << ievt;
          EXPECT_NEAR( cximag( wf[iw6] ), expImag, std::abs( expImag * toleranceXXXs ) )
            << " itest=" << itest << ": " << xxx << "#" << ievt;
#endif
        }
      }
    }
    itest++;
  };
  auto testwf6two = [&]( const cxtype_sv wf[6], const cxtype_sv expwf[6], const char* xxx, int ievt ) {
    if ( testEvents )
    {
      const std::string xxxFull( xxx[0] == 'i' ? "ixxxxx" : "oxxxxx" );
      //std::cout << "Testing " << std::setw(3) << itest << ": ";
      //std::cout << xxx << " #" << ievt << " against " << xxxFull << std::endl;
      ////for ( int iw6 = 0; iw6 < CPPPROCESS_NW6; iw6++ ) std::cout << wf[iw6] << std::endl;
      ////std::cout << "against" << std::endl;
      ////for ( int iw6 = 0; iw6 < CPPPROCESS_NW6; iw6++ ) std::cout << expwf[iw6] << std::endl;
      for ( int iw6 = 0; iw6 < CPPPROCESS_NW6; iw6++ )
      {
        if ( true )
        {
#ifdef MGONGPU_CPPSIMD
          const int ieppM = ievt%neppM; // #event in the current eventpage in this iteration
#ifdef MGONGPU_HAS_CXTYPE_REF
          const fptype expReal = cxreal( expwf[iw6][ieppM] );
          const fptype expImag = cximag( expwf[iw6][ieppM] );
          EXPECT_NEAR( cxreal( wf[iw6][ieppM] ), expReal, std::abs( expReal * toleranceXXXs ) )
            << " itest=" << itest << ": " << xxx << "#" << ievt << " against " << xxxFull;
          EXPECT_NEAR( cximag( wf[iw6][ieppM] ), expImag, std::abs( expImag * toleranceXXXs ) )
            << " itest=" << itest << ": " << xxx << "#" << ievt << " against " << xxxFull;
#else
          const fptype expReal = expwf[iw6].real()[ieppM];
          const fptype expImag = expwf[iw6].imag()[ieppM];
          EXPECT_NEAR( wf[iw6].real()[ieppM], expReal, std::abs( expReal * toleranceXXXs ) )
            << " itest=" << itest << ": " << xxx << "#" << ievt << " against " << xxxFull;
          EXPECT_NEAR( wf[iw6].imag()[ieppM], expImag, std::abs( expImag * toleranceXXXs ) )
            << " itest=" << itest << ": " << xxx << "#" << ievt << " against " << xxxFull;
#endif
#else
          const fptype expReal = cxreal( expwf[iw6] );
          const fptype expImag = cximag( expwf[iw6] );
          EXPECT_NEAR( cxreal( wf[iw6] ), expReal, std::abs( expReal * toleranceXXXs ) )
            << " itest=" << itest << ": " << xxx << "#" << ievt << " against " << xxxFull;
          EXPECT_NEAR( cximag( wf[iw6] ), expImag, std::abs( expImag * toleranceXXXs ) )
            << " itest=" << itest << ": " << xxx << "#" << ievt << " against " << xxxFull;
#endif
        }
      }
    }
  };
  const int nhel = 1;
  cxtype_sv outwfI[6]{}; // last result of ixxxxx (mass==0)
  cxtype_sv outwfO[6]{}; // last result of oxxxxx (mass==0)
  cxtype_sv outwf[6]{};
  for ( auto nsp : { -1, +1 } ) // antifermion/fermion (or initial/final for scalar and vector)
  {
    for ( int ievt = 0; ievt < nevt; ievt++ )
    {
      if ( false )
      {
        std::cout << std::endl;
        for ( int ip4 = 0; ip4 < CPPPROCESS_NP4; ip4++ ) std::cout << par0[ievt * CPPPROCESS_NP4 + ip4] << ", ";
        std::cout << std::endl;
      }
      // Test ixxxxx - NO ASSUMPTIONS
      {
        const fptype fmass = mass0[ievt];
        const int ipagM = ievt/neppM; // #eventpage in this iteration
        ixxxxx( hstMomenta.get(), fmass, nhel, nsp, outwfI, ipagM, ipar );
        testwf6( outwfI, "ixxxxx", ievt, nsp, fmass );
        ixxxxx( hstMomenta.get(), -fmass, nhel, nsp, outwfI, ipagM, ipar );
        testwf6( outwfI, "ixxxxx", ievt, nsp, -fmass );
      }
      // Test ipzxxx - ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
      if ( mass0[ievt] == 0 && !isptgt0[ievt] && ispzgt0[ievt] )
      {
        const int ipagM = ievt/neppM; // #eventpage in this iteration
        ipzxxx( hstMomenta.get(), nhel, nsp, outwf, ipagM, ipar );
        testwf6two( outwf, outwfI, "ipzxxx", ievt );
        testwf6( outwf, "ipzxxx", ievt, nsp, 0 );
      }
      // Test imzxxx - ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
      if ( mass0[ievt] == 0 && !isptgt0[ievt] && ispzlt0[ievt] )
      {
        const int ipagM = ievt/neppM; // #eventpage in this iteration
        imzxxx( hstMomenta.get(), nhel, nsp, outwf, ipagM, ipar );
        testwf6two( outwf, outwfI, "imzxxx", ievt );
        testwf6( outwf, "imzxxx", ievt, nsp, 0 );
      }
      // Test ixzxxx - ASSUMPTIONS: (FMASS == 0) and (PT > 0)
      if ( mass0[ievt] == 0 && isptgt0[ievt] )
      {
        const int ipagM = ievt/neppM; // #eventpage in this iteration
        ixzxxx( hstMomenta.get(), nhel, nsp, outwf, ipagM, ipar );
        testwf6two( outwf, outwfI, "ixzxxx", ievt );
        testwf6( outwf, "ixzxxx", ievt, nsp, 0 );
      }
      // Test vxxxxx - NO ASSUMPTIONS
      {
        const fptype vmass = mass0[ievt];
        const int ipagM = ievt/neppM; // #eventpage in this iteration
        vxxxxx( hstMomenta.get(), vmass, nhel, nsp, outwf, ipagM, ipar );
        testwf6( outwf, "vxxxxx", ievt, nsp, vmass );
        vxxxxx( hstMomenta.get(), -vmass, nhel, nsp, outwf, ipagM, ipar );
        testwf6( outwf, "vxxxxx", ievt, nsp, -vmass );
      }
      // Test sxxxxx - NO ASSUMPTIONS
      {
        const fptype smass = mass0[ievt];
        cxtype_sv outwf3[6]{}; // NB: only 3 are filled by sxxxxx, but 6 are compared!
        const int ipagM = ievt/neppM; // #eventpage in this iteration
        sxxxxx( hstMomenta.get(), smass, nhel, nsp, outwf3, ipagM, ipar );
        testwf6( outwf3, "sxxxxx", ievt, nsp, smass );
        sxxxxx( hstMomenta.get(), -smass, nhel, nsp, outwf3, ipagM, ipar );
        testwf6( outwf3, "sxxxxx", ievt, nsp, -smass );
      }
      // Test oxxxxx - NO ASSUMPTIONS
      {
        const fptype fmass = mass0[ievt];
        const int ipagM = ievt/neppM; // #eventpage in this iteration
        oxxxxx( hstMomenta.get(), fmass, nhel, nsp, outwfO, ipagM, ipar );
        testwf6( outwfO, "oxxxxx", ievt, nsp, fmass );
        oxxxxx( hstMomenta.get(), -fmass, nhel, nsp, outwfO, ipagM, ipar );
        testwf6( outwfO, "oxxxxx", ievt, nsp, -fmass );
      }
      // Test opzxxx - ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
      if ( mass0[ievt] == 0 && !isptgt0[ievt] && ispzgt0[ievt] )
      {
        const int ipagM = ievt/neppM; // #eventpage in this iteration
        opzxxx( hstMomenta.get(), nhel, nsp, outwf, ipagM, ipar );
        testwf6two( outwf, outwfO, "opzxxx", ievt );
        testwf6( outwf, "opzxxx", ievt, nsp, 0 );
      }
      // Test omzxxx - ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
      if ( mass0[ievt] == 0 && !isptgt0[ievt] && ispzlt0[ievt] )
      {
        const int ipagM = ievt/neppM; // #eventpage in this iteration
        omzxxx( hstMomenta.get(), nhel, nsp, outwf, ipagM, ipar );
        testwf6two( outwf, outwfO, "omzxxx", ievt );
        testwf6( outwf, "omzxxx", ievt, nsp, 0 );
      }
      // Test oxzxxx - ASSUMPTIONS: (FMASS == 0) and (PT > 0)
      if ( mass0[ievt] == 0 && isptgt0[ievt] )
      {
        const int ipagM = ievt/neppM; // #eventpage in this iteration
        oxzxxx( hstMomenta.get(), nhel, nsp, outwf, ipagM, ipar );
        testwf6two( outwf, outwfO, "oxzxxx", ievt );
        testwf6( outwf, "oxzxxx", ievt, nsp, 0 );
      }
    }
  }
  if ( dumpEvents )
  {
    dumpFile.close();
    std::cout << "INFO: New reference data dumped to file '" << dumpFileName << "'" << std::endl;
  }
}

//==========================================================================

// This is needed if and only if C++ LTO-like inlining optimizations are used in CPPProcess.cc (issue #229)
#ifdef MGONGPU_INLINE_HELAMPS
#include "../../src/HelAmps_sm.cc"
#endif

//==========================================================================

