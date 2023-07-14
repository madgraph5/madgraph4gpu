// Copyright (C) 2020-2023 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Apr 2021) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2021-2023) for the MG5aMC CUDACPP plugin.
//----------------------------------------------------------------------------
// Use ./runTest.exe --gtest_filter=*xxx to run only testxxx.cc tests
//----------------------------------------------------------------------------

#include "mgOnGpuConfig.h"

#include "CPPProcess.h"
#include "HelAmps_sm.h"
#include "MemoryAccessMomenta.h"
#include "MemoryAccessWavefunctions.h"
#include "MemoryBuffers.h"
#include "epoch_process_id.h"

#include <gtest/gtest.h>

#include <array>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#ifdef __CUDACC__
#define TESTID( s ) s##_GPU_XXX
#else
#define TESTID( s ) s##_CPU_XXX
#endif

#define XTESTID( s ) TESTID( s )

TEST( XTESTID( MG_EPOCH_PROCESS_ID ), testxxx )
{
  constexpr bool dumpEvents = false;       // dump the expected output of the test?
  constexpr bool testEvents = !dumpEvents; // run the test?
  constexpr fptype toleranceXXXs = std::is_same<fptype, double>::value ? 1.E-15 : 1.E-5;
  // Constant parameters
  constexpr int neppM = MemoryAccessMomenta::neppM; // AOSOA layout
  using mgOnGpu::neppV;
  constexpr int np4 = CPPProcess::np4;
  const int nevt = 16;         // 12 independent tests plus 4 duplicates (need a multiple of 8 for floats or for '512z')
  assert( nevt % neppM == 0 ); // nevt must be a multiple of neppM
  assert( nevt % neppV == 0 ); // nevt must be a multiple of neppV
  // Fill in the input momenta
#ifdef __CUDACC__
  mg5amcGpu::PinnedHostBufferMomenta hstMomenta( nevt ); // AOSOA[npagM][npar=4][np4=4][neppM]
#else
  mg5amcCpu::HostBufferMomenta hstMomenta( nevt ); // AOSOA[npagM][npar=4][np4=4][neppM]
#endif /* clang-format off */
  const fptype par0[np4 * nevt] = // AOS[nevt][np4]
    {
      500, 0, 0, 500,      // #0 (m=0 pT=0 E=pz>0)
      500, 0, 0, -500,     // #1 (m=0 pT=0 -E=pz<0)
      500, 300, 400, 0,    // #2 (m=0 pT>0 pz=0)
      500, 180, 240, 400,  // #3 (m=0 pT>0 pz>0)
      500, 180, 240, -400, // #4 (m=0 pT>0 pz<0)
      500, 0, 0, 0,        // #5 (m=50>0 pT=0 pz=0)
      500, 0, 0, 300,      // #6 (m=40>0 pT=0 pz>0)
      500, 0, 0, -300,     // #7 (m=40>0 pT=0 pz<0)
      500, 180, 240, 0,    // #8 (m=40>0 pT>0 pz=0)
      500, -240, -180, 0,  // #9 (m=40>0 pT>0 pz=0)
      500, 180, 192, 144,  // #10 (m=40>0 pT>0 pz>0)
      500, 180, 192, -144, // #11 (m=40>0 pT>0 pz<0)
      500, 0, 0, 500,      // DUPLICATE #12 == #0 (m=0 pT=0 E=pz>0)
      500, 0, 0, -500,     // DUPLICATE #13 == #1 (m=0 pT=0 -E=pz<0)
      500, 300, 400, 0,    // DUPLICATE #14 == #2 (m=0 pT>0 pz=0)
      500, 180, 240, 400   // DUPLICATE #15 == #3 (m=0 pT>0 pz>0)
    }; /* clang-format on */
  // Array initialization: zero-out as "{0}" (C and C++) or as "{}" (C++ only)
  // See https://en.cppreference.com/w/c/language/array_initialization#Notes
  fptype mass0[nevt] = {};
  bool ispzgt0[nevt] = {};
  bool ispzlt0[nevt] = {};
  bool isptgt0[nevt] = {};
  for( int ievt = 0; ievt < nevt; ievt++ )
  {
    const fptype p0 = par0[ievt * np4 + 0];
    const fptype p1 = par0[ievt * np4 + 1];
    const fptype p2 = par0[ievt * np4 + 2];
    const fptype p3 = par0[ievt * np4 + 3];
    mass0[ievt] = sqrt( p0 * p0 - p1 * p1 - p2 * p2 - p3 * p3 );
    ispzgt0[ievt] = ( p3 > 0 );
    ispzlt0[ievt] = ( p3 < 0 );
    isptgt0[ievt] = ( p1 != 0 ) || ( p2 != 0 );
  }
  const int ipar0 = 0; // use only particle0 for this test
  for( int ievt = 0; ievt < nevt; ievt++ )
  {
    for( int ip4 = 0; ip4 < np4; ip4++ )
    {
      MemoryAccessMomenta::ieventAccessIp4Ipar( hstMomenta.data(), ievt, ip4, ipar0 ) = par0[ievt * np4 + ip4]; // AOS to AOSOA
    }
  }
  // Expected output wavefunctions
  std::vector<std::array<fptype, 12>> expwfs;
#include "testxxx_cc_ref.txt" // expwfs.push_back( {...} );
  std::string dumpFileName = "testxxx_cc_ref.txt.new";
  // Compute the output wavefunctions
  // Dump new reference file if requested
  constexpr int nw6 = CPPProcess::nw6; // dimensions of each wavefunction (HELAS KEK 91-11): e.g. 6 for e+ e- -> mu+ mu- (fermions and vectors)
  int itest = 0;                       // index on the expected output vector
  std::ofstream dumpFile;
  if( dumpEvents ) dumpFile.open( dumpFileName, std::ios::trunc );
  auto dumpwf6 = [&]( std::ostream& out, const cxtype_sv wf[6], const char* xxx, int ievt, int nsp, fptype mass )
  {
    out << std::setprecision( 15 ) << std::scientific;
    out << "  expwfs.push_back( {";
    out << "                                   // ---------" << std::endl;
    for( int iw6 = 0; iw6 < nw6; iw6++ )
    {
#ifdef MGONGPU_CPPSIMD
      const int ieppV = ievt % neppV; // #event in the current event vector in this iteration
#ifdef MGONGPU_HAS_CPPCXTYPEV_BRK
      out << std::setw( 26 ) << cxreal( wf[iw6][ieppV] ) << ", ";
      out << std::setw( 22 ) << cximag( wf[iw6][ieppV] );
#else
      out << std::setw( 26 ) << wf[iw6].real()[ieppV] << ", ";
      out << std::setw( 22 ) << wf[iw6].imag()[ieppV];
#endif
#else
      out << std::setw( 26 ) << wf[iw6].real();
      out << ", " << std::setw( 22 ) << wf[iw6].imag();
#endif
      if( iw6 < nw6 - 1 )
        out << ",    ";
      else
        out << " } );";
      out << " // itest=" << itest << ": " << xxx << "#" << ievt;
      out << " nsp=" << nsp << " mass=" << (int)mass << std::endl;
    }
    out << std::defaultfloat;
  };
  auto testwf6 = [&]( const cxtype_sv wf[6], const char* xxx, int ievt, int nsp, fptype mass )
  {
    if( dumpEvents ) dumpwf6( dumpFile, wf, xxx, ievt, nsp, mass );
    if( testEvents )
    {
      std::array<fptype, 12>& expwf = expwfs[itest];
      //std::cout << "Testing " << std::setw(3) << itest << ": " << xxx << " #" << ievt << std::endl;
      ////for ( int iw6 = 0; iw6<nw6; iw6++ ) std::cout << wf[iw6] << std::endl;
      ////std::cout << "against" << std::endl;
      ////for ( int iw6 = 0; iw6<nw6; iw6++ )
      ////  std::cout << "[" << expwf[iw6*2] << "," << expwf[iw6*2+1] << "]" << std::endl; // NB: expwf[iw6*2], expwf[iw6*2+1] are fp
      for( int iw6 = 0; iw6 < nw6; iw6++ )
      {
        const fptype expReal = expwf[iw6 * 2];
        const fptype expImag = expwf[iw6 * 2 + 1];
        if( true )
        {
#ifdef MGONGPU_CPPSIMD
          const int ieppV = ievt % neppV; // #event in the current event vector in this iteration
#ifdef MGONGPU_HAS_CPPCXTYPEV_BRK
          EXPECT_NEAR( cxreal( wf[iw6][ieppV] ), expReal, std::abs( expReal * toleranceXXXs ) )
            << " itest=" << itest << ": " << xxx << "#" << ievt;
          EXPECT_NEAR( cximag( wf[iw6][ieppV] ), expImag, std::abs( expImag * toleranceXXXs ) )
            << " itest=" << itest << ": " << xxx << "#" << ievt;
#else
          EXPECT_NEAR( wf[iw6].real()[ieppV], expReal, std::abs( expReal * toleranceXXXs ) )
            << " itest=" << itest << ": " << xxx << "#" << ievt;
          EXPECT_NEAR( wf[iw6].imag()[ieppV], expImag, std::abs( expImag * toleranceXXXs ) )
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
  auto testwf6two = [&]( const cxtype_sv wf[6], const cxtype_sv expwf[6], const char* xxx, int ievt )
  {
    if( testEvents )
    {
      const std::string xxxFull( xxx[0] == 'i' ? "ixxxxx" : "oxxxxx" );
      //std::cout << "Testing " << std::setw(3) << itest << ": ";
      //std::cout << xxx << " #" << ievt << " against " << xxxFull << std::endl;
      ////for ( int iw6 = 0; iw6<nw6; iw6++ ) std::cout << wf[iw6] << std::endl;
      ////std::cout << "against" << std::endl;
      ////for ( int iw6 = 0; iw6<nw6; iw6++ ) std::cout << expwf[iw6] << std::endl; // NB: expwf[iw6] is cx
      for( int iw6 = 0; iw6 < nw6; iw6++ )
      {
        if( true )
        {
#ifdef MGONGPU_CPPSIMD
          const int ieppV = ievt % neppV; // #event in the current event vector in this iteration
#ifdef MGONGPU_HAS_CPPCXTYPEV_BRK
          const fptype expReal = cxreal( expwf[iw6][ieppV] );
          const fptype expImag = cximag( expwf[iw6][ieppV] );
          EXPECT_NEAR( cxreal( wf[iw6][ieppV] ), expReal, std::abs( expReal * toleranceXXXs ) )
            << " itest=" << itest << ": " << xxx << "#" << ievt << " against " << xxxFull;
          EXPECT_NEAR( cximag( wf[iw6][ieppV] ), expImag, std::abs( expImag * toleranceXXXs ) )
            << " itest=" << itest << ": " << xxx << "#" << ievt << " against " << xxxFull;
#else
          const fptype expReal = expwf[iw6].real()[ieppV];
          const fptype expImag = expwf[iw6].imag()[ieppV];
          EXPECT_NEAR( wf[iw6].real()[ieppV], expReal, std::abs( expReal * toleranceXXXs ) )
            << " itest=" << itest << ": " << xxx << "#" << ievt << " against " << xxxFull;
          EXPECT_NEAR( wf[iw6].imag()[ieppV], expImag, std::abs( expImag * toleranceXXXs ) )
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
  // Array initialization: zero-out as "{0}" (C and C++) or as "{}" (C++ only)
  // See https://en.cppreference.com/w/c/language/array_initialization#Notes
  cxtype_sv outwfI[6] = {}; // last result of ixxxxx (mass==0)
  cxtype_sv outwfO[6] = {}; // last result of oxxxxx (mass==0)
  cxtype_sv outwf[6] = {};
  cxtype_sv outwf3[6] = {};                                // NB: only 3 are filled by sxxxxx, but 6 are compared!
  fptype* fp_outwfI = reinterpret_cast<fptype*>( outwfI ); // proof of concept for using fptype* in the interface
  fptype* fp_outwfO = reinterpret_cast<fptype*>( outwfO ); // proof of concept for using fptype* in the interface
  fptype* fp_outwf = reinterpret_cast<fptype*>( outwf );   // proof of concept for using fptype* in the interface
  fptype* fp_outwf3 = reinterpret_cast<fptype*>( outwf3 ); // proof of concept for using fptype* in the interface
  const int nhel = 1;
  for( auto nsp: { -1, +1 } ) // antifermion/fermion (or initial/final for scalar and vector)
  {
    for( int ievt = 0; ievt < nevt; ievt++ )
    {
#ifdef __CUDACC__
      using namespace mg5amcGpu;
#else
      using namespace mg5amcCpu;
#endif
      if( false )
      {
        std::cout << std::endl;
        for( int ip4 = 0; ip4 < np4; ip4++ ) std::cout << par0[ievt * np4 + ip4] << ", ";
        std::cout << std::endl;
      }
      const int ipagV = ievt / neppV; // #event vector in this iteration
      const fptype* ievt0Momenta = MemoryAccessMomenta::ieventAccessRecordConst( hstMomenta.data(), ipagV * neppV );
      // Test ixxxxx - NO ASSUMPTIONS
      {
        const fptype fmass = mass0[ievt];
        ixxxxx<HostAccessMomenta, HostAccessWavefunctions>( ievt0Momenta, fmass, nhel, nsp, fp_outwfI, ipar0 );
        testwf6( outwfI, "ixxxxx", ievt, nsp, fmass );
        ixxxxx<HostAccessMomenta, HostAccessWavefunctions>( ievt0Momenta, -fmass, nhel, nsp, fp_outwfI, ipar0 );
        testwf6( outwfI, "ixxxxx", ievt, nsp, -fmass );
      }
      // Test ipzxxx - ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
      if( mass0[ievt] == 0 && !isptgt0[ievt] && ispzgt0[ievt] )
      {
        ipzxxx<HostAccessMomenta, HostAccessWavefunctions>( ievt0Momenta, nhel, nsp, fp_outwf, ipar0 );
        testwf6two( outwf, outwfI, "ipzxxx", ievt );
        testwf6( outwf, "ipzxxx", ievt, nsp, 0 );
      }
      // Test imzxxx - ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
      if( mass0[ievt] == 0 && !isptgt0[ievt] && ispzlt0[ievt] )
      {
        imzxxx<HostAccessMomenta, HostAccessWavefunctions>( ievt0Momenta, nhel, nsp, fp_outwf, ipar0 );
        testwf6two( outwf, outwfI, "imzxxx", ievt );
        testwf6( outwf, "imzxxx", ievt, nsp, 0 );
      }
      // Test ixzxxx - ASSUMPTIONS: (FMASS == 0) and (PT > 0)
      if( mass0[ievt] == 0 && isptgt0[ievt] )
      {
        ixzxxx<HostAccessMomenta, HostAccessWavefunctions>( ievt0Momenta, nhel, nsp, fp_outwf, ipar0 );
        testwf6two( outwf, outwfI, "ixzxxx", ievt );
        testwf6( outwf, "ixzxxx", ievt, nsp, 0 );
      }
      // Test vxxxxx - NO ASSUMPTIONS
      {
        const fptype vmass = mass0[ievt];
        vxxxxx<HostAccessMomenta, HostAccessWavefunctions>( ievt0Momenta, vmass, nhel, nsp, fp_outwf, ipar0 );
        testwf6( outwf, "vxxxxx", ievt, nsp, vmass );
        vxxxxx<HostAccessMomenta, HostAccessWavefunctions>( ievt0Momenta, -vmass, nhel, nsp, fp_outwf, ipar0 );
        testwf6( outwf, "vxxxxx", ievt, nsp, -vmass );
      }
      // Test sxxxxx - NO ASSUMPTIONS
      {
        const fptype smass = mass0[ievt];
        sxxxxx<HostAccessMomenta, HostAccessWavefunctions>( ievt0Momenta, nsp, fp_outwf3, ipar0 ); // no mass, no helicity (was "smass>0")
        testwf6( outwf3, "sxxxxx", ievt, nsp, smass );
        sxxxxx<HostAccessMomenta, HostAccessWavefunctions>( ievt0Momenta, nsp, fp_outwf3, ipar0 ); // no mass, no helicity (was "smass<0")
        testwf6( outwf3, "sxxxxx", ievt, nsp, -smass );
      }
      // Test oxxxxx - NO ASSUMPTIONS
      {
        const fptype fmass = mass0[ievt];
        oxxxxx<HostAccessMomenta, HostAccessWavefunctions>( ievt0Momenta, fmass, nhel, nsp, fp_outwfO, ipar0 );
        testwf6( outwfO, "oxxxxx", ievt, nsp, fmass );
        oxxxxx<HostAccessMomenta, HostAccessWavefunctions>( ievt0Momenta, -fmass, nhel, nsp, fp_outwfO, ipar0 );
        testwf6( outwfO, "oxxxxx", ievt, nsp, -fmass );
      }
      // Test opzxxx - ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
      if( mass0[ievt] == 0 && !isptgt0[ievt] && ispzgt0[ievt] )
      {
        opzxxx<HostAccessMomenta, HostAccessWavefunctions>( ievt0Momenta, nhel, nsp, fp_outwf, ipar0 );
        testwf6two( outwf, outwfO, "opzxxx", ievt );
        testwf6( outwf, "opzxxx", ievt, nsp, 0 );
      }
      // Test omzxxx - ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
      if( mass0[ievt] == 0 && !isptgt0[ievt] && ispzlt0[ievt] )
      {
        omzxxx<HostAccessMomenta, HostAccessWavefunctions>( ievt0Momenta, nhel, nsp, fp_outwf, ipar0 );
        testwf6two( outwf, outwfO, "omzxxx", ievt );
        testwf6( outwf, "omzxxx", ievt, nsp, 0 );
      }
      // Test oxzxxx - ASSUMPTIONS: (FMASS == 0) and (PT > 0)
      if( mass0[ievt] == 0 && isptgt0[ievt] )
      {
        oxzxxx<HostAccessMomenta, HostAccessWavefunctions>( ievt0Momenta, nhel, nsp, reinterpret_cast<fptype*>( outwf ), ipar0 );
        testwf6two( outwf, outwfO, "oxzxxx", ievt );
        testwf6( outwf, "oxzxxx", ievt, nsp, 0 );
      }
    }
  }
  if( dumpEvents )
  {
    dumpFile.close();
    std::cout << "INFO: New reference data dumped to file '" << dumpFileName << "'" << std::endl;
  }
}

//==========================================================================
