#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"
#include "mgOnGpuVectors.h"

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
#define TESTID_CPU(s) s##_CPU
#define XTESTID_CPU(s) TESTID_CPU(s)

TEST( XTESTID_CPU(MG_EPOCH_PROCESS_ID), testxxx )
{
  constexpr bool dumpEvents = false; // dump the expected output of the test?
  constexpr bool testEvents = true; // run the test?
  constexpr fptype toleranceXXXs = std::is_same<fptype, double>::value ? 1.E-15 : 1.E-5;
  // Constant parameters
  using mgOnGpu::np4;
  using mgOnGpu::npar;
  using mgOnGpu::neppM;
  const int nevt = 8;
  assert( nevt%neppM == 0 ); // nevt must be a multiple of neppM
  // Fill in the input momenta
  const int nMomenta = np4*npar*nevt;
  auto hstMomenta = hstMakeUnique<fptype_sv>( nMomenta ); // AOSOA[npagM][npar=4][np4=4][neppM]
  const fptype par0[np4 * nevt]                           // AOS[nevt][np4]
    {
      500,  0,    0,    500,   // #0 (m=0 pT=0 E=pz>0)
      500,  0,    0,    -500,  // #1 (m=0 pT=0 E=-pz>0)
      500,  300,  400,  0,     // #2 (m=0 pT>0 pz=0)
      500,  180,  240,  400,   // #3 (m=0 pT>0 pz>0)
      500,  180,  240,  -400,  // #4 (m=0 pT>0 pz<0)
      500., 0.,   0.,   0.,    // #5 (m=50>0 pT=0 pz=0)
      500,  0,    0,    -300,  // #6 (m=40>0 pT=0 pz<0)
      500., 180., 192., -144., // #7 (m=40>0 pT>0 pz<0)
    };
  fptype mass0[nevt]{};
  bool ispzgt0[nevt]{};
  bool ispzlt0[nevt]{};
  bool isptgt0[nevt]{};
  for ( int ievt=0; ievt<nevt; ievt++ )
  {
    const fptype p0 = par0[ievt*np4 + 0];
    const fptype p1 = par0[ievt*np4 + 1];
    const fptype p2 = par0[ievt*np4 + 2];
    const fptype p3 = par0[ievt*np4 + 3];
    mass0[ievt] = sqrt( p0*p0 - p1*p1 - p2*p2 - p3*p3 );
    ispzgt0[ievt] = ( p3 > 0 );
    ispzlt0[ievt] = ( p3 < 0 );
    isptgt0[ievt] = ( p1 != 0 ) || ( p2 != 0 );
  }
  const int ipar=0; // use only particle0 for this test
  for ( int ievt=0; ievt<nevt; ievt++ )
  {
    for ( int ip4=0; ip4<np4; ip4++ )
    {
      const int ipagM = ievt/neppM; // #eventpage in this iteration
      const int ieppM = ievt%neppM; // #event in the current eventpage in this iteration
#ifdef MGONGPU_CPPSIMD
      hstMomenta[ipagM*npar*np4 + ipar*np4 + ip4][ieppM] = par0[ievt*np4 + ip4]; // AOS to AOSOA
#else
      hstMomenta[ipagM*npar*np4*neppM + ipar*np4*neppM + ip4*neppM + ieppM] = par0[ievt*np4 + ip4]; // AOS to AOSOA
#endif
    }
  }
  // Expected output wavefunctions
  std::vector< std::array<fptype,12> > expwfs;
#include "testxxx_cc_ref.txt" // expwfs.push_back( {...} );
  std::string dumpFileName = "testxxx_cc_ref.txt.new";
  // Compute the output wavefunctions
  // Dump new reference file if requested
  using namespace MG5_sm;
  const int nwf6 = 6;
  std::ofstream dumpFile;
  if ( dumpEvents ) dumpFile.open( dumpFileName, std::ios::trunc );
  auto dumpwf6 = [&]( const cxtype_sv wf[6], const char* xxx, int ievt ) {
    dumpFile << std::setprecision(15) << std::scientific;
    dumpFile << "  expwfs.push_back( {";
    dumpFile << "                                   // ---------" << std::endl;
    for ( int iwf6 = 0; iwf6<nwf6; iwf6++ )
    {
#ifdef MGONGPU_CPPSIMD
      const int ieppM = ievt%neppM; // #event in the current eventpage in this iteration
      dumpFile << std::setw(26) << cxreal( wf[iwf6][ieppM] ) << ", ";
      dumpFile << std::setw(22) << cximag( wf[iwf6][ieppM] );
#else
      dumpFile << std::setw(26) << wf[iwf6].real();
      dumpFile << ", " << std::setw(22) << wf[iwf6].imag();
#endif
      if ( iwf6 < nwf6-1 ) dumpFile << ",    ";
      else dumpFile << " } );";
      dumpFile << " // " << xxx << " #" << ievt << std::endl;
    }
    dumpFile << std::defaultfloat;
  };
  int itest = 0; // index on the expected output vector
  auto testwf6 = [&]( const cxtype_sv wf[6], const char* xxx, int ievt ) {
    if ( dumpEvents ) dumpwf6( wf, xxx, ievt );
    if ( testEvents )
    {
      std::cout << "Testing " << std::setw(3) << itest << ": ";
      std::cout << xxx << " #" << ievt << std::endl;
      std::array<fptype,12>& expwf = expwfs[itest];
#ifdef MGONGPU_CPPSIMD
      const int ieppM = ievt%neppM; // #event in the current eventpage in this iteration
      for ( int iwf6 = 0; iwf6<nwf6; iwf6++ )
      {
        EXPECT_NEAR( cxreal( wf[iwf6][ieppM] ), expwf[iwf6*2], std::abs( expwf[iwf6*2] * toleranceXXXs ) );
        EXPECT_NEAR( cximag( wf[iwf6][ieppM] ), expwf[iwf6*2+1], std::abs( expwf[iwf6*2+1] * toleranceXXXs ) );
      }
#else
      for ( int iwf6 = 0; iwf6<nwf6; iwf6++ )
      {
        EXPECT_NEAR( wf[iwf6].real(), expwf[iwf6*2], std::abs( expwf[iwf6*2] * toleranceXXXs ) );
        EXPECT_NEAR( wf[iwf6].imag(), expwf[iwf6*2+1], std::abs( expwf[iwf6*2+1] * toleranceXXXs ) );
      }
#endif
    }
    itest++;
  };
  const int ihel = +1;
  const int nsf = -1;
  //cxtype outwf[6];
  cxtype_sv outwf[6];
  for ( int ievt=0; ievt<nevt; ievt++ )
  {
    if ( false )
    {
      std::cout << std::endl;
      for ( int ip4=0; ip4<np4; ip4++ ) std::cout << par0[ievt*np4 + ip4] << ", ";
      std::cout << std::endl;
    }
    // Test ixxxxx - NO ASSUMPTIONS
    {
      //const fptype fmass = mass0[ievt];
      //ixxxxx( hstMomenta.get(), fmass, ihel, nsf, outwf, ievt, ipar );
      //testwf6( outwf, "ixxxxx", ievt );
      itest++; // SKIP
    }
    // Test ipzxxx - ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
    if ( mass0[ievt] == 0 && ispzgt0[ievt] )
    {
      //ipzxxx( hstMomenta.get(), ihel, nsf, outwf, ievt, ipar );
      //testwf6( outwf, "ipzxxx", ievt );
      itest++; // SKIP
    }
    // Test imzxxx - ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
    if ( mass0[ievt] == 0 && ispzlt0[ievt] )
    {
      const int ipagM = ievt/neppM; // #eventpage in this iteration
      imzxxx( hstMomenta.get(), ihel, nsf, outwf, ipagM, ipar );
      testwf6( outwf, "imzxxx", ievt );
    }
    // Test ixzxxx - ASSUMPTIONS: (FMASS == 0) and (PT > 0)
    if ( mass0[ievt] == 0 && isptgt0[ievt] )
    {
      const int ipagM = ievt/neppM; // #eventpage in this iteration
      ixzxxx( hstMomenta.get(), ihel, nsf, outwf, ipagM, ipar );
      testwf6( outwf, "ixzxxx", ievt );
    }
    // Test vxxxxx - NO ASSUMPTIONS
    {
      //const fptype vmass = mass0[ievt];
      //vxxxxx( hstMomenta.get(), vmass, ihel, nsf, outwf, ievt, ipar );
      //testwf6( outwf, "vxxxxx", ievt );
      itest++; // SKIP
    }
    // Test sxxxxx - NO ASSUMPTIONS
    {
      //const fptype smass = mass0[ievt];
      //sxxxxx( hstMomenta.get(), smass, ihel, nsf, outwf, ievt, ipar );
      //testwf6( outwf, "sxxxxx", ievt );
      itest++; // SKIP
    }
    // Test oxxxxx - NO ASSUMPTIONS
    {
      //const fptype fmass = mass0[ievt];
      //oxxxxx( hstMomenta.get(), fmass, ihel, nsf, outwf, ievt, ipar );
      //testwf6( outwf, "oxxxxx", ievt );
      itest++; // SKIP
    }
    // Test opzxxx - ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
    if ( mass0[ievt] == 0 && ispzgt0[ievt] )
    {
      const int ipagM = ievt/neppM; // #eventpage in this iteration
      opzxxx( hstMomenta.get(), ihel, nsf, outwf, ipagM, ipar );
      testwf6( outwf, "opzxxx", ievt );
    }
    // Test omzxxx - ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
    if ( mass0[ievt] == 0 && ispzlt0[ievt] )
    {
      //omzxxx( hstMomenta.get(), ihel, nsf, outwf, ievt, ipar );
      //testwf6( outwf, "omzxxx", ievt );
      itest++; // SKIP
    }
    // Test oxzxxx - ASSUMPTIONS: (FMASS == 0) and (PT > 0)
    if ( mass0[ievt] == 0 && isptgt0[ievt] )
    {
      const int ipagM = ievt/neppM; // #eventpage in this iteration
      oxzxxx( hstMomenta.get(), ihel, nsf, outwf, ipagM, ipar );
      testwf6( outwf, "oxzxxx", ievt );
    }
  }
  if ( dumpEvents )
  {
    dumpFile.close();
    std::cout << "INFO: New reference data dumped to file '" << dumpFileName << "'" << std::endl;
  }
}
