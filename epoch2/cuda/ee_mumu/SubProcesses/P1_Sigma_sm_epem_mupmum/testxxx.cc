#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"
//#include "mgOnGpuVectors.h"

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
  constexpr bool testEvents = true; // run the test?
  constexpr fptype toleranceXXXs = std::is_same<fptype, double>::value ? 1.E-15 : 1.E-5;
  // Constant parameters
  using mgOnGpu::neppM;
  using mgOnGpu::np4;
  using mgOnGpu::npar;
  const int nevt = 12;
  assert( nevt % neppM == 0 ); // nevt must be a multiple of neppM
  // Fill in the input momenta
  const int nMomenta = np4 * npar * nevt;
  //auto hstMomenta = hstMakeUnique<fptype_sv>( nMomenta ); // AOSOA[npagM][npar=4][np4=4][neppM]
  auto hstMomenta = hstMakeUnique<fptype>( nMomenta ); // AOSOA[npagM][npar=4][np4=4][neppM]
  const fptype par0[np4 * nevt]{
    // AOS[nevt][np4]
    500, 0,    0,    500,  // #0 (m=0 pT=0 E=pz>0)
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
  };
  fptype mass0[nevt]{};
  bool ispzgt0[nevt]{};
  bool ispzlt0[nevt]{};
  bool isptgt0[nevt]{};
  for ( int ievt = 0; ievt < nevt; ievt++ )
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
  const int ipar = 0; // use only particle0 for this test
  for ( int ievt = 0; ievt < nevt; ievt++ )
  {
    for ( int ip4 = 0; ip4 < np4; ip4++ )
    {
      const int ipagM = ievt / neppM; // #eventpage in this iteration
      const int ieppM = ievt % neppM; // #event in the current eventpage in this iteration
      //#ifdef MGONGPU_CPPSIMD
      //hstMomenta[ipagM*npar*np4 + ipar*np4 + ip4][ieppM] = par0[ievt*np4 + ip4]; // AOS to AOSOA
      //#else
      hstMomenta[ipagM * npar * np4 * neppM + ipar * np4 * neppM + ip4 * neppM + ieppM] = par0[ievt * np4 + ip4]; // AOS to AOSOA
      //#endif
    }
  }
  // Expected output wavefunctions
  std::vector<std::array<fptype, 12>> expwfs;
#include "testxxx_cc_ref.txt" // expwfs.push_back( {...} );
  std::string dumpFileName = "testxxx_cc_ref.txt";
  // Compute the output wavefunctions
  // Dump new reference file if requested
  using namespace MG5_sm;
  const int nwf6 = 6;
  int itest = 0; // index on the expected output vector
  std::ofstream dumpFile;
  if ( dumpEvents ) dumpFile.open( dumpFileName, std::ios::trunc );
  auto dumpwf6 = [&]( const cxtype wf[6], const char* xxx, int ievt, int nsp, fptype mass ) {
    dumpFile << std::setprecision( 15 ) << std::scientific;
    dumpFile << "  expwfs.push_back( {";
    dumpFile << "                                   // ---------" << std::endl;
    for ( int iwf6 = 0; iwf6 < nwf6; iwf6++ )
    {
      dumpFile << std::setw( 26 ) << wf[iwf6].real() << ", " << std::setw( 22 ) << wf[iwf6].imag();
      if ( iwf6 < nwf6 - 1 )
        dumpFile << ",    ";
      else
        dumpFile << " } );";
      dumpFile << " // itest=" << itest << ": " << xxx << "#" << ievt;
      dumpFile << " nsp=" << nsp << " mass=" << (int)mass << std::endl;
    }
    dumpFile << std::defaultfloat;
  };
  auto testwf6 = [&]( const cxtype wf[6], const char* xxx, int ievt, int nsp, fptype mass ) {
    if ( dumpEvents ) dumpwf6( wf, xxx, ievt, nsp, mass );
    if ( testEvents )
    {
      //std::cout << "Testing " << itest << ": " << xxx << " #" << ievt << std::endl;
      std::array<fptype, 12>& expwf = expwfs[itest];
      for ( int iwf6 = 0; iwf6 < nwf6; iwf6++ )
      {
        const fptype expReal = expwf[iwf6*2];
        const fptype expImag = expwf[iwf6*2+1];
        EXPECT_NEAR( cxreal( wf[iwf6] ), expReal, std::abs( expReal * toleranceXXXs ) )
          << " itest=" << itest << ": " << xxx << "#" << ievt;
        EXPECT_NEAR( cximag( wf[iwf6] ), expImag, std::abs( expImag * toleranceXXXs ) )
          << " itest=" << itest << ": " << xxx << "#" << ievt;
      }
    }
    itest++;
  };
  const int nhel = 1;
  cxtype outwf[6];
  for ( auto nsp : { -1, +1 } ) // antifermion/fermion (or initial/final for scalar and vector)
  {
    for ( int ievt = 0; ievt < nevt; ievt++ )
    {
      if ( false )
      {
        std::cout << std::endl;
        for ( int ip4 = 0; ip4 < np4; ip4++ ) std::cout << par0[ievt * np4 + ip4] << ", ";
        std::cout << std::endl;
      }
      // Test ixxxxx - NO ASSUMPTIONS
      {
        const fptype fmass = mass0[ievt];
        ixxxxx( hstMomenta.get(), fmass, nhel, nsp, outwf, ievt, ipar );
        testwf6( outwf, "ixxxxx", ievt, nsp, fmass );
        ixxxxx( hstMomenta.get(), -fmass, nhel, nsp, outwf, ievt, ipar );
        testwf6( outwf, "ixxxxx", ievt, nsp, -fmass );
      }
      // Test ipzxxx - ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
      if ( mass0[ievt] == 0 && !isptgt0[ievt] && ispzgt0[ievt] )
      {
        ipzxxx( hstMomenta.get(), nhel, nsp, outwf, ievt, ipar );
        testwf6( outwf, "ipzxxx", ievt, nsp, 0 );
      }
      // Test imzxxx - ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
      if ( mass0[ievt] == 0 && !isptgt0[ievt] && ispzlt0[ievt] )
      {
        imzxxx( hstMomenta.get(), nhel, nsp, outwf, ievt, ipar );
        testwf6( outwf, "imzxxx", ievt, nsp, 0 );
      }
      // Test ixzxxx - ASSUMPTIONS: (FMASS == 0) and (PT > 0)
      if ( mass0[ievt] == 0 && isptgt0[ievt] )
      {
        ixzxxx( hstMomenta.get(), nhel, nsp, outwf, ievt, ipar );
        testwf6( outwf, "ixzxxx", ievt, nsp, 0 );
      }
      // Test vxxxxx - NO ASSUMPTIONS
      {
        const fptype vmass = mass0[ievt];
        vxxxxx( hstMomenta.get(), vmass, nhel, nsp, outwf, ievt, ipar );
        testwf6( outwf, "vxxxxx", ievt, nsp, vmass );
        vxxxxx( hstMomenta.get(), -vmass, nhel, nsp, outwf, ievt, ipar );
        testwf6( outwf, "vxxxxx", ievt, nsp, -vmass );
      }
      // Test sxxxxx - NO ASSUMPTIONS
      {
        const fptype smass = mass0[ievt];
        sxxxxx( hstMomenta.get(), smass, nhel, nsp, outwf, ievt, ipar );
        testwf6( outwf, "sxxxxx", ievt, nsp, smass );
        sxxxxx( hstMomenta.get(), -smass, nhel, nsp, outwf, ievt, ipar );
        testwf6( outwf, "sxxxxx", ievt, nsp, -smass );
      }
      // Test oxxxxx - NO ASSUMPTIONS
      {
        const fptype fmass = mass0[ievt];
        oxxxxx( hstMomenta.get(), fmass, nhel, nsp, outwf, ievt, ipar );
        testwf6( outwf, "oxxxxx", ievt, nsp, fmass );
        oxxxxx( hstMomenta.get(), -fmass, nhel, nsp, outwf, ievt, ipar );
        testwf6( outwf, "oxxxxx", ievt, nsp, -fmass );
      }
      // Test opzxxx - ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
      if ( mass0[ievt] == 0 && !isptgt0[ievt] && ispzgt0[ievt] )
      {
        opzxxx( hstMomenta.get(), nhel, nsp, outwf, ievt, ipar );
        testwf6( outwf, "opzxxx", ievt, nsp, 0 );
      }
      // Test omzxxx - ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
      if ( mass0[ievt] == 0 && !isptgt0[ievt] && ispzlt0[ievt] )
      {
        omzxxx( hstMomenta.get(), nhel, nsp, outwf, ievt, ipar );
        testwf6( outwf, "omzxxx", ievt, nsp, 0 );
      }
      // Test oxzxxx - ASSUMPTIONS: (FMASS == 0) and (PT > 0)
      if ( mass0[ievt] == 0 && isptgt0[ievt] )
      {
        oxzxxx( hstMomenta.get(), nhel, nsp, outwf, ievt, ipar );
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
