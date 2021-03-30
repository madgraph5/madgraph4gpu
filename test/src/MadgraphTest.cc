/*
 * MadgraphTest.cc
 *
 *  Created on: 11.12.2020
 *      Author: shageboeck
 */

#include "MadgraphTest.h"

#include <cmath>
#include <fstream>
#include <sstream>

std::map<unsigned int, ReferenceData> readReferenceData(const std::string& refFileName)
{
  std::ifstream referenceFile(refFileName.c_str());
  EXPECT_TRUE(referenceFile.is_open()) << refFileName;
  std::map<unsigned int, ReferenceData> referenceData;
  unsigned int evtNo;
  unsigned int batchNo;

  for (std::string line; std::getline(referenceFile, line); )
  {
    std::stringstream lineStr(line);
    if (line.empty() || line[0] == '#')
    {
      continue;
    }
    else if (line.find("Event") != std::string::npos)
    {
      std::string dummy;
      lineStr >> dummy >> evtNo >> dummy >> batchNo;
    }
    else if (line.find("ME") != std::string::npos)
    {
      if (evtNo <= referenceData[batchNo].MEs.size())
        referenceData[batchNo].MEs.resize(evtNo + 1);

      std::string dummy;
      lineStr >> dummy >> referenceData[batchNo].MEs[evtNo];
    }
    else
    {
      unsigned int particleIndex;
      lineStr >> particleIndex;

      if (evtNo <= referenceData[batchNo].momenta.size())
        referenceData[batchNo].momenta.resize(evtNo + 1);
      if (particleIndex <= referenceData[batchNo].momenta[evtNo].size())
        referenceData[batchNo].momenta[evtNo].resize(particleIndex + 1);

      auto& fourVec = referenceData[batchNo].momenta[evtNo][particleIndex];
      for (unsigned int i=0; i < fourVec.size(); ++i) {
        EXPECT_TRUE(lineStr.good());
        lineStr >> fourVec[i];
      }
      EXPECT_TRUE(lineStr.eof());
    }
  }
  return referenceData;
}


TEST_P(MadgraphTestDouble, eemumu)
{
  // Set to dump events:
  constexpr bool dumpEvents = false;
  constexpr fptype toleranceMomenta = std::is_same<fptype, double>::value ? 5.E-12 : 1.E-5;
  constexpr fptype toleranceMEs     = std::is_same<fptype, double>::value ? 1.E-6  : 1.E-5;
  constexpr fptype energy = 1500; // historical default, Ecms = 1500 GeV = 1.5 TeV (above the Z peak)

  std::string dumpFileName = std::string("dump_")
      + testing::UnitTest::GetInstance()->current_test_info()->name()
      + ".txt";
  while (dumpFileName.find('/') != std::string::npos) {
    dumpFileName.replace(dumpFileName.find('/'), 1, "_");
  }
  const std::string refFileName = "../../../../../test/eemumu/dump_CPUTest.eemumu.txt";

  std::ofstream dumpFile;
  if ( dumpEvents )
  {
    dumpFile.open(dumpFileName, std::ios::trunc);
  }

  // Read reference data
  std::map<unsigned int, ReferenceData> referenceData = readReferenceData(refFileName);
  ASSERT_FALSE(HasFailure()); // It doesn't make any sense to continue if we couldn't read the reference file.


  // **************************************
  // *** START MAIN LOOP ON #ITERATIONS ***
  // **************************************
  for (unsigned int iiter = 0; iiter < testDriver->niter; ++iiter)
  {
    testDriver->prepareRandomNumbers(iiter);

    testDriver->prepareMomenta(energy);

    testDriver->runSigmaKin(iiter);

    // --- Run checks on all events produced in this iteration
    for (std::size_t ievt = 0; ievt < testDriver->nevt && !HasFailure(); ++ievt)
    {
      if (dumpEvents) {
        ASSERT_TRUE(dumpFile.is_open()) << dumpFileName;
        dumpFile << "Event " << std::setw(8) << ievt << "  "
                 << "Batch " << std::setw(4) << iiter << "\n";
        testDriver->dumpParticles(dumpFile, ievt, testDriver->nparticle, 15, ReferenceData());
        // Dump matrix element
        dumpFile << std::setw(4) << "ME" << std::scientific << std::setw(15+8)
            << testDriver->getMatrixElement(ievt) << "\n" << std::endl << std::defaultfloat;
        continue;
      }


      // Check that we have the required reference data
      ASSERT_GT(referenceData.size(), iiter) << "Don't have enough reference data for iteration " << iiter << ". Ref file:" << refFileName;
      ASSERT_GT(referenceData[iiter].MEs.size(), ievt)     << "Don't have enough reference MEs for iteration " << iiter << " event " << ievt << ".\nRef file: " << refFileName;
      ASSERT_GT(referenceData[iiter].momenta.size(), ievt) << "Don't have enough reference momenta for iteration " << iiter << " event " << ievt << ".\nRef file: " << refFileName;
      ASSERT_GE(referenceData[iiter].momenta[ievt].size(), testDriver->nparticle) << "Don't have enough reference particles for iteration " << iiter << " event " << ievt << ".\nRef file: " << refFileName;


      // This trace will help to understand the event that is being checked.
      // It will only be printed in case of failures:
      std::stringstream eventTrace;
      eventTrace << "In comparing event " << ievt << " from iteration " << iiter << "\n";
      testDriver->dumpParticles(eventTrace, ievt, testDriver->nparticle, 15, referenceData[iiter]);
      eventTrace << std::setw(4) << "ME"   << std::scientific << std::setw(15+8) << testDriver->getMatrixElement(ievt) << "\n"
                 << std::setw(4) << "r.ME" << std::scientific << std::setw(15+8) << referenceData[iiter].MEs[ievt] << std::endl << std::defaultfloat;
      SCOPED_TRACE(eventTrace.str());


      // Compare Momenta
      for (unsigned int ipar = 0; ipar < testDriver->nparticle; ++ipar) {
        std::stringstream momentumErrors;
        for (unsigned int icomp = 0; icomp < testDriver->np4; ++icomp) {
          const double pMadg = testDriver->getMomentum(ievt, ipar, icomp);
          const double pOrig = referenceData[iiter].momenta[ievt][ipar][icomp];
          const double relDelta = fabs( (pMadg - pOrig)/pOrig );
          if (relDelta > toleranceMomenta) {
            momentumErrors << std::setprecision(15) << std::scientific << "\nparticle " << ipar << "\tcomponent " << icomp
                << "\n\t madGraph:  " << std::setw(22) << pMadg
                << "\n\t reference: " << std::setw(22) << pOrig
                << "\n\t rel delta: " << std::setw(22) << relDelta << " exceeds tolerance of " << toleranceMomenta;
          }
        }
        ASSERT_TRUE(momentumErrors.str().empty()) << momentumErrors.str();
      }


      // Compare ME:
      EXPECT_NEAR(testDriver->getMatrixElement(ievt),
          referenceData[iiter].MEs[ievt],
          toleranceMEs * referenceData[iiter].MEs[ievt]);
    }
  }
}

