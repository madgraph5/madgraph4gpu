// Stephan Hageboeck, CERN, 12/2020
#ifndef MADGRAPHTEST_H_
#define MADGRAPHTEST_H_

#include <array>
#include <iomanip>
#include <string>
#include <vector>

#include <gtest/gtest.h>

struct ReferenceData {
  std::vector< std::vector<std::array<double, 4>> > momenta;
  std::vector<double> MEs;
};

std::map<unsigned int, ReferenceData> readReferenceData(const std::string& refFileName);

/**
 * Test driver providing a common interface for testing different implementations.
 * Users need to implement:
 * - Functions to retrieve matrix element and 4-momenta. These are used in the tests.
 * - Driver functions that run the madgraph workflow.
 *
 * Usage:
 * ```
 * class TestImplementation : public BaseTest {
 *   <implement functions>
 * }
 *
 * TEST_F(TestImplementation, <testName>) {
 *   <test code>
 * }
 */
template<typename FpType>
class TestDriverBase {
 public:
  using fptype = FpType;
  unsigned int nparticle = 4;
  static constexpr unsigned int np4 = 4;
  static constexpr unsigned int niter = 2;
  static constexpr unsigned int gpublocks = 2;
  static constexpr unsigned int gputhreads = 128;
  static constexpr unsigned int nevt = gpublocks * gputhreads;

  TestDriverBase() { }
  virtual ~TestDriverBase() { }


  // ------------------------------------------------
  // Interface for retrieving info from madgraph
  // ------------------------------------------------
  virtual double getMomentum(std::size_t evtNo, unsigned int particleNo, unsigned int component) const = 0;
  virtual double getMatrixElement(std::size_t evtNo) const = 0;


  // ------------------------------------------------
  // Interface for steering madgraph run
  // ------------------------------------------------
  virtual void prepareRandomNumbers(unsigned int iiter) = 0;
  virtual void prepareMomenta(fptype energy) = 0;
  virtual void runSigmaKin(std::size_t iiter) = 0;

  // Print the requested event into the stream. If the reference data has enough events, it will be printed as well.
  void dumpParticles(std::ostream& stream, std::size_t ievt, unsigned int numParticles, unsigned int nDigit, const ReferenceData& referenceData)
  {
    const auto width = nDigit + 8;
    for (unsigned int ipar = 0; ipar < numParticles; ipar++)
    {
      // NB: 'setw' affects only the next field (of any type)
      stream << std::scientific // fixed format: affects all floats (default nDigit: 6)
             << std::setprecision(nDigit)
             << std::setw(4) << ipar
             << std::setw(width) << getMomentum(ievt, ipar, 0)
             << std::setw(width) << getMomentum(ievt, ipar, 1)
             << std::setw(width) << getMomentum(ievt, ipar, 2)
             << std::setw(width) << getMomentum(ievt, ipar, 3)
             << "\n";
      if (ievt < referenceData.momenta.size()) {
        stream << "ref" << ipar;
        stream << std::setw(width) << referenceData.momenta[ievt][ipar][0]
            << std::setw(width) << referenceData.momenta[ievt][ipar][1]
            << std::setw(width) << referenceData.momenta[ievt][ipar][2]
            << std::setw(width) << referenceData.momenta[ievt][ipar][3]
            << "\n\n";
      }
      stream << std::flush << std::defaultfloat; // default format: affects all floats
    }
  };
};


// Test class that's using the driver to run the test(s) below.
class MadgraphTestDouble : public testing::TestWithParam<std::function<TestDriverBase<double>*()>> {
protected:
  using fptype = double;
  using TestDriver_t = TestDriverBase<fptype>;
  std::unique_ptr<TestDriverBase<fptype>> testDriver;


public:
  MadgraphTestDouble() :
    TestWithParam(),
    testDriver{ GetParam()() }
  { }
};



#endif /* MADGRAPHTEST_H_ */
