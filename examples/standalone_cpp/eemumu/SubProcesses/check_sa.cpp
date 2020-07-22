#include <algorithm> // perf stats
#include <cstring>
#include <iomanip>
#include <iostream>
#include <numeric> // perf stats

#include "CPPProcess.h"
#include "rambo.h"
#include <CL/sycl.hpp>

bool is_number(const char *s) {
  const char *t = s;
  while (*t != '\0' && isdigit(*t))
    ++t;
  return strlen(s) == t - s;
}

int usage(char* argv0, int ret = 1) {
  std::cout << "Usage: " << argv0 << " [--verbose] numevts" << std::endl;
  return ret;
}

int main(int argc, char **argv) {
  bool verbose = false;
  int numevts = 0;
  for (int argn = 1; argn < argc; ++argn) {
    if (strcmp(argv[argn], "--verbose") == 0)
      verbose = true;
    else if (is_number(argv[argn]))
      numevts = atoi(argv[argn]);
    else
      return usage(argv[0]);
  }
  if (numevts == 0)
    return usage(argv[0]);

  if (verbose)
    std::cout << "num evts: " << numevts << std::endl;

  // Create a process object
  CPPProcess process;
  
  // Used to store timings
  vector<float> t(numevts);

  cl::sycl::range<1> work_items{(unsigned long)numevts};
  cl::sycl::buffer<float> buff_t(t.data(), t.size());
  cl::sycl::queue q; 
 
  q.submit([&](cl::sycl::handler& cgh){

    auto access_t = buff_t.get_access<cl::sycl::access::mode::write>(cgh);
 
    // Put here so that it captures the instance     
    // Read param_card and set parameters
    process.initProc("../../Cards/param_card.dat", verbose);

    cgh.parallel_for<class my_kernel>(work_items,
                                         [=] (cl::sycl::id<1> idx) {
      // Create local (tread) copy (_lc) of the instance
      CPPProcess process_lc = process; 
      double energy = 1500;
      double weight;

      // Get phase space point
      vector<double *> p = 
      get_momenta(process_lc.ninitial, energy, process_lc.getMasses(), weight);

      // Set momenta for this event
      process_lc.setMomenta(p);

      // Evaluate matrix element
      process_lc.sigmaKin(false);

      const double *matrix_elements = process_lc.getMatrixElements();

      // Store the timings back in the orginal instance
      access_t[idx] = process_lc.m_wavetimes[0];

    });
  });
  q.wait();
  process.m_wavetimes = t;
  process.printPerformanceStats();

/* 
   if (verbose) {
      cout << "Momenta:" << endl;
      for (int i = 0; i < process.nexternal; i++)
        cout << setw(4) << i + 1 << setiosflags(ios::scientific) << setw(14)
             << p[i][0] << setiosflags(ios::scientific) << setw(14) << p[i][1]
             << setiosflags(ios::scientific) << setw(14) << p[i][2]
             << setiosflags(ios::scientific) << setw(14) << p[i][3] << endl;
      cout << " ---------------------------------------------------------------"
              "--------------"
           << endl;
    }

    if (verbose) {
      // Display matrix elements
      for (int i = 0; i < process.nprocesses; i++)
        cout << " Matrix element = "
             //	 << setiosflags(ios::fixed) << setprecision(17)
             << matrix_elements[i] << " GeV^" << meGeVexponent << std::endl;

      cout << " ---------------------------------------------------------------"
              "--------------"
           << endl;
    }
   
    for (int i = 0; i < process.nprocesses; i++)
      matrixelementvector.push_back(matrix_elements[i]);

  }
*/

  int num_mes = matrixelementvector.size();
  float sumelem = std::accumulate(matrixelementvector.begin(), matrixelementvector.end(), 0.0);
  float meanelem = sumelem / num_mes;
  float sqselem = std::inner_product(matrixelementvector.begin(), matrixelementvector.end(), 
                                     matrixelementvector.begin(), 0.0);
  float stdelem = std::sqrt(sqselem / num_mes - meanelem * meanelem);
  std::vector<double>::iterator maxelem =
    std::max_element(matrixelementvector.begin(), matrixelementvector.end());
  std::vector<double>::iterator minelem =
    std::min_element(matrixelementvector.begin(), matrixelementvector.end());

  std::cout << "***********************************" << std::endl
            << "NumMatrixElements     = " << num_mes << std::endl
            << std::scientific
            << "MeanMatrixElemValue   = " << meanelem << " GeV^" << meGeVexponent << std::endl
            << "StdErrMatrixElemValue = " << stdelem/sqrt(num_mes) << " GeV^" << meGeVexponent << std::endl
            << "StdDevMatrixElemValue = " << stdelem << " GeV^" << meGeVexponent << std::endl
            << "MinMatrixElemValue    = " << *minelem << " GeV^" << meGeVexponent << std::endl
            << "MaxMatrixElemValue    = " << *maxelem << " GeV^" << meGeVexponent << std::endl;
}
