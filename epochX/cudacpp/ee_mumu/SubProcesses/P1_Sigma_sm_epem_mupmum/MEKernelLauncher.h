//==========================================================================
// This file has been automatically generated for CUDA/C++ standalone by
// MadGraph5_aMC@NLO v. 2.9.5, 2021-08-22
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef MG5_MEKL_Sigma_sm_epem_mupmum_H
#define MG5_MEKL_Sigma_sm_epem_mupmum_H

#include "mgOnGpuConfig.h"
//#include "mgOnGpuTypes.h"

#ifdef __CUDACC__
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{

  //==========================================================================
  // A class for calculating the matrix elements for
  // Process: e+ e- > mu+ mu- WEIGHTED<=4 @1
  //--------------------------------------------------------------------------

  class MEKernelLauncher
  {
  public:

    // Constructor (from command line arguments)
    // Allocates the input and output buffers for the given number of events
#ifdef __CUDACC__
    MEKernelLauncher( int ngpublocks, int ngputhreads );
#else
    MEKernelLauncher( int nevt );
#endif

    // Destructor
    // Deallocates the input and output buffers
    ~MEKernelLauncher();

#ifdef __CUDACC__
    // Compute the output device MEs from the input device momenta
    void computeDevMEs() const;

    // Get the device buffer for the input momenta: AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
    fptype* devMomenta() const { return m_devMomenta; }

    // Get the device buffer for the output MEs: ARRAY[nevt], final |M|^2 averaged over helicities
    const fptype* devMEs() const { return m_devMEs; }
#else
    // Compute the output host MEs from the input host momenta
    void computeHstMEs() const;

    // Get the host buffer for the input momenta: AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
    fptype* hstMomenta() const { return m_hstMomenta; }

    // Get the host buffer for the output MEs: ARRAY[nevt], final |M|^2 averaged over helicities
    const fptype* hstMEs() const { return m_hstMEs; }
#endif

    // Get the number of bytes in the momenta array
    int nbytesMomenta() const { return np4 * npar * m_nevt * sizeof(fptype); }
    
    // Get the number of bytes in the MEs array
    int nbytesMEs() const { return m_nevt * sizeof(fptype); }

  public:

    // Hardcoded parameters (temporarely set them from mgOnGpu; eventually define them only here?)
    static constexpr int npar = mgOnGpu::npar;
    static constexpr int np4 = 4;
    static constexpr int neppM = mgOnGpu::neppM;

    // Hardcoded parameters for this process (constant class variables)
    //static constexpr int ninitial = 2; // mgOnGpu::npari;
    //static constexpr int nexternal = 4; // mgOnGpu::npar (nexternal was nioparticles)
    //static constexpr int nprocesses = 1; // FIXME: assume process.nprocesses == 1
    //static constexpr int nwavefuncs = 6; // mgOnGpu::nwf
    //static constexpr int namplitudes = 2;
    //static constexpr int ncomb = 16; // mgOnGpu::ncomb
    //static constexpr int wrows = 6; // mgOnGpu::nw6;

    // Alignment of the input and output buffers
    // (strictly needed only for the momenta AOSOA to allow using reinterpret_cast with SIMD vectorized C++ code)
#ifndef __CUDACC__
    static constexpr int cppAlign = 64; // alignment requirement for SIMD vectorization (64-byte i.e. 512-bit)
#endif

  private:

    // The number of events
    const int m_nevt;

#ifdef __CUDACC__
    // The number of gpu blocks
    const int m_ngpublocks;

    // The number of gpu threads
    const int m_ngputhreads;
#endif

#ifdef __CUDACC__
    // The device buffer for the input momenta: AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
    fptype* m_devMomenta;

    // The device buffer for the output MEs: ARRAY[nevt], final |M|^2 averaged over helicities
    fptype* m_devMEs;
#else
    // The host buffer for the input momenta: AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
    fptype* m_hstMomenta;

    // The host buffer for the output MEs: ARRAY[nevt], final |M|^2 averaged over helicities
    fptype* m_hstMEs;
#endif

  };

}
#endif // MG5_MEKL_Sigma_sm_epem_mupmum_H
