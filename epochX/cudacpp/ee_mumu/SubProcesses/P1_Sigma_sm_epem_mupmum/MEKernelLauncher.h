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

    // Compute the output MEs from the input momenta
    void computeMEs() const;

    // Get the buffer for the input momenta: AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
    // (host-allocated in C++ and device-allocated in CUDA)
    fptype* momenta() const { return m_momenta; }

    // Get the buffer for the output MEs: ARRAY[nevt], final |M|^2 averaged over helicities
    // (host-allocated in C++ and device-allocated in CUDA)
    const fptype* MEs() const { return m_MEs; }

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

    // The buffer for the input momenta: AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
    // (host-allocated in C++ and device-allocated in CUDA)
    fptype* m_momenta;

    // The buffer for the output MEs: ARRAY[nevt], final |M|^2 averaged over helicities
    // (host-allocated in C++ and device-allocated in CUDA)
    fptype* m_MEs;

  };

}
#endif // MG5_MEKL_Sigma_sm_epem_mupmum_H
