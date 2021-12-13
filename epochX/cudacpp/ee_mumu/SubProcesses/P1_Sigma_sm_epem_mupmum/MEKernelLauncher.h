//M==========================================================================
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
    enum UseHstBuffers { UseNone=0, UseHstMomenta=1, UseHstMEs=2, UseBoth=3 };
    MEKernelLauncher( int ngpublocks, int ngputhreads, UseHstBuffers useHstBuffers = UseNone );
#else
    MEKernelLauncher( int nevt );
#endif

    // Destructor
    // Deallocates the input and output buffers
    ~MEKernelLauncher();

    // === MOMENTA ===

    // Get the host buffer[nevt*npar*4] for the input momenta (CPU or GPU)
    // [NB on GPU, this is a nullptr unless UseHstMomenta or UseBoth are specified]
    fptype* hstMomenta() const { return m_hstMomenta; }

#ifdef __CUDACC__
    // Get the device buffer[nevt*npar*4] for the input momenta (GPU)
    fptype* devMomenta() const { return m_devMomenta; }

    // Copy the input momenta from the device buffer to the host buffer (GPU)
    // [NB on GPU, this throws unless UseHstMomenta or UseBoth are specified]
    void copyDevMomentaToHstMomenta() const;

    // Copy the input momenta from the host buffer to the device buffer (GPU)
    // [NB on GPU, this throws unless UseHstMomenta or UseBoth are specified]
    void copyHstMomentaToDevMomenta() const;
#endif

    // === MATRIX ELEMENTS (final |M|^2 averaged over helicities) ===

#ifndef __CUDACC__
    // Compute the output host MEs from the input host momenta (CPU)
    void computeHstMEsFromHstMomenta() const;
#else
    // Compute the output device MEs from the input device momenta (GPU)
    void computeDevMEsFromDevMomenta() const;
#endif

    // Get the host buffer[nevt] for the output MEs
    // [NB on GPU, this is a nullptr unless UseHstMomenta or UseBoth are specified]
    const fptype* hstMEs() const { return m_hstMEs; }

#ifdef __CUDACC__
    // Get the device buffer[nevt] for the output MEs
    const fptype* devMEs() const { return m_devMEs; }

    // Copy the output MEs from the device buffer to the host buffer
    // [NB this throws unless UseHstMEs or UseBoth are specified]
    void copyDevMEsToHstMEs() const;
#endif

    // === HELICITIES ===

#ifndef __CUDACC__
    // Compute good helicities from the input host momenta (CPU)
    void computeGoodHelFromHstMomenta();
#else
    // Compute good helicities from the input device momenta (GPU)
    void computeGoodHelFromDevMomenta();
#endif

    // Get the host buffer[ncomb] for the helicity mask
    const bool* hstIsGoodHel();

    /*
    // [NOT USED YET - MAY BE USEFUL FOR MULTI THREADING]
    // Set the helicity mask from user input (instead of computing it)
    // [NB throws if good helicities have already been set or computed]
    void setGoodHel( const bool* isGoodHel );
    */

    // Get the number of elements in the momenta buffer
    int nMomenta() const { return np4 * npar * m_nevt; }

    // Get the number of elements in the MEs buffer
    int nMEs() const { return m_nevt; }

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
    static constexpr int ncomb = mgOnGpu::ncomb;
    //static constexpr int wrows = 6; // mgOnGpu::nw6;

    // Alignment of the input and output buffers
    // (strictly needed only for the momenta AOSOA to allow using reinterpret_cast with SIMD vectorized C++ code)
#ifndef __CUDACC__
    static constexpr int cppAlign = 64; // alignment requirement for SIMD vectorization (64-byte i.e. 512-bit)
#endif

  private:

#ifdef __CUDACC__
    // The number of gpu blocks
    const int m_ngpublocks;

    // The number of gpu threads
    const int m_ngputhreads;
#endif

    // The number of events
    const int m_nevt;

    // The host buffer[nevt*npar*4] for the input momenta (CPU or GPU)
    fptype* m_hstMomenta;

#ifdef __CUDACC__
    // The device buffer[nevt*npar*4] for the input momenta (GPU)
    fptype* m_devMomenta;
#endif

    // The host buffer[nevt] for the input MEs (CPU or GPU)
    fptype* m_hstMEs;

#ifdef __CUDACC__
    // The device buffer[nevt] for the input MEs (GPU)
    fptype* m_devMEs;
#endif

    // The host buffer[ncomb] for the helicity mask
    bool* m_hstIsGoodHel;

  };

}
#endif // MG5_MEKL_Sigma_sm_epem_mupmum_H
