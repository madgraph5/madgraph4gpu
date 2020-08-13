#ifndef MGONGPUCONFIG_H
#define MGONGPUCONFIG_H 1

// ** NB1 Throughputs (e.g. 5.00E8) are events/sec for "./gcheck.exe -p 16384 32 12"
// ** NB2 Baseline on b7g47n0002 fluctuates (depends on load?): typically, either ~5.00E8 or ~5.50E8

// Memory layout for momenta (CHOOSE ONLY ONE)
#define MGONGPU_LAYOUT_ASA 1 // default (~5.00E8) 
//#define MGONGPU_LAYOUT_SOA 1 // 7% slower (4.65E8 against 5.00E8)
//#define MGONGPU_LAYOUT_AOS 1 // 4% slower (4.80E8 against 5.00E8)

// Curand random number generation (CHOOSE ONLY ONE)
#define MGONGPU_CURAND_ONDEVICE 1 // default
//#define MGONGPU_CURAND_ONHOST 1

// Use global memory or shared memory for wavefunctions (CHOOSE ONLY ONE)
//#define MGONGPU_WFMEM_LOCAL 1 // default (~5.00E8)
#define MGONGPU_WFMEM_GLOBAL 1 // 30% slower, limited to 256*32 threads (1.18E8 against 1.78E8 for "-p 256 32 12")
//#define MGONGPU_WFMEM_SHARED 1 // 30% slower, limited to 32 threads/block (~3.5E8 against 5.0E8)

// Floating point precision (CHOOSE ONLY ONE)
#define MGONGPU_FPTYPE_DOUBLE 1 // default (~5.00E8)
//#define MGONGPU_FPTYPE_FLOAT 1 // 2.3x faster (~1.14E9 against 5.00E8)

// Complex type in cuda: thrust or cucomplex (CHOOSE ONLY ONE)
#ifdef __CUDACC__
#define MGONGPU_CXTYPE_THRUST 1 // default (~5.00E8)
//#define MGONGPU_CXTYPE_CUCOMPLEX 1 // ~5% slower (4.75E8 against 5.00E8)
#endif

namespace mgOnGpu
{
  // --- Physics process-specific constants that are best declared at compile time

  const int np4 = 4; // the dimension of 4-momenta (E,px,py,pz)
  const int nw6 = 6; // dimension of each wavefunction (see KEK 91-11)

  const int npari = 2; // #particles in the initial state (incoming): e+ e-
  const int nparf = 2; // #particles in the final state (outgoing): mu+ mu-
  const int npar = npari + nparf; // #particles in total (external): e+ e- -> mu+ mu-

  const int nwf = 5; // #wavefunctions: npar (4 external) + 1 (internal, reused for gamma and Z)

  // --- Platform-specific software implementation details

  // Maximum number of blocks per grid
  // ** NB Some arrays of pointers will be allocated statically to fit all these blocks
  // ** (the actual memory for each block will then be allocated dynamically only for existing blocks)
  //const int nbpgMAX = 2048;

  // Maximum number of threads per block
  const int ntpbMAX = 256;
  //const int ntpbMAX = 8; // FOR DEBUGGING!

  // Number of Events Per Page in the AOSOA (ASA) structure
  // ** TODO: this will disappear and become dynamic...
  const int nepp = 32; // choose 32, i.e. the number of threads in a warp
  //const int nepp = ntpbMAX; // choose 256, i.e. the max number of threads in a block
  //const int nepp = 4; // FOR DEBUGGING!

}

#endif // MGONGPUCONFIG_H
