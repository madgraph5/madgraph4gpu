#ifndef MGONGPUCONFIG_H
#define MGONGPUCONFIG_H 1

// Memory layout for momenta
#define MGONGPU_LAYOUT_ASA 1 // default
//#define MGONGPU_LAYOUT_SOA 1
//#define MGONGPU_LAYOUT_AOS 1

// Curand random number generation
#define MGONGPU_CURAND_ONDEVICE 1 // default
//#define MGONGPU_CURAND_ONHOST 1

namespace mgOnGpu
{
  // Number of Events Per Page in the AOSOA (ASA) structure
  const int nepp = 32; // choose 32, like the number of threads in a warp
}

#endif // MGONGPUCONFIG_H
