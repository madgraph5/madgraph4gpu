#ifndef MGONGPUCONFIG_H
#define MGONGPUCONFIG_H 1

#define MGONGPU_LAYOUT_ASA 1
//#define MGONGPU_LAYOUT_SOA 1
//#define MGONGPU_LAYOUT_AOS 1

#if defined MGONGPU_LAYOUT_ASA
namespace mgOnGpu
{
  // Number of Events Per Page in the AOSOA (ASA) structure
  const int nepp = 32; // choose 32, like the number of threads in a warp
}
#endif

#endif // MGONGPUCONFIG_H
