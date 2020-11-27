#include <iostream>
#include <stdio.h>
#include <vector>
#include <list>
#include <utility>
#include <algorithm>
#include <iomanip>


class Properties {

private:
  typedef std::vector<std::pair<std::string, std::string>> PTYPE;
  std::list<PTYPE> allprops;
  PTYPE* theseprops = nullptr;

public:

  Properties& add(const std::string& k, char* v) {
      theseprops->emplace_back(std::make_pair(k, std::string(v)));
      return *this;
  }

  template<typename T>
  Properties& add(const std::string& k, T v) {
    theseprops->emplace_back(std::make_pair(k, std::to_string(v)));
    return *this;
  }

  template<typename T>
  Properties& add(const std::string& k, T v[], int s) {
    std::string temps;
    for (int i = 0; i < s; ++i) {
      temps += std::to_string(v[i]);
      if (i < (s-1)) temps += ", ";
    }
    theseprops->emplace_back(std::make_pair(k, "[" + temps + "]"));
    return *this;
  }

  void fill() {
    int nDevices;
    cudaGetDeviceCount(&nDevices);
    cudaDeviceProp prop;
    for (int i = 0; i < nDevices; i++) {
      allprops.emplace_back(PTYPE());
      theseprops = &allprops.back();
      cudaGetDeviceProperties(&prop, i);
      (*this)
      .add("ECCEnabled", prop.ECCEnabled)
      .add("accessPolicyMaxWindowSize", prop.accessPolicyMaxWindowSize)
      .add("asyncEngineCount", prop.asyncEngineCount)
      .add("canMapHostMemory", prop.canMapHostMemory)
      .add("canUseHostPointerForRegisteredMem",
           prop.canUseHostPointerForRegisteredMem)
      .add("clockRate", prop.clockRate)
      .add("computeMode", prop.computeMode)
      .add("computePreemptionSupported", prop.computePreemptionSupported)
      .add("concurrentKernels", prop.concurrentKernels)
      .add("concurrentManagedAccess", prop.concurrentManagedAccess)
      .add("cooperativeLaunch", prop.cooperativeLaunch)
      .add("cooperativeMultiDeviceLaunch", prop.cooperativeMultiDeviceLaunch)
      .add("deviceOverlap", prop.deviceOverlap)
      .add("directManagedMemAccessFromHost", prop.directManagedMemAccessFromHost)
      .add("globalL1CacheSupported", prop.globalL1CacheSupported)
      .add("hostNativeAtomicSupported", prop.hostNativeAtomicSupported)
      .add("integrated", prop.integrated)
      .add("isMultiGpuBoard", prop.isMultiGpuBoard)
      .add("kernelExecTimeoutEnabled", prop.kernelExecTimeoutEnabled)
      .add("l2CacheSize", prop.l2CacheSize)
      .add("localL1CacheSupported", prop.localL1CacheSupported)
      .add("luid[8]", prop.luid)
      .add("luidDeviceNodeMask", prop.luidDeviceNodeMask)
      .add("major", prop.major)
      .add("managedMemory", prop.managedMemory)
      .add("maxBlocksPerMultiProcessor", prop.maxBlocksPerMultiProcessor)
      .add("maxGridSize[3]", prop.maxGridSize, 3)
      .add("maxSurface1D", prop.maxSurface1D)
      .add("maxSurface1DLayered[2]", prop.maxSurface1DLayered, 2)
      .add("maxSurface2D[2]", prop.maxSurface2D, 2)
      .add("maxSurface2DLayered[3]", prop.maxSurface2DLayered, 3)
      .add("maxSurface3D[3]", prop.maxSurface3D, 3)
      .add("maxSurfaceCubemap", prop.maxSurfaceCubemap)
      .add("maxSurfaceCubemapLayered[2]", prop.maxSurfaceCubemapLayered, 2)
      .add("maxTexture1D", prop.maxTexture1D)
      .add("maxTexture1DLayered[2]", prop.maxTexture1DLayered, 2)
      .add("maxTexture1DLinear", prop.maxTexture1DLinear)
      .add("maxTexture1DMipmap", prop.maxTexture1DMipmap)
      .add("maxTexture2D[2]", prop.maxTexture2D, 2)
      .add("maxTexture2DGather[2]", prop.maxTexture2DGather, 2)
      .add("maxTexture2DLayered[3]", prop.maxTexture2DLayered, 3)
      .add("maxTexture2DLinear[3]", prop.maxTexture2DLinear, 3)
      .add("maxTexture2DMipmap[2]", prop.maxTexture2DMipmap, 3)
      .add("maxTexture3D[3]", prop.maxTexture3D, 3)
      .add("maxTexture3DAlt[3]", prop.maxTexture3DAlt, 3)
      .add("maxTextureCubemap", prop.maxTextureCubemap)
      .add("maxTextureCubemapLayered[2]", prop.maxTextureCubemapLayered, 2)
      .add("maxThreadsDim[3]", prop.maxThreadsDim, 3)
      .add("maxThreadsPerBlock", prop.maxThreadsPerBlock)
      .add("maxThreadsPerMultiProcessor", prop.maxThreadsPerMultiProcessor)
      .add("memPitch", prop.memPitch)
      .add("memoryBusWidth", prop.memoryBusWidth)
      .add("memoryClockRate", prop.memoryClockRate)
      .add("minor", prop.minor)
      .add("multiGpuBoardGroupID", prop.multiGpuBoardGroupID)
      .add("multiProcessorCount", prop.multiProcessorCount)
      .add("name[256]", prop.name)
      .add("pageableMemoryAccess", prop.pageableMemoryAccess)
      .add("pageableMemoryAccessUsesHostPageTables",
           prop.pageableMemoryAccessUsesHostPageTables)
      .add("pciBusID", prop.pciBusID)
      .add("pciDeviceID", prop.pciDeviceID)
      .add("pciDomainID", prop.pciDomainID)
      .add("persistingL2CacheMaxSize", prop.persistingL2CacheMaxSize)
      .add("regsPerBlock", prop.regsPerBlock)
      .add("regsPerMultiprocessor", prop.regsPerMultiprocessor)
      .add("reservedSharedMemPerBlock", prop.reservedSharedMemPerBlock)
      .add("sharedMemPerBlock", prop.sharedMemPerBlock)
      .add("sharedMemPerBlockOptin", prop.sharedMemPerBlockOptin)
      .add("sharedMemPerMultiprocessor", prop.sharedMemPerMultiprocessor)
      .add("singleToDoublePrecisionPerfRatio",
           prop.singleToDoublePrecisionPerfRatio)
      .add("streamPrioritiesSupported", prop.streamPrioritiesSupported)
      .add("surfaceAlignment", prop.surfaceAlignment)
      .add("tccDriver", prop.tccDriver)
      .add("textureAlignment", prop.textureAlignment)
      .add("texturePitchAlignment", prop.texturePitchAlignment)
      .add("totalConstMem", prop.totalConstMem)
      .add("totalGlobalMem", prop.totalGlobalMem)
      .add("unifiedAddressing", prop.unifiedAddressing)
      .add("uuid", (unsigned char*)prop.uuid.bytes, 16)
      .add("warpSize", prop.warpSize);
    }
  }

  void print() {
    int did = 0;
    size_t maxwidth = 0;

    for (auto i : allprops.front())
      maxwidth = std::max(i.first.length(), maxwidth);

    for (auto p : allprops ) {
      std::cout << "DeviceID: " << did << std::endl;
      for (auto i : p) {
        std::cout << "  " << std::left << std::setw(maxwidth+1) << i.first << i.second << std::endl;
      }
      ++did;
    }
  }

};

/* calculate cores per streaming multiprocessor 
https://stackoverflow.com/questions/32530604/how-can-i-get-number-of-cores-in-cuda-device
#include "cuda_runtime_api.h"
// you must first call the cudaGetDeviceProperties() function, then pass
// the devProp structure returned to this function:
int getSPcores(cudaDeviceProp devProp)
{
    int cores = 0;
    int mp = devProp.multiProcessorCount;
    switch (devProp.major){
     case 2: // Fermi
      if (devProp.minor == 1) cores = mp * 48;
      else cores = mp * 32;
      break;
     case 3: // Kepler
      cores = mp * 192;
      break;
     case 5: // Maxwell
      cores = mp * 128;
      break;
     case 6: // Pascal
      if ((devProp.minor == 1) || (devProp.minor == 2)) cores = mp * 128;
      else if (devProp.minor == 0) cores = mp * 64;
      else printf("Unknown device type\n");
      break;
     case 7: // Volta and Turing
      if ((devProp.minor == 0) || (devProp.minor == 5)) cores = mp * 64;
      else printf("Unknown device type\n");
      break;
     case 8: // Ampere
      if (devProp.minor == 0) cores = mp * 64;
      else if (devProp.minor == 6) cores = mp * 128;
      else printf("Unknown device type\n");
      break;
     default:
      printf("Unknown device type\n");
      break;
      }
    return cores;
}
*/


int main() {
  Properties p;
  p.fill();
  p.print();

  exit(0);
}
