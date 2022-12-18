#ifndef COLORAMPS_H
#define COLORAMPS_H 1

namespace mgOnGpu
{

  __device__ constexpr bool icolamp[2][1] = { // FIXME: assume process.nprocesses == 1 for the moment
    { true },
    { true }
  };

}
#endif // COLORAMPS_H
