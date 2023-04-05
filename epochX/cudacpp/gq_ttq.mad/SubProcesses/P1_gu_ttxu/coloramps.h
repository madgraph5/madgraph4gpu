#ifndef COLORAMPS_H
#define COLORAMPS_H 1

namespace mgOnGpu
{

  __device__ constexpr bool icolamp[5][4] = { // FIXME: assume process.nprocesses == 1 for the moment
    { true, false, false, false },
    { true, false, false, false },
    { false, false, false, true },
    { false, false, false, true },
    { true, false, false, true }
  };

}
#endif // COLORAMPS_H
