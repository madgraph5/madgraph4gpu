#ifndef COLORAMPS_H
#define COLORAMPS_H 1

namespace mgOnGpu
{

  __device__ constexpr bool icolamp[5][4] = {
    { false, false, false, true },
    { false, false, false, true },
    { true, false, false, false },
    { true, false, false, false },
    { true, false, false, true }
  };

}
#endif // COLORAMPS_H
