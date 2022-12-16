#ifndef COLORAMPS_H
#define COLORAMPS_H 1

namespace mgOnGpu
{

  __device__ constexpr bool icolamp[3][2] = { // FIXME: assume process.nprocesses == 1 for the moment
    { true, true },
    { true, false },
    { false, true }

  };

}
#endif // COLORAMPS_H
