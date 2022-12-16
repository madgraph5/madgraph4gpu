#ifndef COLORAMPS_H
#define COLORAMPS_H 1

namespace mgOnGpu
{

  __device__ constexpr bool icolamp[15][6] = { // FIXME: assume process.nprocesses == 1 for the moment
    { true, false, true, false, true, true },
    { false, false, false, false, true, true },
    { true, false, true, false, false, false },
    { false, true, false, false, false, false },
    { true, true, false, false, false, false },
    { true, false, false, false, false, false },
    { false, false, false, true, false, false },
    { false, false, false, true, false, true },
    { false, false, false, false, false, true },
    { false, false, true, true, false, false },
    { false, true, false, false, true, false },
    { false, true, true, true, true, false },
    { false, false, true, false, false, false },
    { false, false, false, false, true, false },
    { true, true, false, true, false, true }

  };

}
#endif // COLORAMPS_H
