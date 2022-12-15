#ifndef COLORAMPS_H
#define COLORAMPS_H 1

namespace mgOnGpu
{

  constexpr int ncolor = 2;

  constexpr int nchannel = 3;

  __device__
  constexpr bool icolamp[nchannel][ncolor] = { // FIXME: assume process.nprocesses == 1 for the moment
    { true, true },
    { true, false },
    { false, true }
  };

}
#endif // COLORAMPS_H
