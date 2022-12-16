#ifndef COLORAMPS_H
#define COLORAMPS_H 1

namespace mgOnGpu
{

  __device__ constexpr bool icolamp[%(nb_channel)s][%(nb_color)s] = { // FIXME: assume process.nprocesses == 1 for the moment
%(is_LC)s
  };

}
#endif // COLORAMPS_H
