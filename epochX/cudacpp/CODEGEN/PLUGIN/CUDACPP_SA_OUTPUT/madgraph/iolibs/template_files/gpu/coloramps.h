#ifndef COLORAMPS_H
#define COLORAMPS_H 1

namespace mgOnGpu
{

  __device__ constexpr bool icolamp[%(nb_channel)s][%(nb_color)s] = {
%(is_LC)s
  };

}
#endif // COLORAMPS_H
