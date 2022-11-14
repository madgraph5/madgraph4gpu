#include "RandomNumberKernels.h"

#include "CommonRandomNumbers.h"
#include "MemoryBuffers.h"

#include <cassert>

namespace mg5amcGpu
{

  //--------------------------------------------------------------------------

  CommonRandomNumberKernel::CommonRandomNumberKernel( BufferRandomNumbers& rnarray )
    : RandomNumberKernelBase( rnarray )
    , m_seed( 20211220 )
  {
    if ( m_rnarray.isOnDevice() )
      throw std::runtime_error( "CommonRandomNumberKernel on host with a device random number array" );
  }

  //--------------------------------------------------------------------------

  void CommonRandomNumberKernel::generateRnarray()
  {
    std::vector<double> rnd = CommonRandomNumbers::generate<double>( m_rnarray.size(), m_seed ); // NB: HARDCODED DOUBLE!
    std::copy( rnd.begin(), rnd.end(), m_rnarray.data() ); // NB: this may imply a conversion from double to float
  }

  //--------------------------------------------------------------------------

}
