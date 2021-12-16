#include "RandomNumberKernel.h"

#include "checkCuda.h"

#ifdef __CUDACC__
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  void RandomNumberKernelBase::copyHstRnarrayToDevRnArray() 
  {
    const int nRnarray = np4 * nparf * m_nevt;
    const int nbytesRnarray = nRnarray * sizeof(fptype);
    if ( !m_hstRnarray )
      throw std::logic_error( "Cannot copyHsrRnarrayToDevRnArray unless random numbers are generated on the host" );
    checkCuda( cudaMemcpy( m_devRnarray, m_hstRnarray, nbytesRnarray, cudaMemcpyHostToDevice ) );
  }
#endif

  //--------------------------------------------------------------------------

}
