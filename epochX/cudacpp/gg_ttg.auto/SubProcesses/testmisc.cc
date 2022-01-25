//#include "mgOnGpuConfig.h"

#include <gtest/gtest.h>
#include "epoch_process_id.h"

#ifdef __CUDACC__
#define TESTID( s ) s##_GPU_MISC
#else
#define TESTID( s ) s##_CPU_MISC
#endif

#define XTESTID( s ) TESTID( s )

TEST( XTESTID( MG_EPOCH_PROCESS_ID ), testmisc )
{
  EXPECT_TRUE( true );
}
