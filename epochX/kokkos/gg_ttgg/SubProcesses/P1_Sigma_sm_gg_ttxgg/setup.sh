
export GCC_PATH=/soft/compilers/gcc/9.2.0/linux-rhel7-x86_64
export CMAKE_PATH=/soft/buildtools/cmake/3.22.1
export PATH=/usr/local/bin:/usr/bin:/bin:/usr/lpp/mmfs/bin:/home/jchilders/bin
export LD_LIBRARY_PATH=

CUDA_VERSION_MAJOR=11
CUDA_VERSION_MID=6
CUDA_VERSION_MINOR=2
CUDA_VERSION_BUILD=510.47.03

export CUDA_VERSION=${CUDA_VERSION_MAJOR}.${CUDA_VERSION_MID}.${CUDA_VERSION_MINOR}_${CUDA_VERSION_BUILD}
export CUDA_BASE_PATH=/gpfs/jlse-fs0/projects/AtlasESP/cuda/
export CUDA_PATH=$CUDA_BASE_PATH/cuda_${CUDA_VERSION}_linux

export PATH=$GCC_PATH/bin:$CMAKE_PATH/bin:$CUDA_PATH/bin:$PATH
export LD_LIBRARY_PATH=$GCC_PATH/lib:$CMAKE_LIB/lib:$CUDA_PATH/lib64:$LD_LIBRARY_PATH


