#!/bin/bash

module load gcc/9.3.0/cmake
module load cuda/11.6.1

export MADGRAPH_CUDA_ARCHITECTURE=80

if [ -f /usr/bin/gfortran-8 ]; then
  export FC=/usr/bin/gfortran-8
  export LIBRARY_PATH=/usr/lib/gcc/x86_64-linux-gnu/8/
elif [ -f /usr/bin/gfortran-9 ]; then
  export FC=/usr/bin/gfortran-9
  export LIBRARY_PATH=/usr/lib/gcc/x86_64-linux-gnu/9/
fi

CXX=`which g++`
NVCC=`which nvcc`

echo "Using compilers:"
echo "FC=${FC}"
echo "CXX=${CXX}"
echo "NVCC=${NVCC}"
echo "MADGRAPH_CUDA_ARCHITECTURE=${MADGRAPH_CUDA_ARCHITECTURE}"

export LD_LIBRARY_PATH=../../lib:${LD_LIBRARY_PATH}
