module purge

module use /apps/daint/UES/hackaton/modules/all
module load gcc
module load CUDAcore/11.7.1

export CXX=$(which g++)
export NVCC=$(which nvcc)
echo "Using CXX=${CXX}"
echo "Using NVCC=${NVCC}"

export MADGRAPH_CUDA_ARCHITECTURE=60
echo "Using MADGRAPH_CUDA_ARCHITECTURE=${MADGRAPH_CUDA_ARCHITECTURE}"

export FC=$(which gfortran)
echo "Using FC=${FC}"

if [ "$USER" == "hck03" ]; then # AV
  export PATH=~/ccache/bin:$PATH
  export USECCACHE=1
  export CCACHE_DIR=~/ccache/CCACHE_DIR
fi

