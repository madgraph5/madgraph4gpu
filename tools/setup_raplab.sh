module purge

module load gcc/9.3.0/cmake
module load cuda/11.6.1

if [ "$USER" == "p6asv6bp" ]; then # AV
  export PATH=~/ccache/usr/bin:$PATH
  export USECCACHE=1
  export CCACHE_DIR=~/ccache/CCACHE_DIR
  ###export PATH=~/gfortran-9/usr/bin:$PATH # NOT OK, liblto_plugin.so is missing
fi

export CXX=$(which g++)
export NVCC=$(which nvcc)
echo "Using CXX=${CXX}"
echo "Using NVCC=${NVCC}"

export MADGRAPH_CUDA_ARCHITECTURE=80
echo "Using MADGRAPH_CUDA_ARCHITECTURE=${MADGRAPH_CUDA_ARCHITECTURE}"

if gfortran-9 --version >& /dev/null; then
  export FC=$(which gfortran-9)
  ###export LIBRARY_PATH=$(dirname $(${FC} --print-file-name libgfortran.so)):$LIBRARY_PATH
  echo "Using FC=${FC}"
  echo "Using LIBRARY_PATH=${LIBRARY_PATH}"
else
  echo "ERROR! gfortran-9 is not installed on this node: please use another node"
fi

