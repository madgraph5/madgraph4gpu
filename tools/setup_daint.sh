module purge

module use /apps/daint/UES/hackaton/modules/all
module load CUDAcore/11.7.1
module load gcc

export FC=gfortran
export MADGRAPH_CUDA_ARCHITECTURE=60
