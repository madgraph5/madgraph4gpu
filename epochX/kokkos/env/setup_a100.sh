module load cuda/11.6.2 gcc/11.1.0 cmake

export KOKKOS_HOME=/gpfs/jlse-fs0/projects/AtlasESP/kokkos/install_a100_gcc-11.1.0_cuda-11.6.2_noUVM
#export KOKKOS_HOME=/home/jchilders/git/kokkos/install_v100/
#export KOKKOS_HOME=/gpfs/jlse-fs0/projects/AtlasESP/kokkos/build_v100
export CMAKE_PREFIX_PATH=$KOKKOS_HOME/lib64/cmake/Kokkos
export LD_LIBRARY_PATH=$KOKKOS_HOME/lib64:$LD_LIBRARY_PATH
export PATH=$KOKKOS_HOME/lib64:$PATH
