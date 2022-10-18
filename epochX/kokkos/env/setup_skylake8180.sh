module load gcc/12.1.0 cmake

export KOKKOS_HOME=/gpfs/jlse-fs0/projects/AtlasESP/kokkos/install_skylake_gcc-12.1.0

export CMAKE_PREFIX_PATH=$KOKKOS_HOME/lib64/cmake/Kokkos
export LD_LIBRARY_PATH=$KOKKOS_HOME/lib64:$LD_LIBRARY_PATH
export PATH=$KOKKOS_HOME/lib64:$PATH
