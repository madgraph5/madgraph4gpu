# Building Kokkos

To use the Kokkos plugin, first a Kokkos installation for the specific hard is required.

You can find details here:
[Kokkos Compilation](https://github.com/kokkos/kokkos/wiki/Compiling)

# Building for NVidia A100

```bash
git clone https://github.com/kokkos/kokkos.git
cd kokkos
git checkout 3.5.00
export KOKKOS_BASE_PATH=$PWD
# Need to have CUDA libraries in PATH/LD_LIBRARY_PATH

mkdir build_a100
cd build_a100
cmake .. -DCMAKE_INSTALL_PREFIX=$KOKKOS_BASE_PATH/install_v100 -DKokkos_ENABLE_CUDA=On -DKokkos_ENABLE_CUDA_LAMBDA=On -DKokkos_ENABLE_CUDA_UVM=On -DKokkos_ARCH_AMPERE80=On -DKokkos_ENABLE_OPENMP=On

make -j install
```

This can be used to build an A100 compatible executable.

After using the plugin to generate a process, go to `<process>/SubProcesses/P1.../` then run `KOKKOSPATH_CUDA=$KOKKOS_BASE_PATH/install_v100 CUDA_ARCH_NUM=80 make cuda`. This should build the `ccheck.exe`.
