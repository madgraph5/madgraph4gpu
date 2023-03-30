# Madgraph 4 GPU

[![C/C++ CI](https://github.com/Jooorgen/madgraph4gpu/actions/workflows/c-cpp.yml/badge.svg)](https://github.com/Jooorgen/madgraph4gpu/actions/workflows/c-cpp.yml) [![SYCL CI](https://github.com/Jooorgen/madgraph4gpu/actions/workflows/sycl.yml/badge.svg)](https://github.com/Jooorgen/madgraph4gpu/actions/workflows/sycl.yml) [![CUDA Profiler](https://github.com/Jooorgen/madgraph4gpu/actions/workflows/cudaProfiler.yml/badge.svg)](https://github.com/Jooorgen/madgraph4gpu/actions/workflows/cudaProfiler.yml) [![SYCL Profiler](https://github.com/Jooorgen/madgraph4gpu/actions/workflows/syclProfiler.yml/badge.svg)](https://github.com/Jooorgen/madgraph4gpu/actions/workflows/syclProfiler.yml)

This repository contains code developed in the context of porting the [MadGraph5_aMC@NLO](https://cp3.irmp.ucl.ac.be/projects/madgraph/) event generator software onto GPU platforms and vector instructions on CPUs. MadGraph5_aMC@NLO is able to generate code for various physics processes in different programming languages (Fortran, C, C++). The code generated in this repository in "epochX" of the MadGraph5_aMC@NLO generator allows to also produce source code for those physics processes to run on GPU and CPU platforms. 


