# Building the code
## Setting up the compiler
A [SYCL compatible compiler](https://github.com/intel/llvm/blob/sycl/sycl/doc/GetStartedGuide.md)
is required to build the code. Start by setting the compiler environment
variable `export CXX=/path/to/sycl-compiler` to a SYCL compatible such as
`clang++` or `dpcpp`.
Next set the `SYCLFLAGS` environment variable `export SYCLFLAGS=<...>` for SYCL
specific flags to send to the compiler ([see table below for a few examples](#sycl-flags)).
This set of flags will set the target device (or devices) and other SYCL related
properties.
The SYCL framework can target multiple devices on the same machine if available.
Optional environment variables are [listed at the end of this section](#additional-build-variables).

### SYCL flags
| Device | `SYCLFLAGS` |
| ---: | :--- |
| NVIDIA V100 | `"-fsycl -fsycl-targets=nvptx64-nvidia-cuda -target-backend '--cuda-gpu-arch=sm_70' -fgpu-rdc --cuda-path=</path/to/cuda-toolkit>"` |
| NVIDIA A100 | `"-fsycl -fsycl-targets=nvptx64-nvidia-cuda -target-backend '--cuda-gpu-arch=sm_80' -fgpu-rdc --cuda-path=</path/to/cuda-toolkit>"` |
| AMD MI50 | `"-fsycl -fsycl-targets=amdgcn-amd-amdhsa -Xsycl-target-backend --offload-arch=gfx906 -fgpu-rdc"` |
| AMD MI100 | `"-fsycl -fsycl-targets=amdgcn-amd-amdhsa -Xsycl-target-backend --offload-arch=gfx908 -fgpu-rdc"` |

## Compilation step
Change directories to the desired subprocess directory within the process
directory and run `make` to complete the build.

## Additional build variables
Additional environment variables which have effects on the build:
- `export CXX=/path/to/sycl-compiler` set path to SYCL compatible compiler
- `export SYCLFLAGS=<...>` set appropriate SYCL flags ([see above](#sycl-flags))
- `export USEBUILDDIR=<0|1>` set `1` to use build directory for executable and library files (`default = 0`)
- `export FPTYPE=<d|f>` set to `d` for double precision floating point numbers or `f` for float precision floating point numbers (`default = d`)
- `export HELINL=<0|1>` set to `1` to inline helicity functions (`default = 0`)
- `export NTPBMAX=<maximum-threads-per-block>` set to maximum threads per block for the target device (`default = 1024`)

# Running `madgraph4gpu`
The `madgraph4gpu` app can be run directly from the subprocess directly by
```
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../../lib ./check.exe -j --json_file path/to/output.json --param_card ../../Cards/param_card.dat --device_id <device id> <blocks> <threads> <iterations> 
```
or from anywhere by
```
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/subprocess-dir/../../lib /path/to/subprocess-dir/check.exe -j --json_file path/to/output.json --param_card ../../Cards/param_card.dat --device_id <device_id> <blocks> <threads> <iterations> 
```
The `device_id` can be determined using the `--device_info` flag without the
`--device_id <device_id>` parameter.
Additional or alternate library paths may be required for certain SYCL compiler
installations or when using build directories (see
[example build and run](#example-build-and-run) for more details.
Use the `--help` flag when running for more usage details.

# Example build and run
Here is a set of instructions to run and build for an NVIDIA V100 gpu using
Intel's `clang++` compiler set up using [these instructions](https://github.com/intel/llvm/blob/sycl/sycl/doc/GetStartedGuide.md#build-dpc-toolchain-with-support-for-nvidia-cuda).
```
#set user specific variables
prefix=/path/to/dir-containing-madgraph4gpu-dir
export CUDA_PATH=/path/to/cuda-toolkit
export MG_PROC=gg_ttgg #process to target
export MG_SUBPROC=P1_Sigma_sm_gg_ttxgg #subprocess directory
export DEVICE_ID=0 #if unknown set at the run step after running LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MG_LIBS $MG_EXE --param_card $MG5AMC_CARD_PATH/param_card.dat --device_info 1024 128 10

#set up compiler and compile options
export SYCLFLAGS="-fsycl -fsycl-targets=nvptx64-nvidia-cuda --cuda-include-ptx=sm_70 --cuda-path=$CUDA_PATH"
export USEBUILDDIR=1
export NTPBMAX=1024
export CXX=$DPCPP_HOME/llvm/build/bin/clang++

export WORKSPACE=$HOME/workspace_mg4gpu
mkdir -p $WORKSPACE/mg4gpu/lib
mkdir -p $WORKSPACE/mg4gpu/bin

export MG4GPU_LIB=$WORKSPACE/mg4gpu/lib
export MG4GPU_BIN=$WORKSPACE/mg4gpu/bin

export MG_PROC_DIR=$prefix/madgraph4gpu/epochX/sycl/$MG_PROC
export MG_SP_DIR=$MG_PROC_DIR/SubProcesses/$MG_SUBPROC

export MG_LIBS_DIR="${MG4GPU_LIB}/build_${MG_PROC}_${MG_SUBPROC}_sycl_v100-cuda-11.4"
export MG_LIBS="$DPCPP_HOME/llvm/build/lib:$MG_LIBS_DIR"

export MG_EXE_DIR="${MG4GPU_BIN}/build_${MG_PROC}_${MG_SUBPROC}_sycl_v100-cuda-11.4"
export MG_EXE="$MG_EXE_DIR/check.exe"
export MG5AMC_CARD_PATH=$MG_PROC_DIR/Cards

#build executable
cd $MG_SP_DIR
make
mv ../../lib/build.d_inl0/ $MG_LIBS_DIR
mv build.d_inl0/ $MG_EXE_DIR

#run executable
cd $WORKSPACE
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MG_LIBS $MG_EXE -j --json_file $WORKSPACE/test_v100_sycl-11.4_${MG_PROC}_${MG_SUBPROC}.json --param_card $MG5AMC_CARD_PATH/param_card.dat --device_id $DEVICE_ID 1024 128 10

#view output
vim $WORKSPACE/test_v100_sycl-11.4_${MG_PROC}_${MG_SUBPROC}.json
```
# Detailed example of using Intel's Docker image to build and run

The first part of this example is based upon the instructions for the Intel SYCL-compatible compiler DPC++:

https://github.com/intel/llvm/blob/sycl/sycl/doc/GetStartedGuide.md

This example uses Singularity instead of Docker since that is the container tool that was available
to the author. It should be fairly simple to translate the Singularity commands to their corresponding
Docker versions for Docker users.

These instructions also assume use of CUDA.

First step is to convert the Docker image provided by Intel to Singularity:

`singularity pull docker://ghcr.io/intel/llvm/ubuntu2004_base`

This command creates Singularity image `ubuntu2004_base_latest.sif`.

In your working directory, fetch the Intel compiler code and the madgraph4gpu code:
```
git clone --config core.autocrlf=false https://github.com/intel/llvm -b sycl
git clone https://github.com/madgraph5/madgraph4gpu.git
```

Then set some environment variables:
```
export DPCPP_HOME=/your/working/directory
export CUDA_PATH=/usr/local/cuda-11.6
```
Replace "/your/working/directory" with the full path of your working directory
and "/usr/local/cuda-11.6" with the installed location of CUDA on your system.

Next, execute the Singularity image to run a shell inside the container:
```
singularity exec --writable-tmpfs -B $DPCPP_HOME:/src,$CUDA_PATH -H $DPCPP_HOME ubuntu2004_base_latest.sif /bin/bash
```
Then, in the container, configure and compile the Intel compiler:
```
CUDA_LIB_PATH=$CUDA_PATH/lib64/stubs CC=gcc CXX=g++ python $DPCPP_HOME/llvm/buildbot/configure.py --cuda --cmake-opt="-DCUDA_TOOLKIT_ROOT_DIR=$CUDA_PATH"

CUDA_LIB_PATH=$CUDA_PATH/lib64/stubs CC=gcc CXX=g++ python $DPCPP_HOME/llvm/buildbot/compile.py
```
If DPC++ successfully compiles, then madgraph4gpu can be compiled next.

Many environment variables need to be set:
```
export MG_PROC=gg_ttgg #process to target
export MG_SUBPROC=P1_Sigma_sm_gg_ttxgg #subprocess directory
export DEVICE_ID=0 #if unknown set at the run step after running LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MG_LIBS $MG_EXE --param_card $MG5AMC_CARD_PATH/param_card.dat --device_info 1024 128 10

#set up compiler and compile options
export SYCLFLAGS="-fsycl -fsycl-targets=nvptx64-nvidia-cuda --cuda-include-ptx=sm_75 --cuda-path=$CUDA_PATH"
export USEBUILDDIR=1
export NTPBMAX=1024
export CXX=$DPCPP_HOME/llvm/build/bin/clang++

export WORKSPACE=$DPCPP_HOME/workspace_mg4gpu
mkdir -p $WORKSPACE/mg4gpu/lib
mkdir -p $WORKSPACE/mg4gpu/bin

export MG4GPU_LIB=$WORKSPACE/mg4gpu/lib
export MG4GPU_BIN=$WORKSPACE/mg4gpu/bin

export MG_PROC_DIR=$DPCPP_HOME/madgraph4gpu/epochX/sycl/$MG_PROC
export MG_SP_DIR=$MG_PROC_DIR/SubProcesses/$MG_SUBPROC

export MG_LIBS_DIR="${MG4GPU_LIB}/build_${MG_PROC}_${MG_SUBPROC}_sycl_v100-cuda-11.6"
export MG_LIBS="$DPCPP_HOME/llvm/build/lib:$MG_LIBS_DIR"

export MG_EXE_DIR="${MG4GPU_BIN}/build_${MG_PROC}_${MG_SUBPROC}_sycl_v100-cuda-11.6"
export MG_EXE="$MG_EXE_DIR/check.exe"
export MG5AMC_CARD_PATH=$MG_PROC_DIR/Cards
```
You may need to adjust these values based upon the target system you are compiling for.
These values worked with CUDA and an NVidia T4 GPU.

Next, compile madgraph4gpu:
```
cd $MG_SP_DIR
make
mv ../../lib/build.d_inl0/ $MG_LIBS_DIR
mv build.d_inl0/ $MG_EXE_DIR
```
If successful test it:
```
cd $WORKSPACE
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MG_LIBS $MG_EXE -j --json_file $WORKSPACE/test_v100_sycl-11.6_${MG_PROC}_${MG_SUBPROC}.json --param_card $MG5AMC_CARD_PATH/param_card.dat --device_id $DEVICE_ID 1024 128 10
```
Check the output in `$WORKSPACE/test_v100_sycl-11.6_${MG_PROC}_${MG_SUBPROC}.json`.
