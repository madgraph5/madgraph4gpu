# Documentation

We are currently using [GitHub Actions](https://docs.github.com/en/actions) in conjunction with onsite self-hosted [GitHub Runners](https://docs.github.com/en/actions/hosting-your-own-runners/about-self-hosted-runners) to automate compiling/testing and performance profiling tasks in SYCL and CUDA on A100 and V100s GPUs currently.

## Grafana link: [madgraph4gpu-db.web.cern.ch](https://madgraph4gpu-db.web.cern.ch/)

## Performance Profiling

### Profiling baseline currently used

**GCC - 11.3.0**

**CUDA - 12.0.1**

**Clang - 16**

### GitHub Actions Runner

A [GitHub Runner](https://docs.github.com/en/actions/using-github-hosted-runners/about-github-hosted-runners) is a tool that allows users to automate their workflow by running [actions](https://docs.github.com/en/actions) or tasks in response to specific events on GitHub. This can include tasks such as running tests, building and deploying code, or publishing artifacts. They can be easily configured and managed through the GitHub website, and can help users streamline their development process and ensure that their code is always up-to-date and ready for deployment. In our case we use them to automate CI and nightly performance profiling.

### performanceProfiler.py

This is the main entrypoint for the profiler. It executes the two bash build scripts for SYCL ```buildSYCLProcess.sh``` and CUDA ```buildCUDAProcess.sh``` with the correct ThreadsPerBlock, BlocksPerThread and iteration count.

#### Usage:

Go to the `tools/profiling` directory and run:

```
python3 performanceProfiler.py -l <abstration_layer> -b <branch>
```

The following options are available for this script:

`-l`: This option specifies the abstraction layer to use for profiling. The supported values are "SYCL" and "CUDA". The default value is "SYCL".

`-b`: This option specifies the branch of the madgraph4gpu repository that will be used. The default value is "master".

Example:

Copy code
python script.py
To run the script with a different abstraction layer and branch, you can use the following command:

```
python script.py -l CUDA -b my_branch
```

### buildSYCLProcess.sh

This bash script compiles and executes standalone physics processes using the MadGraph5_aMC@NLO GPU development framework with oneAPI/SYCL.

#### Usage

Go to the `tools/profiling` directory and run:

```
./buildSYCLProcess.sh -n <physics_process> -b <blocks_per_grid> -t <threads_per_block> -i <iterations> [-r <branch>] [-d <device_id/info>]
```

#### Arguments:

* `-n`: Name of the physics process being built and run (e.g., gg_ttgg).

* `-b`: Number of blocks per grid.

* `-t`: Number of threads per block.

* `-i`: Number of iterations.

* `-r`: (Optional) Branch name. Default: not displayed in the report folder prefix.

* `-d`: (Optional) Flag for setting the device ID. Default: "--device_id 2" for oneAPI toolkit runs on GPUs, otherwise "--device_id 0" for LLVM DPCPP compiler. You can also use `-d info` to get the specific device IDs for that host.

#### Example:

```
./buildSYCLProcess.sh -n gg_ttgg -b 1024 -t 128 -i 10 -r master -d 2
```

**Note**:

To also compile to CPUs you need to enable more backends in the DPCPP toolchain (Currently when you follow how to use the LLVM DPCPP compiler for CUDA it does not install the necessary dependencies to see other devices as well on the host). You can read more on how to enable more backends [here](https://intel.github.io/llvm-docs/GetStartedGuide.html#build-dpc-toolchain).

### buildCUDAProcess.sh

This script compiles and executes physics processes using the MadGraph5_aMC@NLO GPU development framework with CUDA.

#### Usage

Go to the `tools/profiling` directory and run:

```
./buildCUDAProcess.sh -n <process_name> -b <blocks_per_grid> -t <threads_per_block> -i <iterations> -r <branch> -m <make_args>
```

#### Arguments:

* `-n`: Name of the physics process being built and run.

* `-b`: Number of blocks per grid.

* `-t`: Number of threads per block.

* `-i`: Number of iterations.

* `-r`: Branch name.

* `-m`: Makefile arguments.

#### Example:

```
./buildCUDAProcess.sh -n gg_ttgg -b 1024 -t 128 -i 10 -r master -m avx2
```

#### Notes

This script assumes that it is run from the profiling directory in the repository.
Make sure to set the correct CUDA path according to your system.
You may need to modify the script to set the correct GPU architecture or compiler options depending on your system.

### sendData.py

#### Usage:

Go to the `tools/profiling` directory and run:

```
python3 sendData.py -r <report_folder_to_upload> -branch <master>
```

The following arguments are available for this script:

* `-r` or `--reportPath```: This argument specifies the path for the reports that will be sent to the database.

* `-f` or `--fields`: This argument specifies the fields in the JSON data that will be sent to the database. The default value is `['EvtsPerSec[MatrixElems] (3)', 'EvtsPerSec[MECalcOnly] (3)']`.

* `-b` or `--branch`: This argument specifies the branch that the profiler data is in. The default value is `master`.

* `-p` or `--profiler`: This argument enables CI profiling defaults. The default value is `False`.

For example, to run the script with the default arguments, you can use the following command:

python3 sendData.py
To run the script with a custom report path and branch, you can use the following command:

python3 sendData.py -r /path/to/reports -b my_branch
Note that some options may not be relevant or may not work as expected in certain situations. For example, the -p option will only work when CI profiling defaults are enabled.

## Known issues:

### Bug in GCC 11.3.0/11.3.1 using the LLVM DPCPP compiler 

There is a [bug](https://bugs.gentoo.org/842405) affecting GCC versions 11.3.0/11.3.1 when compiling the standalone physics processes resulting in two compilation errors `.../fs_path.h:1209:9: error: 'end' is missing exception specification 'noexcept'` and `.../fs_path.h:1217:9: error: 'end' is missing exception specification 'noexcept'`` in the `fs_path.h` file. GCC version 11.2.0 is not affected, and appears to be fixed in later versions (Remains to be tested and cited).

### libmg5amc_common.so: cannot open shared object file: No such file or directory

The libmg5amc_common.so library is not set in the LD_LIBRARY_PATH

### Not linking correctly/Wrong linker version from what you intend to compile with?

If you have problems with wrong linker see which candidate GCC finds with `./sycl_workspace/llvm/build/bin/clang++ -v` and see if it is the correct GCC candidate. If it is not, you can correct this with adding `--gcc-toolchain=/cvmfs/sft.cern.ch/lcg/releases/gcc/11.3.0-ad0f5/x86_64-centos8/lib/gcc/x86_64-pc-linux-gnu/11.3.0` to the `CXXFLAGS`. This will correctly set the GCC candidate to the desired GCC installation. Using `ENABLE_CI_PROFILER=1` automatically adds this in all the standalone physics processes makefiles in SYCL and in CUDA.
