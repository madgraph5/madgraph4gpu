To build and run the code targeting the CPU, the following commands can be used.

  source /cvmfs/projects.cern.ch/intelsw/oneAPI/linux/x86_64/2021/setvars.sh 
  cd madgraph4gpu/epoch2/sycl/gg_ttgg/
  sed -i s/gpu_selector/cpu_selector/ SubProcesses/P1_Sigma_sm_gg_ttxgg/check_sa.cc 
  cmake -B build
  cmake --build  build
  cd SubProcesses/P1_Sigma_sm_gg_ttxgg/
  ./check_sa.exe 1 4 1


Note that by default the build will result in a JIT kernel. To build the kernel AOT, the following command can be used.
  clang++ -v -pthread -fsycl -fsycl-targets=spir64_x86_64-unknown-unknown-sycldevice -Xs "-march=avx2" -fsycl-unnamed-lambda -I src -I SubProcesses/P1_Sigma_sm_gg_ttxgg/ SubProcesses/P1_Sigma_sm_gg_ttxgg/CPPProcess.cc src/rambo.cc src/read_slha.cc src/Parameters_sm.cc SubProcesses/P1_Sigma_sm_gg_ttxgg/check_sa.cc -o SubProcesses/P1_Sigma_sm_gg_ttxgg/check_sa.exe

To target Nvidia platforms, a custom build of the LLVM is required. This has been added to /opt on itscrd02. The following commands can be used to build and run the code.

  export PATH=/opt/sycl_workspace/llvm/build/bin:$PATH
  export LD_LIBRARY_PATH=/opt/sycl_workspace/llvm/build/lib:$LD_LIBRARY_PATH

  cd madgraph4gpu/epoch2/sycl/gg_ttgg/

clang++ -v -ffast-math -pthread -fsycl -fsycl-targets=nvptx64-nvidia-cuda-sycldevice --cuda-include-ptx=sm_70 -fsycl-unnamed-lambda --cuda-path=/usr/local/cuda-10.2/ -I src -I SubProcesses/P1_Sigma_sm_gg_ttxgg/ -I ../../../tools/ src/read_slha.cc src/Parameters_sm.cc SubProcesses/P1_Sigma_sm_gg_ttxgg/check_sa.cc -o SubProcesses/P1_Sigma_sm_gg_ttxgg/gcheck_sa.exe

 cd SubProcesses/P1_Sigma_sm_gg_ttxgg
 chmod +x gcheck_sa.exe
 ./gcheck_sa.exe 1 4 1


