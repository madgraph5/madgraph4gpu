module purge

module load gcc/9.3.0/cmake
module load cuda/11.6.1

export CXX=$(which g++)
export NVCC=$(which nvcc)
echo "Using CXX=${CXX}"
echo "Using NVCC=${NVCC}"

if gfortran-9 --version >& /dev/null; then
  export FC=$(which gfortran-9)
  export LIBRARY_PATH=$(dirname $(${FC} --print-file-name libgfortran.so)):$LIBRARY_PATH
elif gfortran-8 --version >& /dev/null; then
  export FC=$(which gfortran-8)
  export LIBRARY_PATH=$(dirname $(${FC} --print-file-name libgfortran.so)):$LIBRARY_PATH
else
  echo "WARNING! Neither gfortran-9 nor gfortran-8 were found in PATH"
fi
echo "Using FC=${FC}"
echo "Using LIBRARY_PATH=${LIBRARY_PATH}"
