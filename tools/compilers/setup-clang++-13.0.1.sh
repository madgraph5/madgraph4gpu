if [ "$BASH_SOURCE" = "$0" ]; then echo "ERROR! This script ($0) was not sourced"; exit 1; fi
if [ "$BASH_SOURCE" = "" ]; then echo "ERROR! This script was not sourced from bash"; return 1; fi
scrdir=$(cd $(dirname ${BASH_SOURCE}); pwd)
source /cvmfs/sft.cern.ch/lcg/releases/binutils/2.37-355ed/x86_64-centos7/setup.sh
source /cvmfs/sft.cern.ch/lcg/releases/gcc/11.2.0-8a51a/x86_64-centos7/setup.sh
###export ALLOW_UNSUPPORTED_COMPILER_IN_CUDA=1
export CC=$scrdir/mg-lcg-clang++-13.0.1
export CXX=$scrdir/mg-lcg-clang++-13.0.1

