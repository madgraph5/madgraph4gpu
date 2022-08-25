if [ "$BASH_SOURCE" = "$0" ]; then echo "ERROR! This script ($0) was not sourced"; exit 1; fi
if [ "$BASH_SOURCE" = "" ]; then echo "ERROR! This script was not sourced from bash"; return 1; fi
scrdir=$(cd $(dirname ${BASH_SOURCE}); pwd)
###. /cvmfs/sft.cern.ch/lcg/releases/gcc/12.1.0-57c96/x86_64-centos7/setup.sh 
. /cvmfs/sft.cern.ch/lcg/releases/gcc/11.2.0-ad950/x86_64-centos7/setup.sh 
export ALLOW_UNSUPPORTED_COMPILER_IN_CUDA=1
export CC=$scrdir/mg-lcg-clang++-14.0.6
export CXX=$scrdir/mg-lcg-clang++-14.0.6

