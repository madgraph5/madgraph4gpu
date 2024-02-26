if [ "$BASH_SOURCE" = "$0" ]; then echo "ERROR! This script ($0) was not sourced"; exit 1; fi
if [ "$BASH_SOURCE" = "" ]; then echo "ERROR! This script was not sourced from bash"; return 1; fi
scrdir=$(cd $(dirname ${BASH_SOURCE}); pwd)
redrel=$(cat /etc/redhat-release 2> /dev/null)
if [ "${redrel##*release 7}" != "${redrel}" ]; then
  source /cvmfs/sft.cern.ch/lcg/releases/clang/14.0.6-14bdb/x86_64-centos7/setup.sh
elif [ "${redrel##*release 9}" != "${redrel}" ]; then
  source /cvmfs/sft.cern.ch/lcg/releases/clang/14.0.6-14bdb/x86_64-centos9/setup.sh
else
  echo "ERROR! RedHat release ${redrel} is not supported by ${BASH_SOURCE}"
  return 1
fi
###export ALLOW_UNSUPPORTED_COMPILER_IN_CUDA=1
export CC=$scrdir/mg-lcg-clang-14.0.6
export CXX=$scrdir/mg-lcg-clang++-14.0.6

