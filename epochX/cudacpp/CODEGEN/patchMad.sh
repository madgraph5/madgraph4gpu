#!/bin/bash

status=0

scrdir=$(cd $(dirname $0); pwd)

if [ "${1%.madonly}" == "$1" ] && [ "${1%.mad}" == "$1" ]; then
  echo "Usage: $0 <process.[madonly|mad]>"
  exit 1 
elif [ "$2" != "" ]; then
  echo "Usage: $0 <process.[madonly|mad]>"
  exit 1 
fi
dir=$1
###echo "Current dir: $pwd"
###echo "Input dir to patch: $dir"

if [ ! -e ${dir} ]; then echo "ERROR! Directory $dir does not exist"; exit 1; fi

# These two steps are part of "cd Source; make" but they actually are code-generating steps
${dir}/bin/madevent treatcards run
${dir}/bin/madevent treatcards param

# Cleanup
\rm -f ${dir}/crossx.html
\rm -f ${dir}/index.html
\rm -f ${dir}/madevent.tar.gz
\rm -rf ${dir}/bin/internal/__pycache__
\rm -rf ${dir}/bin/internal/ufomodel/__pycache__
echo -e "index.html\n.libs\n.pluginlibs" > ${dir}/.gitignore
touch ${dir}/Events/.keepme
\rm -rf ${dir}/HTML

# Inject C++ counters into the Fortran code
\cp -dpr ${scrdir}/PLUGIN/CUDACPP_SA_OUTPUT/madgraph/iolibs/template_files/.clang-format ${dir}
cd ${dir}/SubProcesses
\cp -dpr ${scrdir}/MG5aMC_patches/fbridge_common.inc .
if ! patch -i ${scrdir}/MG5aMC_patches/patch.makefile; then status=1; fi  
cd -
for p1dir in ${dir}/SubProcesses/P1_*; do
  cd $p1dir
  ln -sf ../fbridge_common.inc .
  \cp -dpr ${scrdir}/MG5aMC_patches/counters.cpp .
  if [ "${dir%.mad}" == "$1" ]; then
    \cp -dpr ${scrdir}/PLUGIN/CUDACPP_SA_OUTPUT/madgraph/iolibs/template_files/gpu/timer.h . # already present through cudacpp in *.mad
  fi
  if ! patch -i ${scrdir}/MG5aMC_patches/patch.driver.f; then status=1; fi
  if ! patch -i ${scrdir}/MG5aMC_patches/patch.matrix1.f; then status=1; fi
  if ! patch -i ${scrdir}/MG5aMC_patches/patch.auto_dsig1.f; then status=1; fi
  cd -
done

exit $status
