#!/bin/bash

status=0

scrdir=$(cd $(dirname $0); pwd)

function usage()
{
  echo "Usage: $0 <process.[madonly|mad]> <vecsize> [--nopatch]"
  exit 1 
}

nopatch=0 # apply patch commands (default)
if [ "${1%.madonly}" == "$1" ] && [ "${1%.mad}" == "$1" ]; then
  usage
elif [ "$2" == "" ]; then
  usage
elif [ "$3" == "--nopatch" ]; then
  if [ "$4" != "" ]; then usage; fi
  nopatch=1 # skip patch commands (produce a new reference)
elif [ "$3" != "" ]; then
  usage
fi
dir=$1
vecsize=$2
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
echo -e "index.html\n.libs\n.cudacpplibs" > ${dir}/.gitignore
touch ${dir}/Events/.keepme
\rm -rf ${dir}/HTML

# Patch the default Fortran code to provide the integration with the cudacpp plugin
# (1) Process-independent patches
\cp -dpr ${scrdir}/PLUGIN/CUDACPP_SA_OUTPUT/madgraph/iolibs/template_files/.clang-format ${dir} # new file
\cp -dpr ${scrdir}/MG5aMC_patches/vector.inc ${dir}/Source # replace default
\cp -dpr ${scrdir}/MG5aMC_patches/fbridge_common.inc ${dir}/SubProcesses # new file
for file in ${dir}/Source/MODEL/rw_para.f; do cat ${file} | sed "s|include 'coupl.inc'|include 'vector.inc'\n      include 'coupl.inc'|" > ${file}.new; \mv ${file}.new ${file}; done
for file in ${dir}/Source/PDF/ElectroweakFlux.inc; do cat ${file} | sed "s|include '../MODEL/coupl.inc'|include 'vector.inc'\n        include 'coupl.inc'|" > ${file}.new; \mv ${file}.new ${file}; done
cd ${dir}/SubProcesses
cd - > /dev/null
if [ "${nopatch}" == "0" ]; then
  cd ${dir}
  if ! patch -p4 -i ${scrdir}/MG5aMC_patches/patch.common; then status=1; fi  
  cd - > /dev/null
fi
for p1dir in ${dir}/SubProcesses/P1_*; do
  cd $p1dir
  echo -e "madevent\n*madevent_cudacpp" > .gitignore # new file
  ln -sf ../fbridge_common.inc . # new file
  \cp -dpr ${scrdir}/MG5aMC_patches/counters.cpp . # new file
  if [ "${dir%.mad}" == "$1" ]; then
    \cp -dpr ${scrdir}/PLUGIN/CUDACPP_SA_OUTPUT/madgraph/iolibs/template_files/gpu/timer.h . # new file, already present via cudacpp in *.mad
  fi
  if [ "${nopatch}" == "0" ]; then
    if ! patch -p6 -i ${scrdir}/MG5aMC_patches/patch.P1; then status=1; fi  
  fi
  cd - > /dev/null
done

# Patch the default Fortran code to provide the integration with the cudacpp plugin
# (2) Process-dependent patches
cd ${dir}/Source/MODEL > /dev/null
echo "c NB vector.inc (defining nb_page_max) must be included before coupl.inc (#458)" > coupl.inc.new
cat coupl.inc | sed "s/($vecsize)/(NB_PAGE_MAX)/g" >> coupl.inc.new
\mv coupl.inc.new coupl.inc
gcs=$(cat coupl_write.inc | awk '{if($1=="WRITE(*,2)") print $NF}') # different printouts for scalar/vector couplings #456
for gc in $gcs; do
  if grep -q "$gc(NB_PAGE_MAX)" coupl.inc; then
    ###echo "DEBUG: Coupling $gc is a vector"
    cat coupl_write.inc | awk -vgc=$gc '{if($1=="WRITE(*,2)" && $NF==gc) print $0"(1)"; else print $0}' > coupl_write.inc.new
    \mv coupl_write.inc.new coupl_write.inc
  ###else
  ###  echo "DEBUG: Coupling $gc is a scalar"
  fi
done
for file in couplings.f couplings1.f couplings2.f; do cat ${file} | sed "s|INCLUDE 'coupl.inc'|include 'vector.inc'\n      include 'coupl.inc'|" > ${file}.new; \mv ${file}.new ${file}; done
cd - > /dev/null
cd ${dir}/SubProcesses > /dev/null
echo "c NB vector.inc (defining nb_page_max) must be included before clusters.inc (#458)" > cluster.inc.new
cat cluster.inc | grep -v "      include 'vector.inc'" | sed "s/nb_page/nb_page_max/g" >> cluster.inc.new
\mv cluster.inc.new cluster.inc
cd - > /dev/null
for p1dir in ${dir}/SubProcesses/P1_*; do
  cd $p1dir > /dev/null
  cat auto_dsig1.f \
    | sed "s|NB_PAGE)|NB_PAGE_MAX)|" \
    | sed "s|1,NB_PAGE|1,NB_PAGE_LOOP|" \
    | sed "s|1, NB_PAGE|1, NB_PAGE_LOOP|" \
    | sed "s|/NB_PAGE|/NB_PAGE_MAX|" \
    > auto_dsig1.f.new
  \mv auto_dsig1.f.new auto_dsig1.f
  if [ "${nopatch}" == "0" ]; then
    if ! patch -p6 -i ${scrdir}/MG5aMC_patches/patch.auto_dsig1.f; then status=1; fi  
  fi
  \rm -f *.orig
  ncolor=$(cat matrix1.f | grep PARAMETER | grep NCOLOR= | sed 's/.*NCOLOR=//' | sed 's/)//')
  cat auto_dsig1.f \
    | sed "s|JAMP2_MULTI(0,IVEC) = 2 ! workaround|JAMP2_MULTI(0,IVEC) = ${ncolor} ! workaround|" \
    > auto_dsig1.f.new
  \mv auto_dsig1.f.new auto_dsig1.f
  cd - > /dev/null
done

exit $status
