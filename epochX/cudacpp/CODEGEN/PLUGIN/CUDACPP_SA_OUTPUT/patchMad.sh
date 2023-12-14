#!/bin/bash
# Copyright (C) 2020-2023 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: A. Valassi (Mar 2022) for the MG5aMC CUDACPP plugin.
# Further modified by: O. Mattelaer, A. Valassi (2022-2023) for the MG5aMC CUDACPP plugin.

set -e # immediate exit on error

status=0

scrdir=$(cd $(dirname $0); pwd)

function usage()
{
  echo "ERROR! Unknown command '$0 $*'"
  echo "Usage: $0 <process_dir> <patch_dir> [--nopatch|--upstream]"
  exit 1 
}

# Patch level
###patchlevel=0 # [--upstream] out of the box codegen from upstream MG5AMC (do not even copy templates)
###patchlevel=1 # [--nopatch] modify upstream MG5AMC but do not apply patch commands (reference to prepare new patches)
patchlevel=2 # [DEFAULT] complete generation of cudacpp .sa/.mad (copy templates and apply patch commands)

if [ "$2" == "" ]; then
  usage $*
elif [ "$3" == "--nopatch" ]; then
  if [ "$4" != "" ]; then usage; fi
  patchlevel=1
elif [ "$3" == "--upstream" ]; then
  if [ "$4" != "" ]; then usage; fi
  patchlevel=0
elif [ "$3" != "" ]; then
  usage $*
fi
dir=$1
dir_patches=$2
###echo "Current dir: $pwd"
###echo "Input dir to patch: $dir"

if [ ! -e ${dir} ]; then echo "ERROR! Directory $dir does not exist"; exit 1; fi

# Exit here for patchlevel 0 (--upstream)
if [ "${patchlevel}" == "0" ]; then exit $status; fi

# Patch the default Fortran code to provide the integration with the cudacpp plugin
# (1) Process-independent patches
touch ${dir}/Events/.keep # this file should already be present (mg5amcnlo copies it from Template/LO/Events/.keep) 
\cp -pr ${scrdir}/MG5aMC_patches/${dir_patches}/fbridge_common.inc ${dir}/SubProcesses # new file
cd ${dir}
echo "DEBUG: standardise ${PWD}/Source/make_opts (fix f2py3 and sort make_opts_variables) before applying patch.common" # AV moved back hare in patchMad.sh from generateAndCompare.sh (see PR #753)
sed -i 's/DEFAULT_F2PY_COMPILER=f2py.*/DEFAULT_F2PY_COMPILER=f2py3/' Source/make_opts 
cat Source/make_opts | sed '/#end/q' | head --lines=-1 | sort > Source/make_opts.new
cat Source/make_opts | sed -n -e '/#end/,$p' >> Source/make_opts.new
\mv Source/make_opts.new Source/make_opts
if [ "${patchlevel}" == "2" ]; then
  echo "DEBUG: cd ${PWD}; patch -p4 -i ${scrdir}/MG5aMC_patches/${dir_patches}/patch.common"
  if ! patch -p4 -i ${scrdir}/MG5aMC_patches/${dir_patches}/patch.common; then status=1; fi  
  \rm -f Source/*.orig
  \rm -f bin/internal/*.orig
fi
cd - > /dev/null
for p1dir in ${dir}/SubProcesses/P*; do
  cd $p1dir
  ln -sf ../fbridge_common.inc . # new file
  cp -pr ${scrdir}/MG5aMC_patches/${dir_patches}/counters.cc . # new file
  cp -pr ${scrdir}/MG5aMC_patches/${dir_patches}/ompnumthreads.cc . # new file
  ###cp -pr ${scrdir}/MG5aMC_patches/${dir_patches}/counters.cc ${dir}/SubProcesses/ # new file (SH)
  ###cp -pr ${scrdir}/MG5aMC_patches/${dir_patches}/ompnumthreads.cc ${dir}/SubProcesses/ # new file (SH)
  ###ln -sf ../counters.cc . # new file (SH)
  ###ln -sf ../ompnumthreads.cc . # new file (SH)
  if [ "${patchlevel}" == "2" ]; then
    echo "DEBUG: cd ${PWD}; patch -p6 -i ${scrdir}/MG5aMC_patches/${dir_patches}/patch.P1"
    if ! patch -p6 -i ${scrdir}/MG5aMC_patches/${dir_patches}/patch.P1; then status=1; fi      
  fi
  \rm -f *.orig
  cd - > /dev/null
done

# Patch the default Fortran code to provide the integration with the cudacpp plugin
# (2) Process-dependent patches
cd ${dir}/Source/MODEL > /dev/null
gcs=$(cat coupl_write.inc | awk '{if($1=="WRITE(*,2)") print $NF}') # different printouts for scalar/vector couplings #456
for gc in $gcs; do
  if grep -q "$gc(VECSIZE_MEMMAX)" coupl.inc; then
    ###echo "DEBUG: Coupling $gc is a vector"
    cat coupl_write.inc | awk -vgc=$gc '{if($1=="WRITE(*,2)" && $NF==gc) print $0"(1)"; else print $0}' > coupl_write.inc.new
    \mv coupl_write.inc.new coupl_write.inc
  ###else
  ###  echo "DEBUG: Coupling $gc is a scalar"
  fi
done
cd - > /dev/null

# Patch the default cudacpp code to fix a bug in coloramps
# ** NEW AUG 2023: DISABLING THE COLORAMPS PATCH FIXES THE LHE COLOR MISMATCH IN GG_TTGG (#655 and #713) **
# (3) Process-dependent patches
#for p1dir in ${dir}/SubProcesses/P*; do
#  cd $p1dir
#  cat coloramps.h | awk -vp=1 '{if (p==1) print $0; if ($1=="__device__") p=0}' > coloramps.h.new
#  cat coloramps.inc | sed 's|)/|)/ {|' | sed 's|/$|}, /|' \
#    | awk -vb= '{if($1~")/"){b=$2}; if($1=="$"){b=b$2}; if($3=="/"){print "    "b}}' \
#    | sed 's/.TRUE./ true/g' | sed 's/.FALSE./ false/g' | sed 's/}/ }/' >> coloramps.h.new
#  truncate -s -2 coloramps.h.new
#  echo "">> coloramps.h.new
#  cat coloramps.h | awk -vp=0 '{if ($1=="};") p=1; if (p==1) print $0}' >> coloramps.h.new
#  \mv coloramps.h.new coloramps.h
#done

exit $status
