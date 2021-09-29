#!/bin/bash

#--------------------------------------------------------------------------------------

function codeGenAndDiff()
{
  proc=$1
  # Process-dependent hardcoded configuration
  echo -e "\n================================================================"
  case "${proc}" in
    ee_mumu)
      cmd="generate e+ e- > mu+ mu-"
      ;;
    gg_tt)
      cmd="generate g g > t t~"
      ;;
    gg_ttg)
      cmd="generate g g > t t~ g"
      ;;
    gg_ttgg)
      cmd="generate g g > t t~ g g"
      ;;
    *)
      echo -e "\nWARNING! Skipping unknown process '$proc'"
      return
      ;;
  esac
  echo -e "\n+++ Generate code for '$proc'\n"
  # Generate code for the specific process
  pushd $MG5AMC_HOME >& /dev/null
  outproc=CODEGEN_${proc}
  \rm -rf ${outproc}*
  ###echo "set stdout_level DEBUG" >> ${outproc}.mg # does not help (log is essentially identical)
  echo "${cmd}" >> ${outproc}.mg
  echo "output standalone_${OUTBCK} ${outproc}" >> ${outproc}.mg
  cat  ${outproc}.mg
  ###{ strace -f -o ${outproc}_strace.txt python3 ./bin/mg5_aMC ${outproc}.mg ; } >& ${outproc}_log.txt
  { time python3 ./bin/mg5_aMC ${outproc}.mg ; } >& ${outproc}_log.txt
  if [ -d ${outproc} ] && ! grep -q "Please report this bug" ${outproc}_log.txt; then
    ###cat ${outproc}_log.txt; exit 0 # FOR DEBUGGING
    cat ${outproc}_log.txt | egrep 'INFO: (Try|Creat|Organiz|Process)'
    mv ${outproc}_log.txt ${outproc}/
  else
    echo "*** ERROR! Code generation failed"
    cat ${outproc}_log.txt
    echo "*** ERROR! Code generation failed"
    exit 1
  fi
  popd >& /dev/null
  # Move the newly generated code to the output source code directory
  rm -rf ${OUTDIR}/${proc}.auto.BKP ${OUTDIR}/${proc}.auto.NEW
  cp -dpr ${MG5AMC_HOME}/${outproc} ${OUTDIR}/${proc}.auto.NEW
  echo -e "\nOutput source code has been copied to ${OUTDIR}/${proc}.auto.NEW"
  # Compare the newly generated code to the existing generated code for the specific process
  pushd ${OUTDIR} >& /dev/null
  echo -e "\n+++ Compare new and old generated code for $proc\n"
  if diff ${BRIEF} -x '*log.txt' -x '*.o' -x '*.o.*' -x '*.a' -x '*.exe' -x 'lib' -r -c ${proc}.auto.NEW ${proc}.auto; then echo "New and old generated codes are identical"; else echo -e "\nWARNING! New and old generated codes differ"; fi
  echo -e "\n+++ Compare new and old code generation log for $proc\n"
  ###diff -c ${proc}.auto.NEW/${outproc}_log.txt ${proc}.auto # context diff
  diff ${proc}.auto.NEW/${outproc}_log.txt ${proc}.auto # normal diff
  popd >& /dev/null
  # Compare the newly generated code to the existing manually developed code for the specific process
  pushd ${OUTDIR} >& /dev/null
  echo -e "\n+++ Compare newly generated code to manually developed code for $proc\n"
  if diff ${BRIEF} -x '*log.txt' -x '*.o' -x '*.o.*' -x '*.a' -x '*.exe' -x 'lib' -r -c ${proc}.auto.NEW ${proc}; then echo "Generated and manual codes are identical"; else echo -e "\nWARNING! Generated and manual codes differ"; fi
  # Replace the existing generated code by the newly generated code if required
  if [ "${REPLACE}" == "1" ]; then
    echo -e "\n+++ Replace existing generated code for $proc (REPLACE=$REPLACE)\n"
    mv ${OUTDIR}/${proc}.auto ${OUTDIR}/${proc}.auto.BKP
    mv ${OUTDIR}/${proc}.auto.NEW ${OUTDIR}/${proc}.auto
    echo -e "Manually developed code is\n  ${OUTDIR}/${proc}"
    echo -e "Old generated code moved to\n  ${OUTDIR}/${proc}.auto.BKP"
    echo -e "New generated code moved to\n  ${OUTDIR}/${proc}.auto"
  else
    echo -e "\n+++ Keep existing generated code for $proc (REPLACE=$REPLACE)\n"
    echo -e "Manually developed code is\n  ${OUTDIR}/${proc}"
    echo -e "Old generated code is\n  ${OUTDIR}/${proc}.auto"
    echo -e "New generated code is\n  ${OUTDIR}/${proc}.auto.NEW"
  fi
}

#--------------------------------------------------------------------------------------

function usage()
{
  echo "Usage: $0 [--replace|--noreplace] [--brief] [<proc1> [... <procN>]]"
  exit 1
}

#--------------------------------------------------------------------------------------

# Replace code directory and create .BKP? (or alternatively keep code directory in .NEW?)
REPLACE=0

# Brief diffs?
BRIEF=

# Process command line arguments (https://unix.stackexchange.com/a/258514)
for arg in "$@"; do
  shift
  if [ "$arg" == "-h" ] || [ "$arg" == "--help" ]; then
    usage; continue; # continue is unnecessary as usage will exit anyway...
  elif [ "$arg" == "--replace" ]; then
    REPLACE=1; continue;
  elif [ "$arg" == "--noreplace" ]; then
    REPLACE=0; continue;
  elif [ "$arg" == "--brief" ]; then
    BRIEF=--brief; continue
  else
    set -- "$@" "$arg"
  fi
done
procs=$@
echo REPLACE=${REPLACE}
echo BRIEF=${BRIEF}
echo procs=${procs}

# Script directory
SCRDIR=$(cd $(dirname $0); pwd)
echo SCRDIR=${SCRDIR}

# Output source code directory for the chosen backend
OUTDIR=$(dirname $SCRDIR) # e.g. epochX/cudacpp if $SCRDIR=epochX/cudacpp/CODEGEN
echo OUTDIR=${OUTDIR}

# Output backend
OUTBCK=$(basename $OUTDIR) # e.g. cudacpp if $OUTDIR=epochX/cudacpp
echo "OUTBCK=${OUTBCK} (uppercase=${OUTBCK^^})"

# Make sure that python3 is installed
if ! python3 --version >& /dev/null; then echo "ERROR! python3 is not installed"; exit 1; fi

# Make sure that $MG5AMC_HOME exists
if [ "$MG5AMC_HOME" == "" ]; then echo "ERROR! MG5AMC_HOME is not defined"; exit 1; fi
echo -e "\nUsing MG5AMC_HOME=$MG5AMC_HOME on $(hostname)\n"
if [ ! -d $MG5AMC_HOME ]; then echo "ERROR! Directory $MG5AMC_HOME does not exist"; exit 1; fi

# Copy MG5AMC patches if any
patches=$(cd $SCRDIR/MG5aMC_patches/2.7.0_gpu; find . -type f -name '*.py')
echo -e "Copy MG5aMC_patches/2.7.0_gpu patches..."
for patch in $patches; do
  echo cp -dpr $SCRDIR/MG5aMC_patches/2.7.0_gpu/$patch $MG5AMC_HOME/$patch
  cp -dpr $SCRDIR/MG5aMC_patches/2.7.0_gpu/$patch $MG5AMC_HOME/$patch
done
echo -e "Copy MG5aMC_patches/2.7.0_gpu patches... done\n"

# Remove MG5aMC fragments from previous runs
rm -rf ${MG5AMC_HOME}/py.py

# Remove and recreate MG5AMC_HOME/PLUGIN
rm -rf ${MG5AMC_HOME}/PLUGIN
mkdir ${MG5AMC_HOME}/PLUGIN
touch ${MG5AMC_HOME}/PLUGIN/__init__.py
cp -dpr ${SCRDIR}/PLUGIN/${OUTBCK^^}_SA_OUTPUT ${MG5AMC_HOME}/PLUGIN/
ls -lR $MG5AMC_HOME/PLUGIN

# Determine the list of processes to generate
###procs="ee_mumu gg_tt gg_ttg gg_ttgg"
if [ "$procs" == "" ] ; then procs=$(cd $OUTDIR; find . -mindepth 1 -maxdepth 1 -type d -name '*.auto' | sed 's/.auto//'); fi

# Iterate through the list of processes to generate
for proc in $procs; do
  if [ -d $OUTDIR/$proc ]; then proc=$(basename $proc); fi
  codeGenAndDiff $proc
done
