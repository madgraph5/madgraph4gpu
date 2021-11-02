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
    ###gg_ttg)
    ###  cmd="generate g g > t t~ g"
    ###  ;;
    gg_ttgg)
      cmd="generate g g > t t~ g g"
      ;;
    *)
      echo -e "\nWARNING! Skipping unknown process '$proc'"
      return
      ;;
  esac
  echo -e "\n+++ Generate code for '$proc'\n"
  ###exit 0 # FOR DEBUGGING
  # Generate code for the specific process
  pushd $MG5AMC_HOME >& /dev/null
  outproc=CODEGEN_${proc}
  \rm -rf ${outproc}*
  echo "set stdout_level DEBUG" >> ${outproc}.mg # does not help (log is essentially identical) but add it anyway
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
  # Replace the existing generated code in the output source code directory by the newly generated code and create a .BKP
  rm -rf ${OUTDIR}/${proc}.auto.BKP
  mv ${OUTDIR}/${proc}.auto ${OUTDIR}/${proc}.auto.BKP
  cp -dpr ${MG5AMC_HOME}/${outproc} ${OUTDIR}/${proc}.auto
  echo -e "\nOutput source code has been copied to ${OUTDIR}/${proc}.auto"
  # Compare the existing generated code to the newly generated code for the specific process
  pushd ${OUTDIR} >& /dev/null
  echo -e "\n+++ Compare old and new code generation log for $proc\n"
  ###diff -c ${proc}.auto.BKP/${outproc}_log.txt ${proc}.auto # context diff
  diff ${proc}.auto.BKP/${outproc}_log.txt ${proc}.auto # normal diff
  echo -e "\n+++ Compare old and new generated code for $proc\n"
  if $SCRDIR/diffCode.sh ${BRIEF} -r -c ${proc}.auto.BKP ${proc}.auto; then echo "Old and new generated codes are identical"; else echo -e "\nWARNING! Old and new generated codes differ"; fi
  popd >& /dev/null
  # Compare the existing manually developed code to the newly generated code for the specific process
  pushd ${OUTDIR} >& /dev/null
  echo -e "\n+++ Compare manually developed code to newly generated code for $proc\n"
  if $SCRDIR/diffCode.sh ${BRIEF} -r -c ${proc} ${proc}.auto; then echo "Manual and generated codes are identical"; else echo -e "\nWARNING! Manual and generated codes differ"; fi
  # Print a summary of the available code
  echo -e "Manually developed code is\n  ${OUTDIR}/${proc}"
  echo -e "Old generated code moved to\n  ${OUTDIR}/${proc}.auto.BKP"
  echo -e "New generated code moved to\n  ${OUTDIR}/${proc}.auto"
}

#--------------------------------------------------------------------------------------

function usage()
{
  echo "Usage: $0 [--nobrief] <proc>" # New: only one process
  exit 1
}

#--------------------------------------------------------------------------------------

function cleanup_MG5AMC_HOME()
{
  # Remove MG5aMC fragments from previous runs
  rm -f ${MG5AMC_HOME}/py.py
  rm -f ${MG5AMC_HOME}/Template/LO/Source/make_opts
  rm -f ${MG5AMC_HOME}/input/mg5_configuration.txt
  rm -f ${MG5AMC_HOME}/models/sm/py3_model.pkl
  # Remove and recreate MG5AMC_HOME/PLUGIN
  rm -rf ${MG5AMC_HOME}/PLUGIN
  mkdir ${MG5AMC_HOME}/PLUGIN
  touch ${MG5AMC_HOME}/PLUGIN/__init__.py
}

#--------------------------------------------------------------------------------------

# Default: brief diffs (use --nobrief to use full diffs)
BRIEF=--brief

# Process command line arguments (https://unix.stackexchange.com/a/258514)
for arg in "$@"; do
  shift
  if [ "$arg" == "-h" ] || [ "$arg" == "--help" ]; then
    usage; continue; # continue is unnecessary as usage will exit anyway...
  elif [ "$arg" == "--nobrief" ]; then
    BRIEF=; continue
  else
    # Keep the possibility to collect more then one process
    # However, require a single process to be chosen (allow full cleanup before/after code generation)
    set -- "$@" "$arg"
  fi
done
###procs=$@
if [ "$1" == "" ] || [ "$2" != "" ]; then usage; fi # New: only one process
proc=$1

echo BRIEF=${BRIEF}
###echo procs=${procs}
echo proc=${proc}

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
if [ "$MG5AMC_HOME" == "" ]; then
  echo "ERROR! MG5AMC_HOME is not defined"
  echo "To download MG5AMC please run 'bzr branch lp:~maddevelopers/mg5amcnlo/2.7.0_gpu'"
  exit 1
fi
echo -e "\nUsing MG5AMC_HOME=$MG5AMC_HOME on $(hostname)\n"
if [ ! -d $MG5AMC_HOME ]; then echo "ERROR! Directory $MG5AMC_HOME does not exist"; exit 1; fi

# Print MG5amc bazaar info if any
# Revert to the appropriate bazaar revision number
if bzr --version >& /dev/null; then
  echo -e "Using $(bzr --version | head -1)"
  echo -e "Retrieving bzr information about MG5AMC_HOME"
  if bzr info ${MG5AMC_HOME} 2> /dev/null | grep parent; then
    revno_patches=$(cat $SCRDIR/MG5aMC_patches/2.7.0_gpu/revision.BZR)
    echo -e "MG5AMC patches in this plugin refer to bzr revno '${revno_patches}'"
    echo -e "Revert MG5AMC_HOME to bzr revno '${revno_patches}'"
    bzr revert ${MG5AMC_HOME} -r ${revno_patches}
    revno_mg5amc=$(bzr revno ${MG5AMC_HOME})
    echo -e "Current bzr revno of MG5AMC_HOME is '${revno_mg5amc}'"
    if [ "${revno_patches}" != "${revno_mg5amc}" ]; then echo -e "\nERROR! bzr revno mismatch!"; exit 1; fi
  else
    echo -e "WARNING! MG5AMC_HOME is not a bzr branch\n"
  fi
else
  echo -e "WARNING! bzr is not installed: cannot retrieve bzr properties of MG5aMC_HOME\n"
fi

# Copy MG5AMC patches if any
patches=$(cd $SCRDIR/MG5aMC_patches/2.7.0_gpu; find . -type f -name '*.py')
echo -e "Copy MG5aMC_patches/2.7.0_gpu patches..."
for patch in $patches; do
  patch=${patch#./}
  echo cp -dpr $SCRDIR/MG5aMC_patches/2.7.0_gpu/$patch $MG5AMC_HOME/$patch
  cp -dpr $SCRDIR/MG5aMC_patches/2.7.0_gpu/$patch $MG5AMC_HOME/$patch
done
echo -e "Copy MG5aMC_patches/2.7.0_gpu patches... done\n"

# Clean up before code generation
cleanup_MG5AMC_HOME

# Print MG5amc bazaar info if any
if bzr --version >& /dev/null; then
  if bzr info ${MG5AMC_HOME} 2> /dev/null | grep parent; then
    echo -e "\n***************** Differences to the current bzr revno [START]"
    if bzr diff ${MG5AMC_HOME}; then echo -e "[No differences]"; fi
    echo -e "***************** Differences to the current bzr revno [END]\n"
  fi
fi

# Copy the new plugin to MG5AMC_HOME
cp -dpr ${SCRDIR}/PLUGIN/${OUTBCK^^}_SA_OUTPUT ${MG5AMC_HOME}/PLUGIN/
ls -l ${MG5AMC_HOME}/PLUGIN

# Generate the chosen process (this will always replace the existing code directory and create a .BKP)
codeGenAndDiff $proc

# Clean up after code generation
cleanup_MG5AMC_HOME
