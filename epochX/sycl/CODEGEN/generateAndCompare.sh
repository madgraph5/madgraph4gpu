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
    gg_ttggg)
      cmd="generate g g > t t~ g g g"
      ;;
    pp_tt)
      cmd="generate p p > t t~"
      ;;
    uu_tt)
      cmd="generate u u~ > t t~"
      ;;
    uu_dd)
      cmd="generate u u~ > d d~"
      ;;
    bb_tt)
      cmd="generate b b~ > t t~"
      ;;
    heft_gg_h)
      cmd="set auto_convert_model T; import model heft; generate g g > h"
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
  outproc=CODEGEN_${OUTBCK}_${proc}
  if [ "${OUTBCK}" == "gridpack" ] && [ "${UNTARONLY}" == "1" ]; then
    echo -e "WARNING! Skip generation of gridpack.tar.gz (--nountaronly was not specified)\n"
  else
    \rm -rf ${outproc} ${outproc}.* ${outproc}_*
    echo "set stdout_level DEBUG" >> ${outproc}.mg # does not help (log is essentially identical) but add it anyway
    echo "${cmd}" >> ${outproc}.mg
    if [ "${OUTBCK}" == "gridpack" ]; then
      if [ "${HELREC}" == "0" ]; then
        echo "output ${outproc} --hel_recycling=False" >> ${outproc}.mg
      else
        echo "output ${outproc}" >> ${outproc}.mg
      fi
      ###echo "!cp -dpr ${outproc} ${outproc}_prelaunch" >> ${outproc}.mg
      echo "launch" >> ${outproc}.mg
      echo "set gridpack True" >> ${outproc}.mg
      echo "set ebeam1 750" >> ${outproc}.mg
      echo "set ebeam2 750" >> ${outproc}.mg
    else
      echo "output standalone_${OUTBCK} ${outproc}" >> ${outproc}.mg
    fi
    cat ${outproc}.mg
    ###{ strace -f -o ${outproc}_strace.txt python3 ./bin/mg5_aMC ${outproc}.mg ; } >& ${outproc}_log.txt
    { time python3 ./bin/mg5_aMC ${outproc}.mg ; } >& ${outproc}_log.txt
    cat ${outproc}_log.txt | egrep -v '(Crash Annotation)' > ${outproc}_log.txt.new # remove firefox 'glxtest: libEGL initialize failed' errors
    \mv ${outproc}_log.txt.new ${outproc}_log.txt
  fi
  if [ -d ${outproc} ] && ! grep -q "Please report this bug" ${outproc}_log.txt; then
    ###cat ${outproc}_log.txt; exit 0 # FOR DEBUGGING
    cat ${MG5AMC_HOME}/${outproc}_log.txt | egrep 'INFO: (Try|Creat|Organiz|Process)'
  else
    echo "*** ERROR! Code generation failed"
    cat ${MG5AMC_HOME}/${outproc}_log.txt
    echo "*** ERROR! Code generation failed"
    exit 1
  fi
  popd >& /dev/null
  # Choose which directory must be copied (for gridpack generation: untar and modify the gridpack)
  if [ "${OUTBCK}" == "gridpack" ]; then
    outprocauto=${MG5AMC_HOME}/${outproc}/run_01_gridpack
    if ! $SCRDIR/untarGridpack.sh ${outprocauto}.tar.gz; then echo "ERROR! untarGridpack.sh failed"; exit 1; fi
  else
    outprocauto=${MG5AMC_HOME}/${outproc}
  fi
  cp -dpr ${MG5AMC_HOME}/${outproc}_log.txt ${outprocauto}/
  # Replace the existing generated code in the output source code directory by the newly generated code and create a .BKP
  rm -rf ${OUTDIR}/${proc}.auto.BKP
  if [ -d ${OUTDIR}/${proc}.auto ]; then mv ${OUTDIR}/${proc}.auto ${OUTDIR}/${proc}.auto.BKP; fi
  cp -dpr ${outprocauto} ${OUTDIR}/${proc}.auto
  echo -e "\nOutput source code has been copied to ${OUTDIR}/${proc}.auto"
  # Compare the existing generated code to the newly generated code for the specific process
  pushd ${OUTDIR} >& /dev/null
  echo -e "\n+++ Compare old and new code generation log for $proc\n"
  ###if diff -c ${proc}.auto.BKP/${outproc}_log.txt ${proc}.auto; then echo "Old and new code generation logs are identical"; fi # context diff
  if diff ${proc}.auto.BKP/${outproc}_log.txt ${proc}.auto; then echo "Old and new code generation logs are identical"; fi # context diff
  echo -e "\n+++ Compare old and new generated code for $proc\n"
  if $SCRDIR/diffCode.sh ${BRIEF} -r -c ${proc}.auto.BKP ${proc}.auto; then echo "Old and new generated codes are identical"; else echo -e "\nWARNING! Old and new generated codes differ"; fi
  popd >& /dev/null
  # Compare the existing manually developed code to the newly generated code for the specific process
  pushd ${OUTDIR} >& /dev/null
  echo -e "\n+++ Compare manually developed code to newly generated code for $proc\n"
  if $SCRDIR/diffCode.sh ${BRIEF} -r -c ${proc} ${proc}.auto; then echo "Manual and generated codes are identical"; else echo -e "\nWARNING! Manual and generated codes differ"; fi
  # Print a summary of the available code
  echo
  echo -e "Manually developed code is\n  ${OUTDIR}/${proc}"
  echo -e "Old generated code moved to\n  ${OUTDIR}/${proc}.auto.BKP"
  echo -e "New generated code moved to\n  ${OUTDIR}/${proc}.auto"
}

#--------------------------------------------------------------------------------------

function usage()
{
  if [ "${OUTBCK}" == "gridpack" ]; then
    echo "Usage: $0 [--nobrief] [--nountaronly] [--nohelrec] <proc>" # New: only one process
  else
    echo "Usage: $0 [--nobrief] <proc>" # New: only one process
  fi
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

# Script directory
SCRDIR=$(cd $(dirname $0); pwd)

# Output source code directory for the chosen backend
OUTDIR=$(dirname $SCRDIR) # e.g. epochX/sycl if $SCRDIR=epochX/sycl/CODEGEN

# Output backend
OUTBCK=$(basename $OUTDIR) # e.g. sycl if $OUTDIR=epochX/sycl

# Default: brief diffs (use --nobrief to use full diffs)
BRIEF=--brief

# Default for gridpacks: untar gridpack.tar.gz but do not regenerate it (use --nountaronly to regenerate it)
UNTARONLY=1

# Default for gridpacks: use helicity recycling (use --nohelrec to disable it)
# (export the value to the untarGridpack.sh script)
export HELREC=1

# Process command line arguments (https://unix.stackexchange.com/a/258514)
for arg in "$@"; do
  shift
  if [ "$arg" == "-h" ] || [ "$arg" == "--help" ]; then
    usage; continue; # continue is unnecessary as usage will exit anyway...
  elif [ "$arg" == "--nobrief" ]; then
    BRIEF=; continue
  elif [ "$arg" == "--nountaronly" ] && [ "${OUTBCK}" == "gridpack" ]; then
    UNTARONLY=0; continue
  elif [ "$arg" == "--nohelrec" ] && [ "${OUTBCK}" == "gridpack" ]; then
    export HELREC=0; continue
  else
    # Keep the possibility to collect more then one process
    # However, require a single process to be chosen (allow full cleanup before/after code generation)
    set -- "$@" "$arg"
  fi
done
###procs=$@
if [ "$1" == "" ] || [ "$2" != "" ]; then usage; fi # New: only one process
proc=$1

echo "SCRDIR=${SCRDIR}"
echo "OUTDIR=${OUTDIR}"
echo "OUTBCK=${OUTBCK} (uppercase=${OUTBCK^^})"

echo "BRIEF=${BRIEF}"
###echo "procs=${procs}"
echo "proc=${proc}"

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
# (NB! 'bzr revert' does not change the output of 'bzr revno': it is NOT like 'git reset --hard'!)
# (See the comments in https://stackoverflow.com/a/37488587)
if bzr --version >& /dev/null; then
  echo -e "Using $(bzr --version | head -1)"
  echo -e "Retrieving bzr information about MG5AMC_HOME"
  if bzr info ${MG5AMC_HOME} > /dev/null; then
    revno_patches=$(cat $SCRDIR/MG5aMC_patches/2.7.0_gpu/revision.BZR)
    echo -e "MG5AMC patches in this plugin refer to bzr revno '${revno_patches}'"
    echo -e "Revert MG5AMC_HOME to bzr revno '${revno_patches}'"
    bzr revert ${MG5AMC_HOME} -r ${revno_patches}
    revno_mg5amc=$(bzr revno ${MG5AMC_HOME} -r ${revno_patches})
    echo -e "Current 'bzr revno -r ${revno_patches}' of MG5AMC_HOME is '${revno_mg5amc}'"
    if [ "${revno_patches}" != "${revno_mg5amc}" ]; then echo -e "\nERROR! bzr revno mismatch!"; exit 1; fi
  else
    ###echo -e "WARNING! MG5AMC_HOME is not a bzr branch\n"
    echo -e "ERROR! MG5AMC_HOME is not a bzr branch\n"; exit 1
  fi
else
  ###echo -e "WARNING! bzr is not installed: cannot retrieve bzr properties of MG5aMC_HOME\n"
  echo -e "ERROR! bzr is not installed: cannot retrieve bzr properties of MG5aMC_HOME\n"; exit 1
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
    echo -e "***************** Differences to the current bzr revno [END]"
  fi
fi

# Copy the new plugin to MG5AMC_HOME (unless this is the gridpack directory)
if [ "${OUTBCK}" != "gridpack" ]; then
  cp -dpr ${SCRDIR}/PLUGIN/${OUTBCK^^}_SA_OUTPUT ${MG5AMC_HOME}/PLUGIN/
  ls -l ${MG5AMC_HOME}/PLUGIN
fi

# For gridpacks, use separate output directories for MG 28x and MG 29x
if [ "${OUTBCK}" == "gridpack" ]; then
  if [ ${revno_patches} -le 365 ]; then
    OUTDIR=${OUTDIR}/28x
  else
    if [ "${HELREC}" == "0" ]; then
      OUTDIR=${OUTDIR}/29x_nohelrec
    else
      OUTDIR=${OUTDIR}/29x
    fi
  fi
  echo "OUTDIR=${OUTDIR} (redefined)"
fi

# Generate the chosen process (this will always replace the existing code directory and create a .BKP)
codeGenAndDiff $proc

# Clean up after code generation
cleanup_MG5AMC_HOME
