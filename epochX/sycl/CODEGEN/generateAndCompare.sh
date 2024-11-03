#!/bin/bash
# Copyright (C) 2020-2023 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: A. Valassi (Sep 2021) for the MG5aMC CUDACPP plugin.
# Further modified by: A. Valassi (2021-2023) for the MG5aMC CUDACPP plugin.
#
# Copyright (C) 2021-2023 Argonne National Laboratory.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Modified by: N. Nichols (2021-2023) for the MG5aMC SYCL plugin.

set -e # fail on error

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
    gg_ttgggg)
      cmd="generate g g > t t~ g g g g"
      ;;
    gg_tt01g)
      cmd="generate g g > t t~; add process g g > t t~ g"
      ;;
    gq_ttq)
      cmd="define q = u c d s u~ c~ d~ s~; generate g q > t t~ q"
      ;;
    gq_ttllq)
      cmd="define q = u c d s u~ c~ d~ s~; generate g q > t t~ l- l+ q"
      ;;
    pp_tt)
      cmd="generate p p > t t~"
      ;;
    pp_tttt)
      cmd="generate p p > t t~ t t~"
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
    smeft_gg_tttt)
      cmd="set auto_convert_model T; import model SMEFTsim_topU3l_MwScheme_UFO -massless_4t; generate g g > t t~ t t~"
      ;;
    susy_gg_tt)
      cmd="import model MSSM_SLHA2; generate g g > t t~"
      ;;
    *)
      echo -e "\nWARNING! Skipping unknown process '$proc'"
      return
      ;;
  esac
  echo -e "\n+++ Generate code for '$proc'\n"
  ###exit 0 # FOR DEBUGGING
  # Vector size for mad/madonly meexporter (VECSIZE_MEMMAX)
  vecsize=16384 # NB THIS IS IGNORED ANYWAY (ALL HARDCODED VALUES ARE REPLACED IN PATCHMAD.SH)...
  # Generate code for the specific process
  pushd $MG5AMC_HOME >& /dev/null
  outproc=CODEGEN_${OUTBCK}_${proc}
  if [ "${SCRBCK}" == "gridpack" ] && [ "${UNTARONLY}" == "1" ]; then
    echo -e "WARNING! Skip generation of gridpack.tar.gz (--nountaronly was not specified)\n"
  else
    \rm -rf ${outproc} ${outproc}.* ${outproc}_*
    if [ "${HELREC}" == "0" ]; then
      helrecopt="--hel_recycling=False"
    else
      helrecopt=
    fi
    echo "set stdout_level DEBUG" >> ${outproc}.mg # does not help (log is essentially identical) but add it anyway
    echo "set zerowidth_tchannel F" >> ${outproc}.mg # workaround for #476: do not use a zero top quark width in fortran (~E-3 effect on physics)
    echo "${cmd}" >> ${outproc}.mg
    if [ "${SCRBCK}" == "gridpack" ]; then # $SCRBCK=$OUTBCK=gridpack
      echo "output ${outproc} ${helrecopt}" >> ${outproc}.mg
      ###echo "!cp -dpr ${outproc} ${outproc}_prelaunch" >> ${outproc}.mg
      echo "launch" >> ${outproc}.mg
      echo "set gridpack True" >> ${outproc}.mg
      echo "set ebeam1 750" >> ${outproc}.mg
      echo "set ebeam2 750" >> ${outproc}.mg
    elif [ "${SCRBCK}" == "alpaka" ]; then # $SCRBCK=$OUTBCK=alpaka
      echo "output standalone_${SCRBCK}_sycl ${outproc}" >> ${outproc}.mg
    elif [ "${OUTBCK}" == "madnovec" ]; then # $SCRBCK=sycl and $OUTBCK=madnovec
      echo "output madevent ${outproc} ${helrecopt}" >> ${outproc}.mg
    elif [ "${OUTBCK}" == "madonly" ]; then # $SCRBCK=sycl and $OUTBCK=madonly
      echo "output madevent ${outproc} ${helrecopt} --vector_size=${vecsize}" >> ${outproc}.mg
    elif [ "${OUTBCK}" == "mad" ]; then # $SCRBCK=sycl and $OUTBCK=mad
      echo "output madevent ${outproc} ${helrecopt} --vector_size=${vecsize} --me_exporter=standalone_sycl" >> ${outproc}.mg
    elif [ "${OUTBCK}" == "madcpp" ]; then # $SCRBCK=sycl and $OUTBCK=madcpp
      echo "output madevent ${outproc} ${helrecopt} --vector_size=32 --me_exporter=standalone_cpp" >> ${outproc}.mg
    elif [ "${OUTBCK}" == "madgpu" ]; then # $SCRBCK=sycl and $OUTBCK=madgpu
      echo "output madevent ${outproc} ${helrecopt} --vector_size=32 --me_exporter=standalone_gpu" >> ${outproc}.mg
    else # $SCRBCK=sycl and $OUTBCK=sycl, cpp or gpu
      echo "output standalone_${OUTBCK} ${outproc}" >> ${outproc}.mg
    fi
    echo "--------------------------------------------------"
    cat ${outproc}.mg
    echo -e "--------------------------------------------------\n"
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
  if [ "${SCRBCK}" == "gridpack" ]; then
    outprocauto=${MG5AMC_HOME}/${outproc}/run_01_gridpack
    if ! $SCRDIR/untarGridpack.sh ${outprocauto}.tar.gz; then echo "ERROR! untarGridpack.sh failed"; exit 1; fi
  else
    outprocauto=${MG5AMC_HOME}/${outproc}
  fi
  cp -dpr ${MG5AMC_HOME}/${outproc}_log.txt ${outprocauto}/
  # Output directories: examples ee_mumu.sa for sycl, eemumu.auto for alpaka and gridpacks, eemumu.cpp or eemumu.gpu for cpp and gpu
  autosuffix=sa
  if [ "${SCRBCK}" == "gridpack" ]; then
    autosuffix=auto
  elif [ "${SCRBCK}" == "alpaka" ]; then
    autosuffix=auto
  elif [ "${OUTBCK}" == "cpp" ]; then
    autosuffix=cpp # no special suffix for the 311 branch any longer
  elif [ "${OUTBCK}" == "gpu" ]; then
    autosuffix=gpu # no special suffix for the 311 branch any longer
  elif [ "${OUTBCK}" == "madnovec" ] || [ "${OUTBCK}" == "madonly" ] || [ "${OUTBCK}" == "mad" ] || [ "${OUTBCK}" == "madcpp" ] || [ "${OUTBCK}" == "madgpu" ]; then
    autosuffix=${OUTBCK}
  fi
  # Replace the existing generated code in the output source code directory by the newly generated code and create a .BKP
  rm -rf ${OUTDIR}/${proc}.${autosuffix}.BKP
  if [ -d ${OUTDIR}/${proc}.${autosuffix} ]; then mv ${OUTDIR}/${proc}.${autosuffix} ${OUTDIR}/${proc}.${autosuffix}.BKP; fi
  cp -dpr ${outprocauto} ${OUTDIR}/${proc}.${autosuffix}
  echo -e "\nOutput source code has been copied to ${OUTDIR}/${proc}.${autosuffix}"
  # Fix build errors which arise because the autogenerated directories are not relocatable (see #400)
  if [ "${OUTBCK}" == "madnovec" ] || [ "${OUTBCK}" == "madonly" ] || [ "${OUTBCK}" == "mad" ] || [ "${OUTBCK}" == "madcpp" ] || [ "${OUTBCK}" == "madgpu" ]; then
    cat ${OUTDIR}/${proc}.${autosuffix}/Cards/me5_configuration.txt | sed 's/mg5_path/#mg5_path/' > ${OUTDIR}/${proc}.${autosuffix}/Cards/me5_configuration.txt.new
    \mv ${OUTDIR}/${proc}.${autosuffix}/Cards/me5_configuration.txt.new ${OUTDIR}/${proc}.${autosuffix}/Cards/me5_configuration.txt
  fi
  # Add a workaround for https://github.com/oliviermattelaer/mg5amc_test/issues/2
  if [ "${OUTBCK}" == "madnovec" ] || [ "${OUTBCK}" == "madonly" ] || [ "${OUTBCK}" == "mad" ] || [ "${OUTBCK}" == "madcpp" ] || [ "${OUTBCK}" == "madgpu" ]; then
    cat ${OUTDIR}/${proc}.${autosuffix}/Cards/ident_card.dat | head -3 > ${OUTDIR}/${proc}.${autosuffix}/Cards/ident_card.dat.new
    cat ${OUTDIR}/${proc}.${autosuffix}/Cards/ident_card.dat | tail -n+4 | sort >> ${OUTDIR}/${proc}.${autosuffix}/Cards/ident_card.dat.new
    \mv ${OUTDIR}/${proc}.${autosuffix}/Cards/ident_card.dat.new ${OUTDIR}/${proc}.${autosuffix}/Cards/ident_card.dat
  fi
  # Use strategy SDE=1 in multichannel mode (see #419)
  if [ "${OUTBCK}" == "mad" ]; then
    cat ${OUTDIR}/${proc}.${autosuffix}/Cards/run_card.dat | sed 's/2  = sde_strategy/1  = sde_strategy/' > ${OUTDIR}/${proc}.${autosuffix}/Cards/run_card.dat.new
    \mv ${OUTDIR}/${proc}.${autosuffix}/Cards/run_card.dat.new ${OUTDIR}/${proc}.${autosuffix}/Cards/run_card.dat
  fi
  # Additional patches for mad directory (integration of Fortran and sycl)
  # [NB: these patches are not applied to madnovec/madonly, which are meant as out-of-the-box references]
  if [ "${OUTBCK}" == "mad" ]; then
    $SCRDIR/patchMad.sh ${OUTDIR}/${proc}.${autosuffix} ${vecsize} ${dir_patches} ${PATCHLEVEL}
  fi
  # Compare the existing generated code to the newly generated code for the specific process
  pushd ${OUTDIR} >& /dev/null
  echo -e "\n+++ Compare old and new code generation log for $proc\n"
  ###if diff -c ${proc}.${autosuffix}.BKP/${outproc}_log.txt ${proc}.${autosuffix}; then echo "Old and new code generation logs are identical"; fi # context diff
  if diff ${proc}.${autosuffix}.BKP/${outproc}_log.txt ${proc}.${autosuffix}; then echo "Old and new code generation logs are identical"; fi # context diff
  echo -e "\n+++ Compare old and new generated code for $proc\n"
  if $SCRDIR/diffCode.sh ${BRIEF} -r -c ${proc}.${autosuffix}.BKP ${proc}.${autosuffix}; then echo "Old and new generated codes are identical"; else echo -e "\nWARNING! Old and new generated codes differ"; fi
  popd >& /dev/null
  # Compare the existing manually developed code to the newly generated code for the specific process
  if [ "${OUTBCK}" == "sycl" ] || [ "${OUTBCK}" == "gridpack" ]; then
    pushd ${OUTDIR} >& /dev/null
    echo -e "\n+++ Compare manually developed code to newly generated code for $proc\n"
    if $SCRDIR/diffCode.sh ${BRIEF} -r -c ${proc} ${proc}.${autosuffix}; then echo "Manual and generated codes are identical"; else echo -e "\nWARNING! Manual and generated codes differ"; fi
    popd >& /dev/null
  fi
  # Print a summary of the available code
  echo
  echo -e "Manually developed code is\n  ${OUTDIR}/${proc}"
  echo -e "Old generated code moved to\n  ${OUTDIR}/${proc}.${autosuffix}.BKP"
  echo -e "New generated code moved to\n  ${OUTDIR}/${proc}.${autosuffix}"
}

#--------------------------------------------------------------------------------------

function usage()
{
  # NB: Generate only one process at a time
  if [ "${SCRBCK}" == "gridpack" ]; then
    # NB: gridpack generation uses the 311 branch by default
    echo "Usage: $0 [--nobrief] [--nountaronly] [--nohelrec] <proc>"
  elif [ "${SCRBCK}" == "alpaka" ]; then
    # NB: alpaka generation uses the 311 branch by default
    echo "Usage: $0 [--nobrief] <proc>"
  else
    # NB: all options with $SCRBCK=sycl use the 311 branch by default and always disable helicity recycling
    echo "Usage: $0 [--nobrief] [--cpp|--gpu|--madnovec|--madonly|--mad|--madgpu] [--nopatch|--upstream] <proc>"
    echo "(NB: a --madcpp option also exists but code generation fails for it)"
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
  # Remove any *~ files in MG5AMC_HOME
  rm -rf $(find ${MG5AMC_HOME} -name '*~')
}

#--------------------------------------------------------------------------------------

# Script directory
SCRDIR=$(cd $(dirname $0); pwd)

# Output source code directory for the chosen backend (generated code will be created as a subdirectory of $OUTDIR)
OUTDIR=$(dirname $SCRDIR) # e.g. epochX/sycl if $SCRDIR=epochX/sycl/CODEGEN

# Script directory backend (sycl, gridpack or alpaka)
SCRBCK=$(basename $OUTDIR) # e.g. sycl if $OUTDIR=epochX/sycl

# Default output backend (in the sycl directory this can be changed using commad line options like --cpp, --gpu or --mad)
OUTBCK=$SCRBCK

# Default: brief diffs (use --nobrief to use full diffs)
BRIEF=--brief

# Default for gridpacks: untar gridpack.tar.gz but do not regenerate it (use --nountaronly to regenerate it)
UNTARONLY=1

# Default: apply all patches in patchMad.sh (--nopatch is ignored unless --mad is also specified)
PATCHLEVEL=

# Default for gridpacks: use helicity recycling (use --nohelrec to disable it)
# (export the value to the untarGridpack.sh script)
# Hardcoded for sycl and alpaka: disable helicity recycling (#400, #279) for the moment
if [ "${SCRBCK}" == "gridpack" ]; then export HELREC=1; else export HELREC=0; fi

# Process command line arguments (https://unix.stackexchange.com/a/258514)
for arg in "$@"; do
  shift
  if [ "$arg" == "-h" ] || [ "$arg" == "--help" ]; then
    usage
  elif [ "$arg" == "--nobrief" ]; then
    BRIEF=
  elif [ "$arg" == "--nopatch" ] && [ "${PATCHLEVEL}" == "" ]; then
    PATCHLEVEL=--nopatch
  elif [ "$arg" == "--upstream" ] && [ "${PATCHLEVEL}" == "" ]; then
    PATCHLEVEL=--upstream
  elif [ "$arg" == "--nountaronly" ] && [ "${SCRBCK}" == "gridpack" ]; then
    UNTARONLY=0
  elif [ "$arg" == "--nohelrec" ] && [ "${SCRBCK}" == "gridpack" ]; then
    export HELREC=0
  elif [ "$arg" == "--cpp" ] && [ "${SCRBCK}" == "sycl" ]; then
    export OUTBCK=${arg#--}
  elif [ "$arg" == "--gpu" ] && [ "${SCRBCK}" == "sycl" ]; then
    export OUTBCK=${arg#--}
  elif [ "$arg" == "--madnovec" ] && [ "${SCRBCK}" == "sycl" ]; then
    export OUTBCK=${arg#--}
  elif [ "$arg" == "--madonly" ] && [ "${SCRBCK}" == "sycl" ]; then
    export OUTBCK=${arg#--}
  elif [ "$arg" == "--mad" ] && [ "${SCRBCK}" == "sycl" ]; then
    export OUTBCK=${arg#--}
  elif [ "$arg" == "--madcpp" ] && [ "${SCRBCK}" == "sycl" ]; then
    export OUTBCK=${arg#--}
  elif [ "$arg" == "--madgpu" ] && [ "${SCRBCK}" == "sycl" ]; then
    export OUTBCK=${arg#--}
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
echo "SCRBCK=${SCRBCK} (uppercase=${SCRBCK^^})"
echo "OUTBCK=${OUTBCK}"

echo "BRIEF=${BRIEF}"
###echo "procs=${procs}"
echo "proc=${proc}"

# Make sure that python3 is installed
if ! python3 --version >& /dev/null; then echo "ERROR! python3 is not installed"; exit 1; fi

# Make sure that $MG5AMC_HOME exists
dir_patches=PROD
branch_patches=$(cat $SCRDIR/MG5aMC_patches/${dir_patches}/branch.GIT)
commit_patches=$(cat $SCRDIR/MG5aMC_patches/${dir_patches}/commit.GIT) # e.g. <commit> or <branch>
if [ "$MG5AMC_HOME" == "" ]; then
  echo "ERROR! MG5AMC_HOME is not defined"
  echo -e "To download MG5AMC please run\n  git clone git@github.com:mg5amcnlo/mg5amcnlo.git\n  cd mg5amcnlo; git checkout ${branch_patches}; git reset --hard ${commit_patches}"
  exit 1
fi
echo -e "\nDefault MG5AMC_HOME=$MG5AMC_HOME on $(hostname)\n"
if [ ! -d $MG5AMC_HOME ]; then
  echo "ERROR! Directory $MG5AMC_HOME does not exist"
  echo -e "To download MG5AMC please run\n  git clone git@github.com:mg5amcnlo/mg5amcnlo.git\n  cd mg5amcnlo; git checkout ${branch_patches}; git reset --hard ${commit_patches}"
  exit 1
fi

# Make sure that $ALPAKA_ROOT and $CUPLA_ROOT exist if alpaka is used
if [ "${SCRBCK}" == "alpaka" ]; then
  if [ "$ALPAKA_ROOT" == "" ]; then
    echo "ERROR! ALPAKA_ROOT is not defined"
    echo "To download ALPAKA please run 'git clone -b 0.8.0 https://github.com/alpaka-group/alpaka.git'"
    exit 1
  fi
  echo -e "Using ALPAKA_ROOT=$ALPAKA_ROOT on $(hostname)\n"
  if [ ! -d $ALPAKA_ROOT ]; then echo "ERROR! Directory $ALPAKA_ROOT does not exist"; exit 1; fi
  if [ "$CUPLA_ROOT" == "" ]; then
    echo "ERROR! CUPLA_ROOT is not defined"
    echo "To download CUPLA please run 'git clone -b 0.3.0 https://github.com/alpaka-group/cupla.git'"
    exit 1
  fi
  echo -e "Using CUPLA_ROOT=$CUPLA_ROOT on $(hostname)\n"
  if [ ! -d $CUPLA_ROOT ]; then echo "ERROR! Directory $CUPLA_ROOT does not exist"; exit 1; fi
fi

# Check that MG5aMC uses the git 311 branch (default for all of sycl, alpaka, gridpack)
# Revert MG5aMC to the appropriate git commit
if ! git --version >& /dev/null; then
  echo -e "ERROR! git is not installed: cannot retrieve git properties of MG5aMC_HOME\n"; exit 1
fi
echo -e "Using $(git --version)"
cd ${MG5AMC_HOME}
echo -e "Retrieving git information about MG5AMC_HOME"
if ! git log -n1 >& /dev/null; then
  echo -e "ERROR! MG5AMC_HOME is not a git clone\n"; exit 1
fi
echo -e "MG5AMC patches in this plugin refer to git branch '${branch_patches}'"
echo -e "Reset MG5AMC_HOME to git commit '${branch_patches}'"
if ! git reset --hard ${branch_patches}; then
  echo -e "ERROR! 'git reset --hard ${branch_patches}' failed\n"; exit 1
fi
echo -e "Check out branch ${branch_patches} in MG5AMC_HOME"
if ! git checkout ${branch_patches}; then
  echo -e "ERROR! 'git checkout ${branch_patches}' failed\n"; exit 1
fi
branch_mg5amc=$(git branch --no-color | \grep ^* | awk '{print $2}')
echo -e "Current git branch of MG5AMC_HOME is '${branch_mg5amc}'"
if [ "${branch_patches}" != "${branch_mg5amc}" ]; then echo -e "\nERROR! git branch mismatch!"; exit 1; fi
commit_patches2=$(git log --oneline -n1 ${commit_patches} | awk '{print $1}') # translate to <commit>
if [ "${commit_patches2}" == "${commit_patches}" ]; then
  echo -e "MG5AMC patches in this plugin refer to git commit '${commit_patches}'"
else
  echo -e "MG5AMC patches in this plugin refer to git commit '${commit_patches}' (i.e. '${commit_patches2}')"
fi  
echo -e "Reset MG5AMC_HOME to git commit '${commit_patches}'"
if ! git reset --hard ${commit_patches}; then
  echo -e "ERROR! 'git reset --hard ${commit_patches}' failed\n"; exit 1
fi
commit_mg5amc=$(git log --oneline -n1 | awk '{print $1}')
echo -e "Current git commit of MG5AMC_HOME is '${commit_mg5amc}'"
if [ "${commit_patches2}" != "${commit_mg5amc}" ]; then echo -e "\nERROR! git commit mismatch!"; exit 1; fi
cd - > /dev/null

# Copy MG5AMC ad-hoc patches if any (unless --upstream is specified)
if [ "${PATCHLEVEL}" != "--upstream" ] && [ "${SCRBCK}" == "sycl" ]; then
  patches=$(cd $SCRDIR/MG5aMC_patches/${dir_patches}; find . -mindepth 2 -type f)
  echo -e "Copy MG5aMC_patches/${dir_patches} patches..."
  for patch in $patches; do
    ###echo "DEBUG: $patch"
    patch=${patch#./} # strip leading "./"
    if [ "${patch}" != "${patch%~}" ]; then continue; fi # skip *~ files
    ###if [ "${patch}" != "${patch%.GIT}" ]; then continue; fi # skip *.GIT files
    ###if [ "${patch}" != "${patch%.py}" ]; then continue; fi # skip *.py files
    ###if [ "${patch}" != "${patch#patch.}" ]; then continue; fi # skip patch.* files
    echo cp -dpr $SCRDIR/MG5aMC_patches/${dir_patches}/$patch $MG5AMC_HOME/$patch
    cp -dpr $SCRDIR/MG5aMC_patches/${dir_patches}/$patch $MG5AMC_HOME/$patch
  done
  echo -e "Copy MG5aMC_patches/${dir_patches} patches... done\n"
fi

# Clean up before code generation
cleanup_MG5AMC_HOME

# Print differences in MG5AMC with respect to git after copying ad-hoc patches
cd ${MG5AMC_HOME}
echo -e "\n***************** Differences to the current git commit ${commit_patches} [START]"
###if [ "$(git diff)" == "" ]; then echo -e "[No differences]"; else git diff; fi
if [ "$(git diff)" == "" ]; then echo -e "[No differences]"; else git diff --name-status; fi
echo -e "***************** Differences to the current git commit ${commit_patches} [END]"
cd - > /dev/null

# Copy the new plugin to MG5AMC_HOME (if the script directory backend is sycl or alpaka)
if [ "${SCRBCK}" == "sycl" ]; then
  if [ "${OUTBCK}" == "no-path-to-this-statement" ]; then
    echo -e "\nWARNING! '${OUTBCK}' mode selected: do not copy the sycl plugin (workaround for #341)"
  else # currently succeeds also for madcpp and madgpu (#341 has been fixed)
    echo -e "\nINFO! '${OUTBCK}' mode selected: copy the sycl plugin\n"
    cp -dpr ${SCRDIR}/PLUGIN/${SCRBCK^^}_SA_OUTPUT ${MG5AMC_HOME}/PLUGIN/
    ls -l ${MG5AMC_HOME}/PLUGIN
  fi
elif [ "${SCRBCK}" == "alpaka" ]; then
  cp -dpr ${SCRDIR}/PLUGIN/${SCRBCK^^}_SYCL_SA_OUTPUT ${MG5AMC_HOME}/PLUGIN/
  ls -l ${MG5AMC_HOME}/PLUGIN
fi

# For gridpacks, use separate output directories for MG 29x and MG 3xx
if [ "${SCRBCK}" == "gridpack" ]; then
  if [ "${HELREC}" == "0" ]; then
    OUTDIR=${OUTDIR}/3xx_nohelrec
  else
    OUTDIR=${OUTDIR}/3xx
  fi
  echo "OUTDIR=${OUTDIR} (redefined)"
fi

# Generate the chosen process (this will always replace the existing code directory and create a .BKP)
codeGenAndDiff $proc

# Clean up after code generation
cleanup_MG5AMC_HOME

# Check formatting in the auto-generated code
#if [ "${OUTBCK}" == "sycl" ]; then
#  echo -e "\n+++ Check code formatting in newly generated code ${proc}.sa\n"
#  if ! $SCRDIR/checkFormatting.sh -q -q ${proc}.sa; then
#    echo "ERROR! Auto-generated code does not respect formatting policies"
#    exit 1
#  fi
#elif [ "${OUTBCK}" == "mad" ]; then
#  echo -e "\n+++ Check code formatting in newly generated code ${proc}.mad\n"
#  if ! $SCRDIR/checkFormatting.sh -q -q ${proc}.mad; then
#    echo "ERROR! Auto-generated code does not respect formatting policies"
#    exit 1
#  fi
#fi
