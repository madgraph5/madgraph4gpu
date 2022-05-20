#!/bin/bash

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
  # Vector size for mad/madonly meexporter (nb_page_max)
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
    echo "${cmd}" >> ${outproc}.mg
    if [ "${SCRBCK}" == "gridpack" ]; then # $SCRBCK=$OUTBCK=gridpack
      echo "output ${outproc} ${helrecopt}" >> ${outproc}.mg
      ###echo "!cp -dpr ${outproc} ${outproc}_prelaunch" >> ${outproc}.mg
      echo "launch" >> ${outproc}.mg
      echo "set gridpack True" >> ${outproc}.mg
      echo "set ebeam1 750" >> ${outproc}.mg
      echo "set ebeam2 750" >> ${outproc}.mg
    elif [ "${SCRBCK}" == "alpaka" ]; then # $SCRBCK=$OUTBCK=alpaka
      echo "output standalone_${SCRBCK}_cudacpp ${outproc}" >> ${outproc}.mg
    elif [ "${OUTBCK}" == "madonly" ]; then # $SCRBCK=cudacpp and $OUTBCK=madonly
      echo "output madevent ${outproc} ${helrecopt} --vector_size=${vecsize}" >> ${outproc}.mg
    elif [ "${OUTBCK}" == "mad" ]; then # $SCRBCK=cudacpp and $OUTBCK=mad
      echo "output madevent ${outproc} ${helrecopt} --vector_size=${vecsize} --me_exporter=standalone_cudacpp" >> ${outproc}.mg
    elif [ "${OUTBCK}" == "madcpp" ]; then # $SCRBCK=cudacpp and $OUTBCK=madcpp
      echo "output madevent ${outproc} ${helrecopt} --vector_size=32 --me_exporter=standalone_cpp" >> ${outproc}.mg
    elif [ "${OUTBCK}" == "madgpu" ]; then # $SCRBCK=cudacpp and $OUTBCK=madgpu
      echo "output madevent ${outproc} ${helrecopt} --vector_size=32 --me_exporter=standalone_gpu" >> ${outproc}.mg
    else # $SCRBCK=cudacpp and $OUTBCK=cudacpp, cpp or gpu
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
  # Output directories: examples ee_mumu.auto for cudacpp and gridpacks, eemumu.cpp270 or eemumu.gpu270 for cpp
  autosuffix=auto
  if [ "${OUTBCK}" == "cpp" ]; then
    if [ "$use270" == "0" ]; then autosuffix=cpp270; else autosuffix=cpp311; fi
  elif [ "${OUTBCK}" == "gpu" ]; then
    if [ "$use270" == "0" ]; then autosuffix=gpu270; else autosuffix=gpu311; fi
  elif [ "${OUTBCK}" == "madonly" ] || [ "${OUTBCK}" == "mad" ] || [ "${OUTBCK}" == "madcpp" ] || [ "${OUTBCK}" == "madgpu" ]; then
    autosuffix=${OUTBCK}
  fi
  # Replace the existing generated code in the output source code directory by the newly generated code and create a .BKP
  rm -rf ${OUTDIR}/${proc}.${autosuffix}.BKP
  if [ -d ${OUTDIR}/${proc}.${autosuffix} ]; then mv ${OUTDIR}/${proc}.${autosuffix} ${OUTDIR}/${proc}.${autosuffix}.BKP; fi
  cp -dpr ${outprocauto} ${OUTDIR}/${proc}.${autosuffix}
  echo -e "\nOutput source code has been copied to ${OUTDIR}/${proc}.${autosuffix}"
  # Fix build errors which arise because the autogenerated directories are not relocatable (see #400)
  if [ "${OUTBCK}" == "madonly" ] || [ "${OUTBCK}" == "mad" ] || [ "${OUTBCK}" == "madcpp" ] || [ "${OUTBCK}" == "madgpu" ]; then
    cat ${OUTDIR}/${proc}.${autosuffix}/Cards/me5_configuration.txt | sed 's/mg5_path/#mg5_path/' > ${OUTDIR}/${proc}.${autosuffix}/Cards/me5_configuration.txt.new
    \mv ${OUTDIR}/${proc}.${autosuffix}/Cards/me5_configuration.txt.new ${OUTDIR}/${proc}.${autosuffix}/Cards/me5_configuration.txt
  fi
  # Add a workaround for https://github.com/oliviermattelaer/mg5amc_test/issues/2
  if [ "${OUTBCK}" == "madonly" ] || [ "${OUTBCK}" == "mad" ] || [ "${OUTBCK}" == "madcpp" ] || [ "${OUTBCK}" == "madgpu" ]; then
    cat ${OUTDIR}/${proc}.${autosuffix}/Cards/ident_card.dat | head -3 > ${OUTDIR}/${proc}.${autosuffix}/Cards/ident_card.dat.new
    cat ${OUTDIR}/${proc}.${autosuffix}/Cards/ident_card.dat | tail -n+4 | sort >> ${OUTDIR}/${proc}.${autosuffix}/Cards/ident_card.dat.new
    \mv ${OUTDIR}/${proc}.${autosuffix}/Cards/ident_card.dat.new ${OUTDIR}/${proc}.${autosuffix}/Cards/ident_card.dat
  fi
  # Additional patches for mad directory (integration of Fortran and cudacpp)
  # [NB: NEW! these patches are no longer applied to madonly, which is now meant as an out-of-the-box reference]
  ###if [ "${OUTBCK}" == "madonly" ] || [ "${OUTBCK}" == "mad" ]; then
  if [ "${OUTBCK}" == "mad" ]; then
    $SCRDIR/patchMad.sh ${OUTDIR}/${proc}.${autosuffix} ${vecsize} ${NOPATCH}
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
  if [ "${OUTBCK}" == "cudacpp" ] || [ "${OUTBCK}" == "gridpack" ]; then
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
    # NB: gridpack generation has been tested only agains the 270 branch so far
    echo "Usage: $0 [--nobrief] [--nountaronly] [--nohelrec] <proc>"
  elif [ "${SCRBCK}" == "alpaka" ]; then
    # NB: alpaka generation has been tested only agains the 270 branch so far
    echo "Usage: $0 [--nobrief] <proc>"
  else
    # NB: all options with $SCRBCK=cudacpp use the 311 branch by default and always disable helicity recycling
    echo "Usage: $0 [--nobrief] [--cpp|--gpu|--madonly|--mad|--madcpp|--madgpu] [--270] [--nopatch] <proc>"
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

# Output source code directory for the chosen backend (generated code will be created as a subdirectory of $OUTDIR)
OUTDIR=$(dirname $SCRDIR) # e.g. epochX/cudacpp if $SCRDIR=epochX/cudacpp/CODEGEN

# Script directory backend (cudacpp, gridpack or alpaka)
SCRBCK=$(basename $OUTDIR) # e.g. cudacpp if $OUTDIR=epochX/cudacpp

# Default output backend (in the cudacpp directory this can be changed using commad line options like --cpp, --gpu or --mad)
OUTBCK=$SCRBCK

# Default: brief diffs (use --nobrief to use full diffs)
BRIEF=--brief

# Default: use the 311 MG5aMC branch (except for alpaka and gridpack)
use270=0
if [ "${SCRBCK}" == "alpaka" ] || [ "${SCRBCK}" == "gridpack" ]; then use270=1; fi

# Default for gridpacks: untar gridpack.tar.gz but do not regenerate it (use --nountaronly to regenerate it)
UNTARONLY=1

# Default: apply all patches in patchMad.sh (this is ignored unless --mad is also specified)
NOPATCH=

# Default for gridpacks: use helicity recycling (use --nohelrec to disable it)
# (export the value to the untarGridpack.sh script)
# Hardcoded for cudacpp and alpaka: disable helicity recycling (#400, #279) for the moment
if [ "${SCRBCK}" == "gridpack" ]; then export HELREC=1; else export HELREC=0; fi

# Process command line arguments (https://unix.stackexchange.com/a/258514)
for arg in "$@"; do
  shift
  if [ "$arg" == "-h" ] || [ "$arg" == "--help" ]; then
    usage
  elif [ "$arg" == "--nobrief" ]; then
    BRIEF=
  elif [ "$arg" == "--nopatch" ]; then
    NOPATCH=--nopatch
  elif [ "$arg" == "--270" ]; then
    use270=1
  elif [ "$arg" == "--nountaronly" ] && [ "${SCRBCK}" == "gridpack" ]; then
    UNTARONLY=0
  elif [ "$arg" == "--nohelrec" ] && [ "${SCRBCK}" == "gridpack" ]; then
    export HELREC=0
  elif [ "$arg" == "--cpp" ] && [ "${SCRBCK}" == "cudacpp" ]; then
    export OUTBCK=${arg#--}
  elif [ "$arg" == "--gpu" ] && [ "${SCRBCK}" == "cudacpp" ]; then
    export OUTBCK=${arg#--}
  elif [ "$arg" == "--madonly" ] && [ "${SCRBCK}" == "cudacpp" ]; then
    export OUTBCK=${arg#--}
  elif [ "$arg" == "--mad" ] && [ "${SCRBCK}" == "cudacpp" ]; then
    export OUTBCK=${arg#--}
  elif [ "$arg" == "--madcpp" ] && [ "${SCRBCK}" == "cudacpp" ]; then
    export OUTBCK=${arg#--}
  elif [ "$arg" == "--madgpu" ] && [ "${SCRBCK}" == "cudacpp" ]; then
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
mg5amc270=2.7.0_gpu
mg5amc311=3.1.1_lo_vectorization
mg5amcBrn=${mg5amc311}
if [ "${OUTBCK}" == "alpaka" ]; then
  revno_patches=370
else
  revno_patches=$(cat $SCRDIR/MG5aMC_patches/${mg5amcBrn}/revision.BZR)
fi
if [ "$MG5AMC_HOME" == "" ]; then
  echo "ERROR! MG5AMC_HOME is not defined"
  echo -e "To download MG5AMC please run\n  bzr branch lp:~maddevelopers/mg5amcnlo/${mg5amcBrn} -r ${revno_patches}"
  exit 1
fi
echo -e "\nDefault MG5AMC_HOME=$MG5AMC_HOME on $(hostname)\n"
if [ ! -d $MG5AMC_HOME ]; then
  echo "ERROR! Directory $MG5AMC_HOME does not exist"
  echo -e "To download MG5AMC please run\n  bzr branch lp:~maddevelopers/mg5amcnlo/${mg5amcBrn} -r ${revno_patches}"
  exit 1
fi
if [ "$(basename ${MG5AMC_HOME})" != "${mg5amcBrn}" ]; then
  echo "ERROR! MG5AMC_HOME basename is not ${mg5amcBrn}"
  exit 1
fi

# Redefine MG5AMC_HOME to use the 270 branch if required
if [ "$use270" == "1" ]; then
  if [ "${OUTBCK}" != "cpp" ] && [ "${OUTBCK}" != "gpu" ]; then
    echo "ERROR! Option --270 is only supported in --cpu or --gpu mode"
    exit 1
  fi
  mg5amcBrn=${mg5amc270}
  export MG5AMC_HOME=$(dirname ${MG5AMC_HOME})/${mg5amcBrn}
  echo -e "Using non-default MG5AMC_HOME=$MG5AMC_HOME on $(hostname)\n"
  if [ ! -d $MG5AMC_HOME ]; then echo "ERROR! Directory $MG5AMC_HOME does not exist"; exit 1; fi
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

# Print MG5amc bazaar info if any
# Revert to the appropriate bazaar revision number
# (NB! 'bzr revert' does not change the output of 'bzr revno': it is NOT like 'git reset --hard'!)
# (See the comments in https://stackoverflow.com/a/37488587)
if bzr --version >& /dev/null; then
  echo -e "Using $(bzr --version | head -1)"
  echo -e "Retrieving bzr information about MG5AMC_HOME"
  if bzr info ${MG5AMC_HOME} > /dev/null; then
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
if [ "${SCRBCK}" == "cudacpp" ]; then
  ###patches=$(cd $SCRDIR/MG5aMC_patches/${mg5amcBrn}; find . -type f -name '*.py')
  patches=$(cd $SCRDIR/MG5aMC_patches/${mg5amcBrn}; find . -type f ! -name '*.BZR')
  echo -e "Copy MG5aMC_patches/${mg5amcBrn} patches..."
  for patch in $patches; do
    patch=${patch#./}
    echo cp -dpr $SCRDIR/MG5aMC_patches/${mg5amcBrn}/$patch $MG5AMC_HOME/$patch
    cp -dpr $SCRDIR/MG5aMC_patches/${mg5amcBrn}/$patch $MG5AMC_HOME/$patch
  done
  echo -e "Copy MG5aMC_patches/${mg5amcBrn} patches... done\n"
fi

# Clean up before code generation
cleanup_MG5AMC_HOME

# Print MG5amc bazaar info if any
if bzr --version >& /dev/null; then
  if bzr info ${MG5AMC_HOME} 2> /dev/null | grep parent; then
    echo -e "\n***************** Differences to the current bzr revno ${revno_patches} [START]"
    if bzr diff ${MG5AMC_HOME} -r ${revno_patches}; then echo -e "[No differences]"; fi
    echo -e "***************** Differences to the current bzr revno ${revno_patches} [END]"
  fi
fi

# Copy the new plugin to MG5AMC_HOME (if the script directory backend is cudacpp or alpaka)
if [ "${SCRBCK}" == "cudacpp" ]; then
  if [ "${OUTBCK}" == "no-path-to-this-statement" ]; then
    echo -e "\nWARNING! '${OUTBCK}' mode selected: do not copy the cudacpp plugin (workaround for #341)"
  else # currently succeeds also for madcpp and madgpu (#341 has been fixed)
    echo -e "\nINFO! '${OUTBCK}' mode selected: copy the cudacpp plugin\n"
    cp -dpr ${SCRDIR}/PLUGIN/${SCRBCK^^}_SA_OUTPUT ${MG5AMC_HOME}/PLUGIN/
    ls -l ${MG5AMC_HOME}/PLUGIN
  fi
elif [ "${SCRBCK}" == "alpaka" ]; then
  cp -dpr ${SCRDIR}/PLUGIN/${SCRBCK^^}_CUDACPP_SA_OUTPUT ${MG5AMC_HOME}/PLUGIN/
  ls -l ${MG5AMC_HOME}/PLUGIN
fi

# For gridpacks, use separate output directories for MG 28x and MG 29x
if [ "${SCRBCK}" == "gridpack" ]; then
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

# Check formatting in the auto-generated code
if [ "${OUTBCK}" == "cudacpp" ]; then
  echo -e "\n+++ Check code formatting in newly generated code for $proc\n"
  if ! $SCRDIR/checkFormatting.sh -q -q ${proc}.auto; then
    echo "ERROR! Auto-generated code does not respect formatting policies"
    exit 1
  fi
fi
