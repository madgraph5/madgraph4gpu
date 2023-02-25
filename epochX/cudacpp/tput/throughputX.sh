#!/bin/bash
# Copyright (C) 2020-2023 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: A. Valassi (Apr 2021) for the MG5aMC CUDACPP plugin.

set +x # not verbose
set -e # fail on error

scrdir=$(cd $(dirname $0); pwd)
bckend=$(basename $(cd $scrdir; cd ..; pwd)) # cudacpp or alpaka
topdir=$(cd $scrdir; cd ../../..; pwd)

function usage()
{
  echo "Usage: $0 <processes [-eemumu][-ggtt][-ggttg][-ggttgg][-ggttggg][-ggttgggg][-gqttq][-heftggh]> [-nocpp|[-avxall][-nocuda][-noneonly][-sse4only][-avx2only][-512yonly][-512zonly]] [-sa] [-noalpaka] [-flt|-fltonly|-mix|-mixonly] [-inl|-inlonly] [-hrd|-hrdonly] [-common|-curhst] [-rmbhst|-bridge] [-omp] [-makeonly|-makeclean|-makecleanonly|-dryrun] [-makej] [-3a3b] [-div] [-req] [-detailed] [-gtest] [-nofpe] [-v] [-dlp <dyld_library_path>]"
  exit 1
}

##########################################################################
# PART 0 - decode command line arguments
##########################################################################

procs=
eemumu=0
ggtt=0
ggttg=0
ggttgg=0
ggttggg=0
ggttgggg=0
gqttq=0
heftggh=0

suffs=".mad/"

omp=0
simds=
cpp=1
cuda=1
alpaka=1

fptypes="d"
helinls="0"
hrdcods="0"
rndgen=""
rmbsam=""

maketype=
makej=

ab3=0
div=0
req=0
detailed=0
gtest=0
nofpe=0
verbose=0

dlp=

# Optional hack to build only the cudacpp plugin (without building the madevent code) in .mad directories
makef=
###makef="-f Makefile"

# Workaround to allow 'make avxall' when '-avxall' is specified #536
simdsall="none sse4 avx2 512y 512z"

if [ "$bckend" != "alpaka" ]; then alpaka=0; fi # alpaka mode is only available in the alpaka directory

while [ "$1" != "" ]; do
  if [ "$1" == "-eemumu" ]; then
    if [ "$eemumu" == "0" ]; then procs+=${procs:+ }$1; fi
    eemumu=1
    shift
  elif [ "$1" == "-ggtt" ]; then
    if [ "$ggtt" == "0" ]; then procs+=${procs:+ }$1; fi
    ggtt=1
    shift
  elif [ "$1" == "-ggttg" ]; then
    if [ "$ggttg" == "0" ]; then procs+=${procs:+ }$1; fi
    ggttg=1
    shift
  elif [ "$1" == "-ggttgg" ]; then
    if [ "$ggttgg" == "0" ]; then procs+=${procs:+ }$1; fi
    ggttgg=1
    shift
  elif [ "$1" == "-ggttggg" ]; then
    if [ "$ggttggg" == "0" ]; then procs+=${procs:+ }$1; fi
    ggttggg=1
    shift
  elif [ "$1" == "-ggttgggg" ]; then
    if [ "$ggttgggg" == "0" ]; then procs+=${procs:+ }$1; fi
    ggttgggg=1
  elif [ "$1" == "-gqttq" ]; then
    if [ "$gqttq" == "0" ]; then procs+=${procs:+ }$1; fi
    gqttq=1
    shift
  elif [ "$1" == "-heftggh" ]; then
    if [ "$heftggh" == "0" ]; then procs+=${procs:+ }$1; fi
    heftggh=1
    shift
  elif [ "$1" == "-sa" ]; then
    suffs=".sa/"
    shift
  elif [ "$1" == "-nocuda" ]; then
    if [ "${cpp}" == "0" ]; then echo "ERROR! Options -nocuda and -nocpp are incompatible"; usage; fi
    cuda=0
    shift
  elif [ "$1" == "-nocpp" ]; then
    if [ "${omp}" == "1" ]; then echo "ERROR! Options -omp and -nocpp are incompatible"; usage; fi
    if [ "${simds}" != "" ]; then echo "ERROR! Incompatible option $1: SIMDs are already defined as '$simds'"; usage; fi
    if [ "${cuda}" == "0" ]; then echo "ERROR! Options -nocuda and -nocpp are incompatible"; usage; fi
    cpp=0
    shift
  elif [ "$1" == "-omp" ]; then
    if [ "${cpp}" == "0" ]; then echo "ERROR! Options -omp and -nocpp are incompatible"; usage; fi
    omp=1
    shift
  elif [ "$1" == "-avxall" ]; then
    if [ "${cpp}" == "0" ]; then echo "ERROR! Options $1 and -nocpp are incompatible"; usage; fi
    if [ "${simds}" != "" ]; then echo "ERROR! Incompatible option $1: SIMDs are already defined as '$simds'"; usage; fi
    simds="${simdsall}"
    shift
  elif [ "$1" == "-noneonly" ]; then
    if [ "${cpp}" == "0" ]; then echo "ERROR! Options $1 and -nocpp are incompatible"; usage; fi
    if [ "${simds}" != "" ]; then echo "ERROR! Incompatible option $1: SIMDs are already defined as '$simds'"; usage; fi
    simds="none"
    cuda=0; echo "WARNING! Option $1 implies -nocuda"
    shift
  elif [ "$1" == "-sse4only" ]; then
    if [ "${cpp}" == "0" ]; then echo "ERROR! Options $1 and -nocpp are incompatible"; usage; fi
    if [ "${simds}" != "" ]; then echo "ERROR! Incompatible option $1: SIMDs are already defined as '$simds'"; usage; fi
    simds="sse4"
    cuda=0; echo "WARNING! Option $1 implies -nocuda"
    shift
  elif [ "$1" == "-avx2only" ]; then
    if [ "${cpp}" == "0" ]; then echo "ERROR! Options $1 and -nocpp are incompatible"; usage; fi
    if [ "${simds}" != "" ]; then echo "ERROR! Incompatible option $1: SIMDs are already defined as '$simds'"; usage; fi
    simds="avx2"
    cuda=0; echo "WARNING! Option $1 implies -nocuda"
    shift
  elif [ "$1" == "-512yonly" ]; then
    if [ "${cpp}" == "0" ]; then echo "ERROR! Options $1 and -nocpp are incompatible"; usage; fi
    if [ "${simds}" != "" ]; then echo "ERROR! Incompatible option $1: SIMDs are already defined as '$simds'"; usage; fi
    simds="512y"
    cuda=0; echo "WARNING! Option $1 implies -nocuda"
    shift
  elif [ "$1" == "-512zonly" ]; then
    if [ "${cpp}" == "0" ]; then echo "ERROR! Options $1 and -nocpp are incompatible"; usage; fi
    if [ "${simds}" != "" ]; then echo "ERROR! Incompatible option $1: SIMDs are already defined as '$simds'"; usage; fi
    simds="512z"
    cuda=0; echo "WARNING! Option $1 implies -nocuda"
    shift
  elif [ "$1" == "-noalpaka" ]; then
    alpaka=0
    shift
  elif [ "$1" == "-flt" ]; then
    if [ "${fptypes}" != "d" ] && [ "${fptypes}" != "d f" ]; then echo "ERROR! Options -flt, -fltonly, -mix and -mixonly are incompatible"; usage; fi
    fptypes="d f"
    shift
  elif [ "$1" == "-fltonly" ]; then
    if [ "${fptypes}" != "d" ] && [ "${fptypes}" != "f" ]; then echo "ERROR! Options -flt, -fltonly, -mix and -mixonly are incompatible"; usage; fi
    fptypes="f"
    shift
  elif [ "$1" == "-mix" ]; then
    if [ "${fptypes}" != "d" ] && [ "${fptypes}" != "d f m" ]; then echo "ERROR! Options -flt, -fltonly, -mix and -mixonly are incompatible"; usage; fi
    fptypes="d f m"
    shift
  elif [ "$1" == "-mixonly" ]; then
    if [ "${fptypes}" != "d" ] && [ "${fptypes}" != "m" ]; then echo "ERROR! Options -flt, -fltonly, -mix and -mixonly are incompatible"; usage; fi
    fptypes="m"
    shift
  elif [ "$1" == "-inl" ]; then
    if [ "${helinls}" == "1" ]; then echo "ERROR! Options -inl and -inlonly are incompatible"; usage; fi
    helinls="0 1"
    shift
  elif [ "$1" == "-inlonly" ]; then
    if [ "${helinls}" == "0 1" ]; then echo "ERROR! Options -inl and -inlonly are incompatible"; usage; fi
    helinls="1"
    shift
  elif [ "$1" == "-hrd" ]; then
    if [ "${hrdcods}" == "1" ]; then echo "ERROR! Options -hrd and -hrdonly are incompatible"; usage; fi
    hrdcods="0 1"
    shift
  elif [ "$1" == "-hrdonly" ]; then
    if [ "${hrdcods}" == "0 1" ]; then echo "ERROR! Options -hrd and -hrdonly are incompatible"; usage; fi
    hrdcods="1"
    shift
  elif [ "$1" == "-common" ]; then
    rndgen=" -${1}"
    shift
  elif [ "$1" == "-curhst" ]; then
    rndgen=" -${1}"
    shift
  elif [ "$1" == "-rmbhst" ]; then
    rmbsmp=" -${1}"
    shift
  elif [ "$1" == "-bridge" ]; then
    rmbsmp=" -${1}"
    shift
  elif [ "$1" == "-makeonly" ] || [ "$1" == "-makeclean" ] || [ "$1" == "-makecleanonly" ] || [ "$1" == "-dryrun" ]; then
    if [ "${maketype}" != "" ] && [ "${maketype}" != "$1" ]; then
      echo "ERROR! Options -makeonly, -makeclean, -makecleanonly and -dryrun are incompatible"; usage
    fi
    maketype="$1"
    shift
  elif [ "$1" == "-makej" ]; then
    ###makej=-j
    makej=-j5 # limit build parallelism to avoid "cudafe++ died due to signal 9" (#639)
    shift
  elif [ "$1" == "-3a3b" ]; then
    ab3=1
    shift
  elif [ "$1" == "-div" ]; then
    div=1
    shift
  elif [ "$1" == "-req" ]; then
    req=1
    shift
  elif [ "$1" == "-detailed" ]; then
    detailed=1
    shift
  elif [ "$1" == "-gtest" ]; then
    # For simplicity a gtest runTest.exe is executed for each build where check.exe is executed
    if [ "${cpp}" == "0" ]; then echo "ERROR! Options -gtest and -nocpp are incompatible"; usage; fi
    gtest=1
    shift
  elif [ "$1" == "-nofpe" ]; then
    nofpe=1
    shift
  elif [ "$1" == "-v" ]; then
    verbose=1
    shift
  elif [ "$1" == "-dlp" ] && [ "$2" != "" ]; then
    dlp="$2"
    shift
    shift
  else
    echo "ERROR! Invalid option '$1'"
    usage
  fi
done
###exit 1

# Workaround for MacOS SIP (SystemIntegrity Protection): set DYLD_LIBRARY_PATH In subprocesses
if [ "${dlp}" != "" ]; then
  echo "export DYLD_LIBRARY_PATH=$dlp"
  export DYLD_LIBRARY_PATH=$dlp
fi

# Enable FPEs in check.exe by default (see #733)
if [ "${nofpe}" == "0" ]; then
  echo "export CUDACPP_RUNTIME_ENABLEFPE=on"
  export CUDACPP_RUNTIME_ENABLEFPE=on
else
  echo "unset CUDACPP_RUNTIME_ENABLEFPE"
  unset CUDACPP_RUNTIME_ENABLEFPE
fi

# Check that at least one process has been selected
if [ "${eemumu}" == "0" ] && [ "${ggtt}" == "0" ] && [ "${ggttg}" == "0" ] && [ "${ggttgg}" == "0" ] && [ "${ggttggg}" == "0" ] && [ "${ggttgggg}" == "0" ] && [ "${gqttq}" == "0" ] && [ "${heftggh}" == "0" ]; then usage; fi

# Check that heftggh does not run in .mad mode
if [ "${heftggh}" == "1" ] && [ "${suffs/.mad\/}" != "${suffs}" ]; then
  echo "ERROR! Invalid option -heftggh for .mad directories"; exit 1
fi

# Define the default simds if none are defined
if [ "${simds}" == "" ]; then simds="none 512y"; fi

# Use only the .auto process directories in the alpaka directory
if [ "$bckend" == "alpaka" ]; then
  echo "WARNING! alpaka directory: using .auto process directories only"
  suffs=".auto/"
fi

# Use only HRDCOD=0 in the alpaka directory (old epochX-golden2 code base)
if [ "$bckend" == "alpaka" ]; then
  echo "WARNING! alpaka directory: using HRDCOD=0 only"
  hrdcods="0"
fi

# Check whether Alpaka should and can be run
if [ "${alpaka}" == "1" ]; then
  if [ "${CUPLA_ROOT}" == "" ]; then echo "ERROR! CUPLA_ROOT is not set!"; exit 1; fi
  echo "CUPLA_ROOT=$CUPLA_ROOT"
  if [ ! -d "${CUPLA_ROOT}" ]; then echo "ERROR! $CUPLA_ROOT does not exist!"; exit 1; fi
  if [ "${ALPAKA_ROOT}" == "" ]; then echo "ERROR! ALPAKA_ROOT is not set!"; exit 1; fi
  echo "ALPAKA_ROOT=$ALPAKA_ROOT"
  if [ ! -d "${ALPAKA_ROOT}" ]; then echo "ERROR! $ALPAKA_ROOT does not exist!"; exit 1; fi
  if [ "${BOOSTINC}" == "" ]; then echo "ERROR! BOOSTINC is not set!"; exit 1; fi
  echo "BOOSTINC=$BOOSTINC"
  if [ ! -d "${BOOSTINC}" ]; then echo "ERROR! $BOOSTINC does not exist!"; exit 1; fi  
else
  export CUPLA_ROOT=none
fi

# Determine the O/S and the processor architecture for late decisions
unames=$(uname -s)
unamep=$(uname -p)

# Determine the working directory below topdir based on suff, bckend and proc
function showdir()
{
  if [ "${suff}" == ".mad/" ]; then
    if [ "${proc}" == "-eemumu" ]; then 
      dir=$topdir/epochX/${bckend}/ee_mumu${suff}SubProcesses/P1_epem_mupmum
    elif [ "${proc}" == "-ggtt" ]; then 
      dir=$topdir/epochX/${bckend}/gg_tt${suff}SubProcesses/P1_gg_ttx
    elif [ "${proc}" == "-ggttg" ]; then 
      dir=$topdir/epochX/${bckend}/gg_ttg${suff}SubProcesses/P1_gg_ttxg
    elif [ "${proc}" == "-ggttgg" ]; then 
      dir=$topdir/epochX/${bckend}/gg_ttgg${suff}SubProcesses/P1_gg_ttxgg
    elif [ "${proc}" == "-ggttggg" ]; then 
      dir=$topdir/epochX/${bckend}/gg_ttggg${suff}SubProcesses/P1_gg_ttxggg
    elif [ "${proc}" == "-ggttgggg" ]; then 
      dir=$topdir/epochX/${bckend}/gg_ttgggg${suff}SubProcesses/P1_gg_ttxgggg
    elif [ "${proc}" == "-gqttq" ]; then 
      ###dir=$topdir/epochX/${bckend}/gq_ttq${suff}SubProcesses/P1_gu_ttxu
      dir=$topdir/epochX/${bckend}/gq_ttq${suff}SubProcesses/P1_gux_ttxux # only 1 out of 2 for now
    elif [ "${proc}" == "-heftggh" ]; then 
      echo "ERROR! Options -mad and -madonly are not supported with -heftggh"; exit 1
    fi
  else
    if [ "${proc}" == "-eemumu" ]; then 
      dir=$topdir/epochX/${bckend}/ee_mumu${suff}SubProcesses/P1_Sigma_sm_epem_mupmum
    elif [ "${proc}" == "-ggtt" ]; then 
      dir=$topdir/epochX/${bckend}/gg_tt${suff}SubProcesses/P1_Sigma_sm_gg_ttx
    elif [ "${proc}" == "-ggttg" ]; then 
      dir=$topdir/epochX/${bckend}/gg_ttg${suff}SubProcesses/P1_Sigma_sm_gg_ttxg
    elif [ "${proc}" == "-ggttgg" ]; then 
      dir=$topdir/epochX/${bckend}/gg_ttgg${suff}SubProcesses/P1_Sigma_sm_gg_ttxgg
    elif [ "${proc}" == "-ggttggg" ]; then 
      dir=$topdir/epochX/${bckend}/gg_ttggg${suff}SubProcesses/P1_Sigma_sm_gg_ttxggg
    elif [ "${proc}" == "-ggttgggg" ]; then 
      dir=$topdir/epochX/${bckend}/gg_ttgggg${suff}SubProcesses/P1_Sigma_sm_gg_ttxgggg
    elif [ "${proc}" == "-gqttq" ]; then 
      ###dir=$topdir/epochX/${bckend}/gq_ttq${suff}SubProcesses/P1_Sigma_sm_gu_ttxu
      dir=$topdir/epochX/${bckend}/gq_ttq${suff}SubProcesses/P1_Sigma_sm_gux_ttxux # only 1 out of 2 for now
    elif [ "${proc}" == "-heftggh" ]; then 
      dir=$topdir/epochX/${bckend}/heft_gg_h${suff}/SubProcesses/P1_Sigma_heft_gg_h
    fi
  fi
  echo $dir
}

###echo -e "\n********************************************************************************\n"
printf "\n"

##########################################################################
# PART 1 - compile the list of the executables which should be run
##########################################################################

dirs=
for proc in $procs; do
  for suff in $suffs; do
    dirs+=${dirs:+ }$(showdir)
  done
done
if [ "$dirs" == "" ]; then echo "ERROR! no valid directories found?"; exit 1; fi  

exes=

for dir in $dirs; do
  
  #=====================================
  # CUDA   (epochX - manual/mad)
  # ALPAKA (epochX - manual/auto)
  #=====================================
  if [ "${cuda}" == "1" ]; then
    for hrdcod in $hrdcods; do
      hrdsuf=_hrd${hrdcod}
      if [ "$bckend" == "alpaka" ]; then hrdsuf=""; fi
      for helinl in $helinls; do
        for fptype in $fptypes; do
          exes="$exes ${dir}/build.none_${fptype}_inl${helinl}${hrdsuf}/gcheck.exe"
          if [ "${alpaka}" == "1" ]; then
            exes="$exes ${dir}/build.none_${fptype}_inl${helinl}${hrdsuf}/alpcheck.exe"
          fi
        done
      done
    done
  fi
  
  #=====================================
  # C++ (eemumu/epochX - manual/mad)
  #=====================================
  if [ "${cpp}" == "1" ]; then 
    for hrdcod in $hrdcods; do
      hrdsuf=_hrd${hrdcod}
      if [ "$bckend" == "alpaka" ]; then hrdsuf=""; fi
      for helinl in $helinls; do
        for fptype in $fptypes; do
          for simd in none sse4 avx2 512y 512z; do
            if [ "${simds}" == "${simdsall}" ] || [ "${simds/${simd}}" != "${simds}" ]; then 
              exes="$exes $dir/build.${simd}_${fptype}_inl${helinl}${hrdsuf}/check.exe"
            fi
          done
        done
      done
    done
  fi

done

##########################################################################
# PART 2 - build the executables which should be run
##########################################################################

unset GTEST_ROOT
unset LOCALGTEST

if [ "${maketype}" == "-dryrun" ]; then

  printf "DRYRUN: SKIP MAKE\n\n"

else

  # Iterate over all directories (the first one will build googletest)
  gtestlibs=0
  for dir in $dirs; do
    export USEBUILDDIR=1
    pushd $dir >& /dev/null
    echo "Building in $(pwd)"
    if [ "${maketype}" != "-makecleanonly" ] && [ "${gtestlibs}" == "0" ]; then
      # Build googletest once and for all to avoid issues in parallel builds
      # NB1: $topdir/test is NO LONGER RELEVANT and googletest must be built from one specific process
      # NB2: CXXNAMESUFFIX must be set by cudacpp.mk, so googletest must be built from one P1 directory
      gtestlibs=1
      make -f cudacpp.mk gtestlibs
    fi
    if [ "${maketype}" == "-makeclean" ]; then make cleanall; echo; fi
    if [ "${maketype}" == "-makecleanonly" ]; then make cleanall; echo; continue; fi
    for hrdcod in $hrdcods; do
      export HRDCOD=$hrdcod
      for helinl in $helinls; do
	export HELINL=$helinl
	for fptype in $fptypes; do
          export FPTYPE=$fptype
          if [ "${simds}" == "${simdsall}" ]; then
            make ${makef} ${makej} avxall; echo # allow 'make avxall' again #536
          else
            if [ "${cuda}" == "1" ] || [ "${simds/none}" != "${simds}" ]; then 
              make ${makef} ${makej} AVX=none; echo
            fi
            for simd in ${simds/none}; do
              make ${makef} ${makej} AVX=${simd}; echo
            done
          fi
	done
      done
    done
    popd >& /dev/null
    export USEBUILDDIR=
    export HRDCOD=
    export HELINL=
    export FPTYPE=
  done

  if [ "${maketype}" == "-makecleanonly" ]; then printf "MAKE CLEANALL COMPLETED\n"; exit 0; fi
  if [ "${maketype}" == "-makeonly" ]; then printf "MAKE COMPLETED\n"; exit 0; fi

fi
  
##########################################################################
# PART 3 - run all the executables which should be run
##########################################################################

if [ "${maketype}" != "-dryrun" ]; then
  printf "DATE: $(date '+%Y-%m-%d_%H:%M:%S')\n\n"
fi

function runExe() {
  exe=$1
  args="$2"
  args="$args$rndgen$rmbsmp"
  echo "runExe $exe $args OMP=$OMP_NUM_THREADS"
  if [ "${maketype}" == "-dryrun" ]; then return; fi
  pattern="Process|fptype_sv|OMP threads|EvtsPerSec\[MECalc|MeanMatrix|FP precision|TOTAL       :"
  # Optionally add other patterns here for some specific configurations (e.g. clang)
  if [ "${exe%%/gcheck*}" != "${exe}" ]; then pattern="${pattern}|EvtsPerSec\[Matrix"; fi
  pattern="${pattern}|Workflow"
  ###pattern="${pattern}|CUCOMPLEX"
  ###pattern="${pattern}|COMMON RANDOM|CURAND HOST \(CUDA"
  pattern="${pattern}|ERROR"
  pattern="${pattern}|WARNING"
  pattern="${pattern}|EvtsPerSec\[Rmb" # TEMPORARY! for rambo timing tests
  pattern="${pattern}|EvtsPerSec\[Matrix" # TEMPORARY! OLD C++/CUDA CODE
  if [ "${ab3}" == "1" ]; then pattern="${pattern}|3a|3b"; fi
  if [ "${req}" == "1" ]; then pattern="${pattern}|memory layout"; fi
  if perf stat -d date >& /dev/null; then # perf exists and CAP_SYS_ADMIN allows users to use it
    # -- Newer version using perf stat
    pattern="${pattern}|instructions|cycles"
    pattern="${pattern}|elapsed"
    if [ "${detailed}" == "1" ]; then pattern="${pattern}|#"; fi
    if [ "${verbose}" == "1" ]; then set -x; fi
    ###perf stat -d $exe $args 2>&1 | grep -v "Performance counter stats"
    perf stat -d $exe $args 2>&1 | egrep "(${pattern})" | grep -v "Performance counter stats"
    set +x
  else
    # -- Older version using time
    # For TIMEFORMAT see https://www.gnu.org/software/bash/manual/html_node/Bash-Variables.html
    if [ "${verbose}" == "1" ]; then set -x; fi
    TIMEFORMAT=$'real\t%3lR' && time $exe $args 2>&1 | egrep "(${pattern})"
    set +x
  fi
}

function cmpExe() {
  exe=$1
  exef=${exe/\/check//fcheck}
  exef=${exef/\/gcheck//fgcheck}
  argsf="2 64 2"
  args="--common -p ${argsf}"
  echo "cmpExe $exe $args"
  echo "cmpExe $exef $argsf"
  if [ "${maketype}" == "-dryrun" ]; then return; fi
  tmp=$(mktemp)
  me1=$(${exe} ${args} 2>${tmp} | grep MeanMatrix | awk '{print $4}'); cat ${tmp}
  me2=$(${exef} ${argsf} 2>${tmp} | grep Average | awk '{print $4}'); cat ${tmp}
  if [ "${exe%%/gcheck*}" != "${exe}" ]; then tag="/CUDA)"; else tag="/C++) "; fi
  echo -e "Avg ME (C++${tag}   = ${me1}\nAvg ME (F77${tag}   = ${me2}"
  if [ "${me2}" == "NaN" ]; then
    echo "ERROR! Fortran calculation (F77${tag} returned NaN"
  elif [ "${me2}" == "" ]; then
    echo "ERROR! Fortran calculation (F77${tag} crashed"
  else
    # NB skip python comparison if Fortran returned NaN or crashed, otherwise python returns an error status and the following tests are not executed
    python3 -c "me1=${me1}; me2=${me2}; reldif=abs((me2-me1)/me1); print('Relative difference =', reldif); ok = reldif <= 5E-3; print ( '%s (relative difference %s 5E-3)' % ( ('OK','<=') if ok else ('ERROR','>') ) )"
  fi
}

# Profile #registers and %divergence only
function runNcu() {
  if [ "${maketype}" == "-dryrun" ]; then return; fi
  exe=$1
  args="$2"
  args="$args$rndgen$rmbsmp"
  echo "runNcu $exe $args"
  if [ "${verbose}" == "1" ]; then set -x; fi
  #$(which ncu) --metrics launch__registers_per_thread,sm__sass_average_branch_targets_threads_uniform.pct --target-processes all --kernel-id "::sigmaKin:" --kernel-name-base function $exe $args | egrep '(sigmaKin|registers| sm)' | tr "\n" " " | awk '{print $1, $2, $3, $15, $17; print $1, $2, $3, $18, $20$19}'
  set +e # do not fail on error
  out=$($(which ncu) --metrics launch__registers_per_thread,sm__sass_average_branch_targets_threads_uniform.pct --target-processes all --kernel-id "::sigmaKin:" --kernel-name-base function $exe $args)
  echo "$out" | egrep '(ERROR|WARNING)' # NB must escape $out in between quotes
  set -e # fail on error (after ncu and after egrep!)
  out=$(echo "${out}" | egrep '(sigmaKin|registers| sm)' | tr "\n" " ") # NB must escape $out in between quotes
  echo $out | awk -v key1="launch__registers_per_thread" '{val1="N/A"; for (i=1; i<=NF; i++){if ($i==key1 && $(i+1)!="(!)") val1=$(i+2)}; print $1, $2, $3, key1, val1}'
  echo $out | awk -v key1="sm__sass_average_branch_targets_threads_uniform.pct" '{val1="N/A"; for (i=1; i<=NF; i++){if ($i==key1 && $(i+1)!="(!)") val1=$(i+2)$(i+1)}; print $1, $2, $3, key1, val1}'
  set +x
}

# Profile divergence metrics more in detail
# See https://www.pgroup.com/resources/docs/18.10/pdf/pgi18profug.pdf
# See https://docs.nvidia.com/gameworks/content/developertools/desktop/analysis/report/cudaexperiments/kernellevel/branchstatistics.htm
# See https://docs.nvidia.com/gameworks/content/developertools/desktop/analysis/report/cudaexperiments/sourcelevel/divergentbranch.htm
function runNcuDiv() {
  if [ "${maketype}" == "-dryrun" ]; then return; fi
  exe=$1
  args="-p 1 32 1"
  args="$args$rndgen$rmbsmp"
  ###echo "runNcuDiv $exe $args"
  if [ "${verbose}" == "1" ]; then set -x; fi
  ###$(which ncu) --query-metrics $exe $args
  ###$(which ncu) --metrics regex:.*branch_targets.* --target-processes all --kernel-id "::sigmaKin:" --kernel-name-base function $exe $args
  ###$(which ncu) --metrics regex:.*stalled_barrier.* --target-processes all --kernel-id "::sigmaKin:" --kernel-name-base function $exe $args
  ###$(which ncu) --metrics sm__sass_average_branch_targets_threads_uniform.pct,smsp__warps_launched.sum,smsp__sass_branch_targets.sum,smsp__sass_branch_targets_threads_divergent.sum,smsp__sass_branch_targets_threads_uniform.sum --target-processes all --kernel-id "::sigmaKin:" --kernel-name-base function $exe $args | egrep '(sigmaKin| sm)' | tr "\n" " " | awk '{printf "%29s: %-51s %s\n", "", $18, $19; printf "%29s: %-51s %s\n", "", $22, $23; printf "%29s: %-51s %s\n", "", $20, $21; printf "%29s: %-51s %s\n", "", $24, $26}'
  #$(which ncu) --metrics sm__sass_average_branch_targets_threads_uniform.pct,smsp__warps_launched.sum,smsp__sass_branch_targets.sum,smsp__sass_branch_targets_threads_divergent.sum,smsp__sass_branch_targets_threads_uniform.sum,smsp__sass_branch_targets.sum.per_second,smsp__sass_branch_targets_threads_divergent.sum.per_second,smsp__sass_branch_targets_threads_uniform.sum.per_second --target-processes all --kernel-id "::sigmaKin:" --kernel-name-base function $exe $args | egrep '(sigmaKin| sm)' | tr "\n" " " | awk '{printf "%29s: %-51s %-10s %s\n", "", $18, $19, $22$21; printf "%29s: %-51s %-10s %s\n", "", $28, $29, $32$31; printf "%29s: %-51s %-10s %s\n", "", $23, $24, $27$26; printf "%29s: %-51s %s\n", "", $33, $35}'
  out=$($(which ncu) --metrics sm__sass_average_branch_targets_threads_uniform.pct,smsp__warps_launched.sum,smsp__sass_branch_targets.sum,smsp__sass_branch_targets_threads_divergent.sum,smsp__sass_branch_targets_threads_uniform.sum,smsp__sass_branch_targets.sum.per_second,smsp__sass_branch_targets_threads_divergent.sum.per_second,smsp__sass_branch_targets_threads_uniform.sum.per_second --target-processes all --kernel-id "::sigmaKin:" --kernel-name-base function $exe $args | egrep '(sigmaKin| sm)' | tr "\n" " ")
  ###echo $out
  echo $out | awk -v key1="smsp__sass_branch_targets.sum" '{key2=key1".per_second"; val1="N/A"; val2=""; for (i=1; i<=NF; i++){if ($i==key1 && $(i+1)!="(!)") val1=$(i+1); if ($i==key2 && $(i+1)!="(!)") val2=$(i+2)$(i+1)}; printf "%29s: %-51s %-10s %s\n", "", key1, val1, val2}'
  echo $out | awk -v key1="smsp__sass_branch_targets_threads_uniform.sum" '{key2=key1".per_second"; val1="N/A"; val2=""; for (i=1; i<=NF; i++){if ($i==key1 && $(i+1)!="(!)") val1=$(i+1); if ($i==key2 && $(i+1)!="(!)") val2=$(i+2)$(i+1)}; printf "%29s: %-51s %-10s %s\n", "", key1, val1, val2}'
  echo $out | awk -v key1="smsp__sass_branch_targets_threads_divergent.sum" '{key2=key1".per_second"; val1="N/A"; val2=""; for (i=1; i<=NF; i++){if ($i==key1 && $(i+1)!="(!)") val1=$(i+1); if ($i==key2 && $(i+1)!="(!)") val2=$(i+2)$(i+1)}; printf "%29s: %-51s %-10s %s\n", "", key1, val1, val2}'
  echo $out | awk -v key="smsp__warps_launched.sum" '{val1="N/A"; for (i=1; i<=NF; i++){if ($i==key && $(i+1)!="(!)") val1=$(i+2)}; printf "%29s: %-51s %s\n", "", key, val1}'
  set +x
}

# Profiles sectors and requests
function runNcuReq() {
  if [ "${maketype}" == "-dryrun" ]; then return; fi
  exe=$1
  ncuArgs="$2"
  ncuArgs="$ncuArgs$rndgen$rmbsmp"
  if [ "${verbose}" == "1" ]; then set -x; fi
  for args in "-p 1 1 1" "-p 1 4 1" "-p 1 8 1" "-p 1 32 1" "$ncuArgs"; do
    ###echo "runNcuReq $exe $args"
    # NB This will print nothing if $args are invalid (eg "-p 1 4 1" when neppR=8)
    $(which ncu) --metrics l1tex__t_sectors_pipe_lsu_mem_global_op_ld.sum,l1tex__t_requests_pipe_lsu_mem_global_op_ld.sum,launch__registers_per_thread,sm__sass_average_branch_targets_threads_uniform.pct --target-processes all --kernel-id "::sigmaKin:" --kernel-name-base function $exe $args | egrep '(sigmaKin|registers| sm|l1tex)' | tr "\n" " " | awk -vtag="[$args]" '{print $1, $2, $3, $16"s", $17";", $19"s", $20, tag}'
  done
  set +x
}

if nvidia-smi -L > /dev/null 2>&1; then gpuTxt="$(nvidia-smi -L | wc -l)x $(nvidia-smi -L | awk '{print $3,$4}' | sort -u)"; else gpuTxt=none; fi
if [ "${unames}" == "Darwin" ]; then 
  cpuTxt=$(sysctl -h machdep.cpu.brand_string)
  cpuTxt=${cpuTxt/machdep.cpu.brand_string: }
elif [ "${unamep}" == "ppc64le" ]; then 
  cpuTxt=$(cat /proc/cpuinfo | grep ^machine | awk '{print substr($0,index($0,"Power"))", "}')$(cat /proc/cpuinfo | grep ^cpu | head -1 | awk '{print substr($0,index($0,"POWER"))}')
else
  cpuTxt=$(cat /proc/cpuinfo | grep '^model name' |& head -1 | awk '{i0=index($0,"Intel"); if (i0==0) i0=index($0,"AMD"); i1=index($0," @"); if (i1>0) {print substr($0,i0,i1-i0)} else {print substr($0,i0)}}')
fi
echo -e "On $HOSTNAME [CPU: $cpuTxt] [GPU: $gpuTxt]:"

# These two settings are needed by BMK containers: do not change them
BMKEXEARGS="" # if BMKEXEARGS is set, exeArgs is set equal to BMKEXEARGS, while exeArgs2 is set to ""
BMKMULTIPLIER=1 # the pre-defined numbers of iterations (including those in BMKEXEARGS) are multiplied by BMKMULTIPLIER

lastExe=
###echo "exes=$exes"
for exe in $exes; do
  ###echo EXE=$exe; continue
  exeArgs2=""
  if [ "$(basename $exe)" != "$lastExe" ]; then
    echo "========================================================================="
    lastExe=$(basename $exe)
  else
    echo "-------------------------------------------------------------------------"
  fi
  if [ ! -f $exe ]; then echo "Not found: $exe"; continue; fi
  if [ "${unamep}" != "x86_64" ]; then
    if [ "${exe/build.avx2}" != "${exe}" ]; then echo "$exe is not supported on ${unamep}"; continue; fi
    if [ "${exe/build.512y}" != "${exe}" ]; then echo "$exe is not supported on ${unamep}"; continue; fi
    if [ "${exe/build.512z}" != "${exe}" ]; then echo "$exe is not supported on ${unamep}"; continue; fi
  elif [ "${unames}" == "Darwin" ]; then
    if [ "${exe/build.512y}" != "${exe}" ]; then echo "$exe is not supported on ${unames}"; continue; fi
    if [ "${exe/build.512z}" != "${exe}" ]; then echo "$exe is not supported on ${unames}"; continue; fi
  elif [ "$(grep -m1 -c avx512vl /proc/cpuinfo)" != "1" ]; then
    if [ "${exe/build.512y}" != "${exe}" ]; then echo "$exe is not supported (no avx512vl in /proc/cpuinfo)"; continue; fi
    if [ "${exe/build.512z}" != "${exe}" ]; then echo "$exe is not supported (no avx512vl in /proc/cpuinfo)"; continue; fi
  fi
  if [ "${exe%%/gcheck*}" != "${exe}" ] && [ "$gpuTxt" == "none" ]; then continue; fi
  if [ "${exe%%/heft_gg_h*}" != "${exe}" ]; then 
    # For heftggh: 2->1 process, hence all events are identical and random numbers are ignored, use bare minimum 1 8 1
    exeArgs="-p 1 8 1"
    ncuArgs="-p 1 8 1"
  elif [ "${exe%%/gg_ttgggg*}" != "${exe}" ]; then 
    # For ggttgggg: initially try the same settings as ggttggg
    exeArgs="-p 1 256 2"
    ncuArgs="-p 1 256 1"
    # For ggttgggg: initially try the same settings as ggttggg
    exeArgs2="-p 64 256 1"
  elif [ "${exe%%/gq_ttq*}" != "${exe}" ]; then 
    # For gqttq, use the same settings as for ggttg
    exeArgs="-p 64 256 10"
    ncuArgs="-p 64 256 1"
    # For gqttq, use the same settings as for ggttg
    exeArgs2="-p 2048 256 1"
  elif [ "${exe%%/gg_ttggg*}" != "${exe}" ]; then 
    # For ggttggg: this is far too little for GPU (4.8E2), but it keeps the CPU to a manageble level (1sec with 512y)
    ###exeArgs="-p 1 256 1" # too short! see https://its.cern.ch/jira/browse/BMK-1056
    exeArgs="-p 1 256 2"
    ncuArgs="-p 1 256 1"
    # For ggttggg: on GPU test also "64 256" to reach the plateau (only ~5% lower than "2048 256": 1.18E4 vs 1.25E4 on cuda116/gcc102)
    exeArgs2="-p 64 256 1"
  elif [ "${exe%%/gg_ttgg*}" != "${exe}" ]; then 
    # For ggttgg (OLD): this is a good GPU middle point: tput is 1.5x lower with "32 256 1", only a few% higher with "128 256 1"
    exeArgs="-p 64 256 1"
    ncuArgs="-p 64 256 1"
    # For ggttgg (NEW): on GPU test both "64 256" and "2048 256" for ggttgg as the latter gives ~10% higher throughput on cuda110/gcc92
    exeArgs2="-p 2048 256 1"
  elif [ "${exe%%/gg_ttg*}" != "${exe}" ]; then 
    # For ggttg, as on ggttgg: this is a good GPU middle point: tput is 1.5x lower with "32 256 1", only a few% higher with "128 256 1"
    ###exeArgs="-p 64 256 1" # too short! see https://its.cern.ch/jira/browse/BMK-1056
    exeArgs="-p 64 256 10"
    ncuArgs="-p 64 256 1"
    # For ggttg, as on ggttgg: on GPU test both "64 256" and "2048 256" for ggttg as the latter gives ~10% higher throughput on cuda110/gcc92
    exeArgs2="-p 2048 256 1"
  elif [ "${exe%%/gg_tt*}" != "${exe}" ]; then 
    ###exeArgs="-p 2048 256 1" # too short! see https://its.cern.ch/jira/browse/BMK-1056
    exeArgs="-p 2048 256 2"
    ncuArgs="-p 2048 256 1"
  else # eemumu
    exeArgs="-p 2048 256 12"
    ncuArgs="-p 2048 256 1"
  fi
  if [ "$BMKEXEARGS" ]; then exeArgs="$BMKEXEARGS"; exeArgs2=""; fi
  exeArgs="${exeArgs% *} $((${exeArgs##* }*BMKMULTIPLIER))"
  if [ "${exeArgs2}" != "" ]; then exeArgs2="${exeArgs2% *} $((${exeArgs2##* }*BMKMULTIPLIER))"; fi  
  exeDir=$(dirname $exe)
  cd $exeDir/.. # workaround for reading '../../Cards/param_card.dat' without setting MG5AMC_CARD_PATH
  unset OMP_NUM_THREADS
  runExe $exe "$exeArgs"
  if [ "${exe%%/check*}" != "${exe}" ]; then 
    if [ "${maketype}" != "-dryrun" ]; then
      obj=${exe%%/check*}/CPPProcess.o; $scrdir/simdSymSummary.sh -stripdir ${obj} -dumptotmp # comment out -dumptotmp to keep full objdump
    fi
    if [ "${omp}" == "1" ]; then 
      echo "-------------------------------------------------------------------------"
      export OMP_NUM_THREADS=$(nproc --all)
      runExe $exe "$exeArgs"
      unset OMP_NUM_THREADS
    fi
    if [ "${gtest}" == "1" ]; then
      echo "-------------------------------------------------------------------------"
      exe2=${exe/check/runTest}
      echo "runExe $exe2"
      $exe2 2>&1 | tail -1
      if [ ${PIPESTATUS[0]} -ne "0" ]; then exit 1; fi 
    fi
  elif [ "${exe%%/gcheck*}" != "${exe}" ] ||  [ "${exe%%/alpcheck*}" != "${exe}" ]; then 
    runNcu $exe "$ncuArgs"
    if [ "${div}" == "1" ]; then runNcuDiv $exe; fi
    if [ "${req}" == "1" ]; then runNcuReq $exe "$ncuArgs"; fi
    if [ "${exeArgs2}" != "" ]; then echo "........................................................................."; runExe $exe "$exeArgs2"; fi
  fi
  if [ "${bckend}" != "alpaka" ]; then
    echo "-------------------------------------------------------------------------"
    cmpExe $exe
  fi
done
echo "========================================================================="

# Workaround for reading of data files
printf "\nTEST COMPLETED\n"
