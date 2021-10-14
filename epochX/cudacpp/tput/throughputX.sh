#!/bin/bash

set +x # not verbose
set -e # fail on error

omp=0
avxall=0
cpp=1
cuda=1
ab3=0
eemumu=0
ggtt=0
ggttgg=0
div=0
req=0
fptypes="d"
helinls="0"
suffs="/"
makeonly=0
detailed=0
verbose=0

function usage()
{
  echo "Usage: $0 [-nocpp|[-omp][-avxall][-nocuda]] [-3a3b] [-eemumu] [-ggtt] [-ggttgg] [-div] [-req] [-flt|-fltonly] [-inl|-inlonly] [-auto|-autoonly] [-makeonly] [-detailed] [-v]"
  exit 1
}

##########################################################################
# PART 0 - decode command line arguments
##########################################################################

scrdir=$(cd $(dirname $0); pwd)
topdir=$(cd $scrdir; cd ../../..; pwd)

while [ "$1" != "" ]; do
  if [ "$1" == "-omp" ]; then
    if [ "${cpp}" == "0" ]; then echo "ERROR! Options -omp and -nocpp are incompatible"; usage; fi
    omp=1
    shift
  elif [ "$1" == "-avxall" ]; then
    if [ "${cpp}" == "0" ]; then echo "ERROR! Options -avxall and -nocpp are incompatible"; usage; fi
    avxall=1
    shift
  elif [ "$1" == "-nocuda" ]; then
    if [ "${cpp}" == "0" ]; then echo "ERROR! Options -nocuda and -nocpp are incompatible"; usage; fi
    cuda=0
    shift
  elif [ "$1" == "-nocpp" ]; then
    if [ "${omp}" == "1" ]; then echo "ERROR! Options -omp and -nocpp are incompatible"; usage; fi
    if [ "${avxall}" == "1" ]; then echo "ERROR! Options -avxall and -nocpp are incompatible"; usage; fi
    if [ "${cuda}" == "1" ]; then echo "ERROR! Options -nocuda and -nocpp are incompatible"; usage; fi
    cpp=0
    shift
  elif [ "$1" == "-3a3b" ]; then
    ab3=1
    shift
  elif [ "$1" == "-eemumu" ]; then
    eemumu=1
    shift
  elif [ "$1" == "-ggtt" ]; then
    ggtt=1
    shift
  elif [ "$1" == "-ggttgg" ]; then
    ggttgg=1
    shift
  elif [ "$1" == "-div" ]; then
    div=1
    shift
  elif [ "$1" == "-req" ]; then
    req=1
    shift
  elif [ "$1" == "-flt" ]; then
    if [ "${fptypes}" == "f" ]; then echo "ERROR! Options -flt and -fltonly are incompatible"; usage; fi
    fptypes="d f"
    shift
  elif [ "$1" == "-fltonly" ]; then
    if [ "${fptypes}" == "d f" ]; then echo "ERROR! Options -flt and -fltonly are incompatible"; usage; fi
    fptypes="f"
    shift
  elif [ "$1" == "-inl" ]; then
    if [ "${helinls}" == "1" ]; then echo "ERROR! Options -inl and -inlonly are incompatible"; usage; fi
    helinls="0 1"
    shift
  elif [ "$1" == "-inlonly" ]; then
    if [ "${helinls}" == "0 1" ]; then echo "ERROR! Options -inl and -inlonly are incompatible"; usage; fi
    helinls="1"
    shift
  elif [ "$1" == "-auto" ]; then
    if [ "${suffs}" == ".auto/" ]; then echo "ERROR! Options -auto and -autoonly are incompatible"; usage; fi
    suffs="/ .auto/"
    shift
  elif [ "$1" == "-autoonly" ]; then
    if [ "${suffs}" == "/ .auto/" ]; then echo "ERROR! Options -auto and -autoonly are incompatible"; usage; fi
    suffs=".auto/"
    shift
  elif [ "$1" == "-makeonly" ]; then
    makeonly=1
    shift
  elif [ "$1" == "-detailed" ]; then
    detailed=1
    shift
  elif [ "$1" == "-v" ]; then
    verbose=1
    shift
  else
    usage
  fi
done
###exit 1

# New default: run both eemumu and ggttgg
###if [ "${eemumu}" == "0" ] && [ "${ggttgg}" == "0" ]; then 
###  eemumu=1
###  ggttgg=1
###fi

###echo -e "\n********************************************************************************\n"
printf "\n"

##########################################################################
# PART 1 - compile the list of the executables which should be run
##########################################################################

exes=

for suff in $suffs; do

  #=====================================
  # CUDA (NEW eemumu/epochX - manual)
  #=====================================
  if [ "${cuda}" == "1" ]; then
    if [ "${eemumu}" == "1" ] && [ "${suff}" != ".auto/" ]; then 
      for helinl in $helinls; do
        for fptype in $fptypes; do
          exes="$exes $topdir/epochX/cudacpp/ee_mumu${suff}SubProcesses/P1_Sigma_sm_epem_mupmum/build.none_${fptype}_inl${helinl}/gcheck.exe"
        done
      done
    fi
  fi
  
  #=====================================
  # C++ (NEW eemumu/epochX - manual)
  #=====================================
  if [ "${cpp}" == "1" ]; then 
    if [ "${eemumu}" == "1" ] && [ "${suff}" != ".auto/" ]; then 
      for helinl in $helinls; do
        for fptype in $fptypes; do
          exes="$exes $topdir/epochX/cudacpp/ee_mumu${suff}SubProcesses/P1_Sigma_sm_epem_mupmum/build.none_${fptype}_inl${helinl}/check.exe"
          if [ "${avxall}" == "1" ]; then 
            exes="$exes $topdir/epochX/cudacpp/ee_mumu${suff}SubProcesses/P1_Sigma_sm_epem_mupmum/build.sse4_${fptype}_inl${helinl}/check.exe"
            exes="$exes $topdir/epochX/cudacpp/ee_mumu${suff}SubProcesses/P1_Sigma_sm_epem_mupmum/build.avx2_${fptype}_inl${helinl}/check.exe"
          fi
          if [ "$(grep -m1 -c avx512vl /proc/cpuinfo)" == "1" ]; then 
            exes="$exes $topdir/epochX/cudacpp/ee_mumu${suff}SubProcesses/P1_Sigma_sm_epem_mupmum/build.512y_${fptype}_inl${helinl}/check.exe"
            if [ "${avxall}" == "1" ]; then 
              exes="$exes $topdir/epochX/cudacpp/ee_mumu${suff}SubProcesses/P1_Sigma_sm_epem_mupmum/build.512z_${fptype}_inl${helinl}/check.exe"
            fi
          fi
        done
      done
    fi
  fi

  #=====================================
  # CUDA (OLD eemumu/epochX - auto)
  #=====================================
  if [ "${cuda}" == "1" ]; then
    if [ "${eemumu}" == "1" ] && [ "${suff}" == ".auto/" ]; then 
      exes="$exes $topdir/epochX/cudacpp/ee_mumu${suff}SubProcesses/P1_Sigma_sm_epem_mupmum/gcheck.exe"
    fi
  fi

  #=====================================
  # C++ (OLD eemumu/epochX - auto)
  #=====================================
  if [ "${cpp}" == "1" ]; then 
    if [ "${eemumu}" == "1" ] && [ "${suff}" == ".auto/" ]; then 
      exes="$exes $topdir/epochX/cudacpp/ee_mumu${suff}SubProcesses/P1_Sigma_sm_epem_mupmum/check.exe"
    fi
  fi

  #=====================================
  # CUDA (ggtt/epochX)
  #=====================================
  if [ "${cuda}" == "1" ]; then
    if [ "${ggtt}" == "1" ]; then 
      exes="$exes $topdir/epochX/cudacpp/gg_tt${suff}SubProcesses/P1_Sigma_sm_gg_ttx/gcheck.exe"
    fi
  fi

  #=====================================
  # C++ (ggtt/epochX)
  #=====================================
  if [ "${cpp}" == "1" ]; then 
    if [ "${ggtt}" == "1" ]; then 
      exes="$exes $topdir/epochX/cudacpp/gg_tt${suff}SubProcesses/P1_Sigma_sm_gg_ttx/check.exe"
    fi
  fi

  #=====================================
  # CUDA (ggttgg/epochX)
  #=====================================
  if [ "${cuda}" == "1" ]; then
    if [ "${ggttgg}" == "1" ]; then 
      exes="$exes $topdir/epochX/cudacpp/gg_ttgg${suff}SubProcesses/P1_Sigma_sm_gg_ttxgg/gcheck.exe"
    fi
  fi

  #=====================================
  # C++ (ggttgg/epochX)
  #=====================================
  if [ "${cpp}" == "1" ]; then 
    if [ "${ggttgg}" == "1" ]; then 
      exes="$exes $topdir/epochX/cudacpp/gg_ttgg${suff}SubProcesses/P1_Sigma_sm_gg_ttxgg/check.exe"
    fi
  fi

done

##########################################################################
# PART 2 - build the executables which should be run
##########################################################################

###echo "exes=$exes"

for suff in $suffs; do

  if [ "${eemumu}" == "1" ] && [ "${suff}" != ".auto/" ]; then 
    export USEBUILDDIR=1
    pushd $topdir/epochX/cudacpp/ee_mumu${suff}SubProcesses/P1_Sigma_sm_epem_mupmum >& /dev/null
    pwd
    for helinl in $helinls; do
      export HELINL=$helinl
      for fptype in $fptypes; do
        export FPTYPE=$fptype
        make AVX=none; echo
        if [ "${avxall}" == "1" ]; then make AVX=sse4; echo; fi
        if [ "${avxall}" == "1" ]; then make AVX=avx2; echo; fi
        if [ "$(grep -m1 -c avx512vl /proc/cpuinfo)" == "1" ]; then # skip avx512 if not supported!
          if [ "${cpp}" == "1" ]; then make AVX=512y; echo; fi # use 512y as C++ ref even if avx2 is faster on clang
          if [ "${avxall}" == "1" ]; then make AVX=512z; echo; fi
        fi
      done
    done
    popd >& /dev/null
    export USEBUILDDIR=
    export HELINL=
    export FPTYPE=
  fi

  if [ "${eemumu}" == "1" ] && [ "${suff}" == ".auto/" ]; then 
    pushd $topdir/epochX/cudacpp/ee_mumu${suff}SubProcesses/P1_Sigma_sm_epem_mupmum >& /dev/null
    pwd
    ###make; echo
    make AVX=none; echo
    popd >& /dev/null
  fi

  if [ "${ggtt}" == "1" ]; then 
    pushd $topdir/epochX/cudacpp/gg_tt${suff}SubProcesses/P1_Sigma_sm_gg_ttx >& /dev/null
    pwd
    make; echo
    popd >& /dev/null
  fi

  if [ "${ggttgg}" == "1" ]; then 
    pushd $topdir/epochX/cudacpp/gg_ttgg${suff}SubProcesses/P1_Sigma_sm_gg_ttxgg >& /dev/null
    pwd
    make; echo
    popd >& /dev/null
  fi

done

if [ "${makeonly}" == "1" ]; then printf "BUILD COMPLETED\n"; exit 0; fi

##########################################################################
# PART 3 - run all the executables which should be run
##########################################################################

printf "DATE: $(date '+%Y-%m-%d_%H:%M:%S')\n\n"

function runExe() {
  exe=$1
  args="$2"
  echo "runExe $exe $args OMP=$OMP_NUM_THREADS"
  pattern="Process|fptype_sv|OMP threads|EvtsPerSec\[MECalc|MeanMatrix|FP precision|TOTAL       :"
  # Optionally add other patterns here for some specific configurations (e.g. clang)
  if [ "${exe%%/gcheck*}" != "${exe}" ]; then pattern="${pattern}|EvtsPerSec\[Matrix"; fi
  pattern="${pattern}|CUCOMPLEX"
  pattern="${pattern}|COMMON RANDOM"
  pattern="${pattern}|ERROR"
  # TEMPORARY! OLD C++/CUDA CODE (START)
  pattern="${pattern}|EvtsPerSec\[Matrix"
  # TEMPORARY! OLD C++/CUDA CODE (END)
  if [ "${ab3}" == "1" ]; then pattern="${pattern}|3a|3b"; fi
  if [ "${req}" == "1" ]; then pattern="${pattern}|memory layout"; fi
  if perf --version >& /dev/null; then
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

# Profile #registers and %divergence only
function runNcu() {
  exe=$1
  args="$2"
  ###echo "runNcu $exe $args"
  if [ "${verbose}" == "1" ]; then set -x; fi
  #$(which ncu) --metrics launch__registers_per_thread,sm__sass_average_branch_targets_threads_uniform.pct --target-processes all --kernel-id "::sigmaKin:" --kernel-base mangled $exe $args | egrep '(sigmaKin|registers| sm)' | tr "\n" " " | awk '{print $1, $2, $3, $15, $17; print $1, $2, $3, $18, $20$19}'
  out=$($(which ncu) --metrics launch__registers_per_thread,sm__sass_average_branch_targets_threads_uniform.pct --target-processes all --kernel-id "::sigmaKin:" --kernel-base mangled $exe $args | egrep '(sigmaKin|registers| sm)' | tr "\n" " ")
  ###echo $out
  echo $out | awk -v key1="launch__registers_per_thread" '{val1="N/A"; for (i=1; i<=NF; i++){if ($i==key1 && $(i+1)!="(!)") val1=$(i+2)}; print $1, $2, $3, key1, val1}'
  echo $out | awk -v key1="sm__sass_average_branch_targets_threads_uniform.pct" '{val1="N/A"; for (i=1; i<=NF; i++){if ($i==key1 && $(i+1)!="(!)") val1=$(i+2)$(i+1)}; print $1, $2, $3, key1, val1}'
  set +x
}

# Profile divergence metrics more in detail
# See https://www.pgroup.com/resources/docs/18.10/pdf/pgi18profug.pdf
# See https://docs.nvidia.com/gameworks/content/developertools/desktop/analysis/report/cudaexperiments/kernellevel/branchstatistics.htm
# See https://docs.nvidia.com/gameworks/content/developertools/desktop/analysis/report/cudaexperiments/sourcelevel/divergentbranch.htm
function runNcuDiv() {
  exe=$1
  args="-p 1 32 1"
  ###echo "runNcuDiv $exe $args"
  if [ "${verbose}" == "1" ]; then set -x; fi
  ###$(which ncu) --query-metrics $exe $args
  ###$(which ncu) --metrics regex:.*branch_targets.* --target-processes all --kernel-id "::sigmaKin:" --kernel-base mangled $exe $args
  ###$(which ncu) --metrics regex:.*stalled_barrier.* --target-processes all --kernel-id "::sigmaKin:" --kernel-base mangled $exe $args
  ###$(which ncu) --metrics sm__sass_average_branch_targets_threads_uniform.pct,smsp__warps_launched.sum,smsp__sass_branch_targets.sum,smsp__sass_branch_targets_threads_divergent.sum,smsp__sass_branch_targets_threads_uniform.sum --target-processes all --kernel-id "::sigmaKin:" --kernel-base mangled $exe $args | egrep '(sigmaKin| sm)' | tr "\n" " " | awk '{printf "%29s: %-51s %s\n", "", $18, $19; printf "%29s: %-51s %s\n", "", $22, $23; printf "%29s: %-51s %s\n", "", $20, $21; printf "%29s: %-51s %s\n", "", $24, $26}'
  #$(which ncu) --metrics sm__sass_average_branch_targets_threads_uniform.pct,smsp__warps_launched.sum,smsp__sass_branch_targets.sum,smsp__sass_branch_targets_threads_divergent.sum,smsp__sass_branch_targets_threads_uniform.sum,smsp__sass_branch_targets.sum.per_second,smsp__sass_branch_targets_threads_divergent.sum.per_second,smsp__sass_branch_targets_threads_uniform.sum.per_second --target-processes all --kernel-id "::sigmaKin:" --kernel-base mangled $exe $args | egrep '(sigmaKin| sm)' | tr "\n" " " | awk '{printf "%29s: %-51s %-10s %s\n", "", $18, $19, $22$21; printf "%29s: %-51s %-10s %s\n", "", $28, $29, $32$31; printf "%29s: %-51s %-10s %s\n", "", $23, $24, $27$26; printf "%29s: %-51s %s\n", "", $33, $35}'
  out=$($(which ncu) --metrics sm__sass_average_branch_targets_threads_uniform.pct,smsp__warps_launched.sum,smsp__sass_branch_targets.sum,smsp__sass_branch_targets_threads_divergent.sum,smsp__sass_branch_targets_threads_uniform.sum,smsp__sass_branch_targets.sum.per_second,smsp__sass_branch_targets_threads_divergent.sum.per_second,smsp__sass_branch_targets_threads_uniform.sum.per_second --target-processes all --kernel-id "::sigmaKin:" --kernel-base mangled $exe $args | egrep '(sigmaKin| sm)' | tr "\n" " ")
  ###echo $out
  echo $out | awk -v key1="smsp__sass_branch_targets.sum" '{key2=key1".per_second"; val1="N/A"; val2=""; for (i=1; i<=NF; i++){if ($i==key1 && $(i+1)!="(!)") val1=$(i+1); if ($i==key2 && $(i+1)!="(!)") val2=$(i+2)$(i+1)}; printf "%29s: %-51s %-10s %s\n", "", key1, val1, val2}'
  echo $out | awk -v key1="smsp__sass_branch_targets_threads_uniform.sum" '{key2=key1".per_second"; val1="N/A"; val2=""; for (i=1; i<=NF; i++){if ($i==key1 && $(i+1)!="(!)") val1=$(i+1); if ($i==key2 && $(i+1)!="(!)") val2=$(i+2)$(i+1)}; printf "%29s: %-51s %-10s %s\n", "", key1, val1, val2}'
  echo $out | awk -v key1="smsp__sass_branch_targets_threads_divergent.sum" '{key2=key1".per_second"; val1="N/A"; val2=""; for (i=1; i<=NF; i++){if ($i==key1 && $(i+1)!="(!)") val1=$(i+1); if ($i==key2 && $(i+1)!="(!)") val2=$(i+2)$(i+1)}; printf "%29s: %-51s %-10s %s\n", "", key1, val1, val2}'
  echo $out | awk -v key="smsp__warps_launched.sum" '{val1="N/A"; for (i=1; i<=NF; i++){if ($i==key && $(i+1)!="(!)") val1=$(i+2)}; printf "%29s: %-51s %s\n", "", key, val1}'
  set +x
}

# Profiles sectors and requests
function runNcuReq() {
  exe=$1
  ncuArgs="$2"
  if [ "${verbose}" == "1" ]; then set -x; fi
  for args in "-p 1 1 1" "-p 1 4 1" "-p 1 8 1" "-p 1 32 1" "$ncuArgs"; do
    ###echo "runNcuReq $exe $args"
    # NB This will print nothing if $args are invalid (eg "-p 1 4 1" when neppR=8)
    $(which ncu) --metrics l1tex__t_sectors_pipe_lsu_mem_global_op_ld.sum,l1tex__t_requests_pipe_lsu_mem_global_op_ld.sum,launch__registers_per_thread,sm__sass_average_branch_targets_threads_uniform.pct --target-processes all --kernel-id "::sigmaKin:" --kernel-base mangled $exe $args | egrep '(sigmaKin|registers| sm|l1tex)' | tr "\n" " " | awk -vtag="[$args]" '{print $1, $2, $3, $16"s", $17";", $19"s", $20, tag}'
  done
  set +x
}

if nvidia-smi -L > /dev/null 2>&1; then gpuTxt="$(nvidia-smi -L | wc -l)x $(nvidia-smi -L | awk '{print $3,$4}' | sort -u)"; else gpuTxt=none; fi
unamep=$(uname -p)
if [ "${unamep}" == "ppc64le" ]; then 
  cpuTxt=$(cat /proc/cpuinfo | grep ^machine | awk '{print substr($0,index($0,"Power"))", "}')$(cat /proc/cpuinfo | grep ^cpu | head -1 | awk '{print substr($0,index($0,"POWER"))}')
else
  cpuTxt=$(cat /proc/cpuinfo | grep '^model name' | head -1 | awk '{i0=index($0,"Intel"); if (i0==0) i0=index($0,"AMD"); i1=index($0," @"); if (i1>0) {print substr($0,i0,i1-i0)} else {print substr($0,i0)}}')
fi
echo -e "On $HOSTNAME [CPU: $cpuTxt] [GPU: $gpuTxt]:"

# Workaround for reading of data files
pushd $topdir/epochX/cudacpp/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum >& /dev/null

lastExe=
for exe in $exes; do
  ###if [ ! -f $exe ]; then continue; fi
  if [ ! -f $exe ]; then echo "Not found: $exe"; continue; fi
  if [ "${exe%%/gcheck*}" != "${exe}" ] && [ "$gpuTxt" == "none" ]; then continue; fi
  if [ "${exe%%/gg_ttgg*}" != "${exe}" ]; then 
    # This is a good GPU middle point: tput is 1.5x lower with "32 256 1", only a few% higher with "128 256 1"
    exeArgs="-p 64 256 1"
    ncuArgs="-p 64 256 1"
  elif [ "${exe%%/gg_tt*}" != "${exe}" ]; then 
    exeArgs="-p 2048 256 1"
    ncuArgs="-p 2048 256 1"
  else # eemumu
    exeArgs="-p 2048 256 12"
    ncuArgs="-p 2048 256 1"
  fi
  if [ "$(basename $exe)" != "$lastExe" ]; then
    echo "========================================================================="
    lastExe=$(basename $exe)
  else
    echo "-------------------------------------------------------------------------"
  fi
  unset OMP_NUM_THREADS
  runExe $exe "$exeArgs"
  if [ "${exe%%/check*}" != "${exe}" ]; then 
    obj=${exe%%/check*}/CPPProcess.o; $scrdir/simdSymSummary.sh -stripdir ${obj}
    if [ "${omp}" == "1" ]; then 
      echo "-------------------------------------------------------------------------"
      export OMP_NUM_THREADS=$(nproc --all)
      runExe $exe "$exeArgs"
    fi
  elif [ "${exe%%/gcheck*}" != "${exe}" ]; then 
    runNcu $exe "$ncuArgs"
    if [ "${div}" == "1" ]; then runNcuDiv $exe; fi
    if [ "${req}" == "1" ]; then runNcuReq $exe "$ncuArgs"; fi
  fi
done
echo "========================================================================="

# Workaround for reading of data files
popd >& /dev/null
printf "\nTEST COMPLETED\n"

