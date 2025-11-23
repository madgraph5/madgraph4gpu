#!/bin/bash
# Copyright (C) 2020-2025 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: A. Valassi (Oct 2025) for the MG5aMC CUDACPP plugin.
# Further modified by: A. Valassi (2025) for the MG5aMC CUDACPP plugin.

set -e # exit on error

OUTFILE=""
scrdir=$(cd $(dirname ${0}); pwd -P)

if [ "${HOSTNAME}" == "itscrd-a100.cern.ch" ]; then
  node=a100
elif [ "${HOSTNAME}" == "itgold91.cern.ch" ]; then
  node=gold
else
  node=rd90
fi

function runDirFpBld()
{
  if [ "$3" == "" ] || [ "$5" != "" ]; then echo "Usage: $0 <dir> <fptype> <bld> [<checkexe_args>]"; exit 1; fi
  dir=$1
  fp=$2
  bld0=$3
  arg0=$4
  cd $dir
  tmp1=colortimer_TMP1.txt
  tmp2=colortimer_TMP2.txt
  # Configure DCDIAG and BLAS in CUDA
  unset CUDACPP_RUNTIME_BLASCOLORSUM
  unset CUDACPP_RUNTIME_CUBLASTF32TENSOR
  if [ "${bld0}" == "dcd0-blas-TC" ]; then
    bld=cuda; dcd=_dcd0; export CUDACPP_RUNTIME_BLASCOLORSUM=1; export CUDACPP_RUNTIME_CUBLASTF32TENSOR=1
  elif [ "${bld0}" == "dcd1-blas-TC" ]; then
    bld=cuda; dcd=_dcd1; export CUDACPP_RUNTIME_BLASCOLORSUM=1; export CUDACPP_RUNTIME_CUBLASTF32TENSOR=1
  elif [ "${bld0}" == "dcd0-blas" ]; then
    bld=cuda; dcd=_dcd0; export CUDACPP_RUNTIME_BLASCOLORSUM=1
  elif [ "${bld0}" == "dcd1-blas" ]; then
    bld=cuda; dcd=_dcd1; export CUDACPP_RUNTIME_BLASCOLORSUM=1
  elif [ "${bld0}" == "dcd0" ]; then
    bld=cuda; dcd=_dcd0
  elif [ "${bld0}" == "dcd1" ]; then
    bld=cuda; dcd=_dcd1
  else # C++ build
    bld=${bld0}; dcd=
  fi
  # Check.exe arguments (NB use grid size where fptype=f reaches ~peak throughput)
  proc=$(basename $(cd $(pwd -P)/../..; pwd -P))
  proc=${proc/.sa}
  if [ "${arg0}" != "" ]; then
    argCpu="${arg0}"
    argGpu="${arg0}"
  elif [ "${proc#gg_ttgggg}" != "${proc}" ]; then
    argCpu="1 32 1" # 16 32 1?
    argGpu="1 32 1" # 16 32 10?
  else
    echo "ERROR! Unknown proc ${proc}"; exit 1
  fi
  if [ "${bld}" == "cuda" ]; then arg=${argGpu}; else arg=${argCpu}; fi
  # Check.exe command
  if [ "${bld}" == "cuda" ]; then cc=cuda; else cc=cpp; fi
  cmd="./build.${bld}_${fp}_inl0_hrd0${dcd}/check_${cc}.exe -p ${arg}"
  # Use common random numbers to allow comparisons between CUDA/a100 and SIMD/gold
  cmd="${cmd} --common"
  # Skip helicity filtering (this test is designed for gg_ttgggg+)
  export CUDACPP_RUNTIME_GOODHELICITIES=ALL
  # Banner
  echo "PROC=${proc} FPTYPE=${fp} BLD=${bld0} (ARG='${arg}')"
  # Run without timer (check timer overhead)
  unset CUDACPP_RUNTIME_COLORTIMER
  ${cmd} > ${tmp1}
  sk0=$(cat ${tmp1} | awk '/SigmaKin/{print $4}')
  # Run with timer
  export CUDACPP_RUNTIME_COLORTIMER=1
  ${cmd} > ${tmp2}
  ne=$(cat ${tmp2} | awk '/TotalEventsComputed/{print $3}')
  tp=$(cat ${tmp2} | awk '/EvtsPerSec\[MECalcOnly\]/{print $5}')
  sk=$(cat ${tmp2} | awk '/SigmaKin/{print $4}')
  me=$(cat ${tmp2} | awk '/TOTALMEKCMES/{print $3}')
  ja=$(cat ${tmp2} | awk '/CALCJAMPS/{print $4}')
  cs=$(cat ${tmp2} | awk '/23  COLORSUM/{print $4}')
  # Dump timer overhead
  if [ -z ${CUDACPP_RUNTIME_USECHRONOTIMERS+x} ]; then ch=0; else ch=1; fi # check if set even if empty (see https://stackoverflow.com/a/13864829)
  python3 -c "sk=${sk}; sk0=${sk0}; ch=${ch}; print('-> SK with / without timers: %6f / %6f (x%6.4f) [chronotimers=%i]'%(sk,sk0,sk/sk0,ch))"
  # Dump throughput
  python3 -c "sk=${sk}; ne=${ne}; tp=${tp}; print('-> Throughput: %6f/s ( %3i / %6fs )'%(tp,ne,sk))"
  # Dump colortimer results
  python3 -c "me=${me}; ja=${ja}; cs=${cs}; print('-> Jamps    / MEs : %6f / %6f (%7.4f%%)'%(ja,me,ja/me*100))"
  python3 -c "me=${me}; ja=${ja}; cs=${cs}; print('-> ColorSum / MEs : %6f / %6f (%7.4f%%)'%(cs,me,cs/me*100))"
  # Dump physics results
  cat ${tmp2} | awk '/MeanMatrixElemValue/{print "->", $1, ":", $4}'
  # Save colortimer results to file
  if [ "${OUTFILE}" != "" ]; then
    cspct=$(python3 -c "me=${me}; cs=${cs}; print('%7.4f'%(cs/me*100))")
    varg=(${arg})
    printf "%-8s  %-1s  %-12s  %4s %3s %3s  %7s\n" ${proc} ${fp} ${bld0} ${varg[0]} ${varg[1]} ${varg[2]} ${cspct} >> ${OUTFILE}
  fi   
  # Clean up
  unset CUDACPP_RUNTIME_CUBLASTF32TENSOR
  unset CUDACPP_RUNTIME_BLASCOLORSUM
  unset CUDACPP_RUNTIME_COLORTIMER
  \rm ${tmp1}
  \rm ${tmp2}
}

function runggtt4gCuda()
{
  if [ "$1" != "" ]; then echo "Usage $0"; exit 1; fi
  fp=m
  ###OUTFILE=${scrdir}/cs_${node}_ggttgggg_scan_${fp}.txt; \rm -f ${OUTFILE} # save results to file
  dir100=${scrdir}/../../gg_ttgggg.dpg100dpf100.sa/SubProcesses/P1_Sigma_sm_gg_ttxgggg
  dir1000=${scrdir}/../../gg_ttgggg.dpg1000dpf1000.sa/SubProcesses/P1_Sigma_sm_gg_ttxgggg
  ###for carg in "1 32 1"; do # QUICK TEST
  for carg in "1 32 1" "4 32 1" "16 32 1"; do
    ###for pair in "$dir100 dcd1"; do # QUICK TEST
    for pair in "$dir1000 dcd0" "$dir100 dcd0" "$dir100 dcd1" "$dir100 dcd1-blas" "$dir100 dcd1-blas-TC"; do 
      dir=$(echo $pair | cut -d' ' -f1)
      bld0=$(echo $pair | cut -d' ' -f2)
      cd $dir
      runDirFpBld . ${fp} ${bld0} "${carg}"
    done
    echo
  done
  ###if [ "${OUTFILE}" != "" ]; then echo; echo "Result file: ${OUTFILE}"; cat ${OUTFILE}; fi
}

function runggtt4gSimd()
{
  if [ "$1" != "" ]; then echo "Usage $0"; exit 1; fi
  fp=m
  ###OUTFILE=${scrdir}/cs_${node}_ggttgggg_scan_${fp}.txt; \rm -f ${OUTFILE} # save results to file
  dir100=${scrdir}/../../gg_ttgggg.dpg100dpf100.sa/SubProcesses/P1_Sigma_sm_gg_ttxgggg
  dir1000=${scrdir}/../../gg_ttgggg.dpg1000dpf1000.sa/SubProcesses/P1_Sigma_sm_gg_ttxgggg
  ###for dir in $dir100; do # QUICK TEST
  for dir in $dir100 $dir1000; do
    ###for bld0 in 512z; do # QUICK TEST
    for bld0 in none sse4 avx2 512y 512z; do
      cd $dir
      runDirFpBld . ${fp} ${bld0} "1 32 1"
    done
    echo
  done
  ###if [ "${OUTFILE}" != "" ]; then echo; echo "Result file: ${OUTFILE}"; cat ${OUTFILE}; fi
}

# TEST INDIVIDUAL COMPONENTS
# See https://unix.stackexchange.com/a/129077
###runDirFpBld "$@"

# FOR THE PAPER: GGTTGGGG CUDA
#runggtt4gCuda |& tee ${scrdir}/cs_${node}_ggtt4gCuda.txt

# FOR THE PAPER: GGTTGGGG SIMD
#runggtt4gSimd |& tee ${scrdir}/cs_${node}_ggtt4gSimd.txt
