#!/bin/bash
# Copyright (C) 2020-2025 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: A. Valassi (Oct 2025) for the MG5aMC CUDACPP plugin.
# Further modified by: A. Valassi (2025) for the MG5aMC CUDACPP plugin.

set -e # exit on error

OUTFILE=""
scrdir=$(cd $(dirname ${0}); pwd -P)

function runDirFpBld()
{
  if [ "$3" == "" ] || [ "$5" != "" ]; then echo "Usage $0 <dir> <fptype> <bld> [<checkexe_args>]"; exit 1; fi
  dir=$1
  fp=$2
  bld0=$3
  arg0=$4
  cd $1
  tmp=colortimer_TMP.txt
  # Enable BLAS in CUDA?
  unset CUDACPP_RUNTIME_BLASCOLORSUM
  unset CUDACPP_RUNTIME_CUBLASTF32TENSOR
  if [ "${bld0}" == "cuda-blas-TC" ]; then
    bld=cuda; export CUDACPP_RUNTIME_BLASCOLORSUM=1; export CUDACPP_RUNTIME_CUBLASTF32TENSOR=1
  elif [ "${bld0}" == "cuda-blas" ]; then
    bld=cuda; export CUDACPP_RUNTIME_BLASCOLORSUM=1
  else
    bld=${bld0}
  fi
  # Check.exe arguments (NB use grid size where fptype=f reaches ~peak throughput)
  proc=$(basename $(cd $(pwd -P)/../..; pwd -P))
  proc=${proc/.mad}
  if [ "${arg0}" != "" ]; then
    argCpu="${arg0}"
    argGpu="${arg0}"
  elif [ "${proc}" == "gg_tt" ]; then
    argCpu="2048 32 1"
    argGpu="2048 32 10"
  elif [ "${proc}" == "gg_ttg" ]; then
    argCpu="1024 32 1"
    argGpu="1024 32 10"
  elif [ "${proc}" == "gg_ttgg" ]; then
    argCpu="256 32 1"
    argGpu="256 32 10"
  elif [ "${proc}" == "gg_ttggg" ]; then
    if [ "${skipCuda}" == "" ]; then
      argCpu="16 32 1"
    else
      argCpu="4 32 1"
    fi
    ###argGpu="4 32 10" # blas always loses
    ###argGpu="8 32 10" # blas always loses
    argGpu="16 32 10" # blas beats kernel for fptype=d (NB for fptype=f, "4 32 10" has much lower tput!)
  else
    echo "ERROR! Unknown proc ${proc}"; exit 1
  fi
  if [ "${bld}" == "cuda" ]; then arg=${argGpu}; else arg=${argCpu}; fi
  # Check.exe command
  if [ "${bld}" == "cuda" ]; then cc=cuda; else cc=cpp; fi
  cmd="./build.${bld}_${fp}_inl0_hrd0/check_${cc}.exe -p ${arg}"
  # Banner
  echo
  echo "PROC=${proc} FPTYPE=${fp} BLD=${bld0} (ARG='${arg}')"
  # Run without timer (check timer overhead)
  unset CUDACPP_RUNTIME_COLORTIMER
  ${cmd} > ${tmp}
  sk0=$(cat ${tmp} | awk '/SigmaKin/{print $4}')
  # Run with timer
  export CUDACPP_RUNTIME_COLORTIMER=1
  ${cmd} > ${tmp}
  sk=$(cat ${tmp} | awk '/SigmaKin/{print $4}')
  me=$(cat ${tmp} | awk '/TOTALMEKCMES/{print $3}')
  ja=$(cat ${tmp} | awk '/CALCJAMPS/{print $4}')
  cs=$(cat ${tmp} | awk '/23  COLORSUM/{print $4}')
  # Dump timer overhead
  if [ -z ${CUDACPP_RUNTIME_USECHRONOTIMERS+x} ]; then ch=0; else ch=1; fi # check if set even if empty (see https://stackoverflow.com/a/13864829)
  python3 -c "sk=${sk}; sk0=${sk0}; ch=${ch}; print('-> SK with / without timers: %6f / %6f (x%6.4f) [chronotimers=%i]'%(sk,sk0,sk/sk0,ch))"
  # Dump colortimer results
  python3 -c "me=${me}; ja=${ja}; cs=${cs}; print('-> Jamps    / MEs : %6f / %6f (%7.4f%%)'%(ja,me,ja/me*100))"
  python3 -c "me=${me}; ja=${ja}; cs=${cs}; print('-> ColorSum / MEs : %6f / %6f (%7.4f%%)'%(cs,me,cs/me*100))"
  # Dump physics results
  cat ${tmp} | awk '/MeanMatrixElemValue/{print "->", $1, ":", $4}'
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
  \rm ${tmp}
}

function runDirFp()
{
  if [ "$2" == "" ] || [ "$3" != "" ]; then echo "Usage $0 <dir> <fptype>"; exit 1; fi
  dir=$1
  fp=$2
  cd $1
  if [ "${skipCuda}" == "" ]; then
    if [ "${HOSTNAME}" == "itscrd-a100.cern.ch" ]; then
      runDirFpBld . ${fp} cuda-blas-TC
    fi  
    runDirFpBld . ${fp} cuda-blas
    runDirFpBld . ${fp} cuda
  fi
  runDirFpBld . ${fp} none
  runDirFpBld . ${fp} sse4
  runDirFpBld . ${fp} avx2
  if [ "${HOSTNAME}" != "itscrd-a100.cern.ch" ]; then
    runDirFpBld . ${fp} 512y
    runDirFpBld . ${fp} 512z
  fi
}

function runDir()
{
  if [ "$1" == "" ] || [ "$2" != "" ]; then echo "Usage $0 <dir>"; exit 1; fi
  dir=$1
  cd $1
  runDirFp . m
  runDirFp . d
  runDirFp . f
}

function runAll()
{
  if [ "${HOSTNAME}" == "itscrd-a100.cern.ch" ]; then node=a100; else node=rd90; fi
  OUTFILE=${scrdir}/cs_${node}_allproc_dmf.txt; \rm -f ${OUTFILE} # save results to file
  runDir ${scrdir}/../gg_tt.mad/SubProcesses/P1_gg_ttx
  runDir ${scrdir}/../gg_ttg.mad/SubProcesses/P1_gg_ttxg
  runDir ${scrdir}/../gg_ttgg.mad/SubProcesses/P1_gg_ttxgg
  runDir ${scrdir}/../gg_ttggg.mad/SubProcesses/P1_gg_ttxggg
  if [ "${OUTFILE}" != "" ]; then echo; echo "Result file: ${OUTFILE}"; cat ${OUTFILE}; fi
}

function buildDir()
{
  if [ "$1" == "" ] || [ "$2" != "" ]; then echo "Usage $0 <dir>"; exit 1; fi
  dir=$1
  cd $1
  make -j -f cudacpp.mk cleanall
  make -j -f cudacpp.mk bldall FPTYPE=m
  make -j -f cudacpp.mk bldall FPTYPE=d
  make -j -f cudacpp.mk bldall FPTYPE=f
}

function buildAll()
{
  cd ${scrdir}
  buildDir ${scrdir}/../gg_tt.mad/SubProcesses/P1_gg_ttx
  buildDir ${scrdir}/../gg_ttg.mad/SubProcesses/P1_gg_ttxg
  buildDir ${scrdir}/../gg_ttgg.mad/SubProcesses/P1_gg_ttxgg
  buildDir ${scrdir}/../gg_ttggg.mad/SubProcesses/P1_gg_ttxggg
}

function runggttgggFp()
{
  if [ "$1" == "" ] || [ "$2" != "" ]; then echo "Usage $0 <fp>"; exit 1; fi
  fp=$1
  if [ "${HOSTNAME}" == "itscrd-a100.cern.ch" ]; then node=a100; else node=rd90; fi
  OUTFILE=${scrdir}/cs_${node}_ggttggg_scan_${fp}.txt; \rm -f ${OUTFILE} # save results to file
  dir=${scrdir}/../gg_ttggg.mad/SubProcesses/P1_gg_ttxggg
  cd $dir
  ###for carg in "4 32 1"; do # QUICK TEST
  ###for carg in "4 32 16" "8 32 8" "16 32 4" "32 32 2" "64 32 1"; do
  ###for carg in "4 32 32" "8 32 16" "16 32 8" "32 32 4" "64 32 2" "128 32 1"; do
  ###for carg in "4 32 64" "8 32 32" "16 32 16" "32 32 8" "64 32 4" "128 32 2" "256 32 1"; do
  for carg in "4 32 128" "8 32 64" "16 32 32" "32 32 16" "64 32 8" "128 32 4" "256 32 2" "512 32 1"; do
    if [ "${HOSTNAME}" == "itscrd-a100.cern.ch" ]; then
      runDirFpBld . ${fp} cuda-blas-TC "${carg}"
    fi
    runDirFpBld . ${fp} cuda-blas "${carg}"
    runDirFpBld . ${fp} cuda "${carg}"
  done
  if [ "${OUTFILE}" != "" ]; then echo; echo "Result file: ${OUTFILE}"; cat ${OUTFILE}; fi
}

# SKIP CUDA?
skipCuda=

# TEST INDIVIDUAL COMPONENTS
###buildDir $*
###runDirFpBld $*
###runDirFp $*
###runDir $*

# FOR THE PAPER: BUILD ALL PROCESSES
###buildAll

# FOR THE PAPER: ALL PROCESSES
#runAll

# FOR THE PAPER: GGTTGGG SCANS
#runggttgggFp f
#runggttgggFp m
#runggttgggFp d

# FOR THE PAPER: GGTTGGG/SIMD
skipCuda=1; cd ${scrdir}/../gg_ttggg.mad/SubProcesses/P1_gg_ttxggg; runDir . | tee ${scrdir}/simd_gold91_raw.txt; cd -
${scrdir}/simdparser.py ${scrdir}/simd_gold91_raw.txt | tee  ${scrdir}/simd_gold91_summary.txt
