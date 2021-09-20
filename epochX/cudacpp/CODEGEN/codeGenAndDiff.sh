#!/bin/bash

#--------------------------------------------------------------------------------------

function codeGenAndDiff()
{
  proc=$1
  if [ "${proc}" == "$(basename $SCRDIR)" ]; then return; fi # e.g. skip CODEGEN  
  if [ "${proc}" != "${proc%.BKP}" ]; then return; fi # e.g. skip ee_mumu.BKP
  if [ "${proc}" != "${proc%.NEW}" ]; then return; fi # e.g. skip ee_mumu.NEW
  echo -e "\n================================================================"
  echo -e "\n+++ Generate code for $proc\n"
  # Process-dependent hardcoded configuration
  case "${proc}" in
    ee_mumu)
      cmd="generate e+ e- > mu+ mu-"
      ;;
    *)
      echo "WARNING! Skipping unknown process $proc"
      return
      ;;
  esac
  # Generate code for the specific process
  pushd $MG5AMC_HOME >& /dev/null
  outproc=CODEGEN_${proc}
  \rm -rf ${outproc}*
  echo "${cmd}" >> ${outproc}.mg
  echo "output standalone_${OUTBCK} ${outproc}" >> ${outproc}.mg
  cat  ${outproc}.mg
  { time python3 ./bin/mg5_aMC ${outproc}.mg ; } >& ${outproc}_log.txt
  cat ${outproc}_log.txt | egrep 'INFO: (Try|Creat|Organiz|Process)'
  mv ${outproc}_log.txt ${outproc}/
  popd >& /dev/null
  # Move the newly generated code to the output source code directory
  rm -rf ${OUTDIR}/${proc}.BKP ${OUTDIR}/${proc}.NEW
  cp -dpr ${MG5AMC_HOME}/${outproc} ${OUTDIR}/${proc}.NEW
  echo -e "\nOutput source code has been copied to ${OUTDIR}/${proc}.NEW"
  # Compare the newly generated code to the existing one for the specific process
  echo -e "\n+++ Compare code for $proc\n"
  pushd ${OUTDIR} >& /dev/null
  diff -rs ${proc}.NEW ${proc}
  popd >& /dev/null 
  # Replace the existing code by the newly generated code if required
  if [ "${REPLACE}" == "1" ]; then
    echo -e "\n+++ Replace existing code for $proc (REPLACE=$REPLACE)\n"
    mv ${OUTDIR}/${proc} ${OUTDIR}/${proc}.BKP
    mv ${OUTDIR}/${proc}.NEW ${OUTDIR}/${proc}
    echo -e "Old code moved to ${OUTDIR}/${proc}.BKP"
    echo -e "New code moved to ${OUTDIR}/${proc}"
  fi
}

#--------------------------------------------------------------------------------------

function usage()
{
  echo "Usage: $0 [--noreplace] [<proc1> [... <procN>]]"
  exit 1
}

#--------------------------------------------------------------------------------------

# Print command line usage if required
if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then usage; fi

# Script directory
SCRDIR=$(cd $(dirname $0); pwd)
echo SCRDIR=${SCRDIR}

# Output source code directory for the chosen backend
OUTDIR=$(dirname $SCRDIR) # e.g. epochX/cudacpp if $SCRDIR=epochX/cudacpp/CODEGEN
echo OUTDIR=${OUTDIR}

# Output backend
OUTBCK=$(basename $OUTDIR) # e.g. cudacpp if $OUTDIR=epochX/cudacpp
echo "OUTBCK=${OUTBCK} (uppercase=${OUTBCK^^})"

# Replace code directory and create .BKP? (or alternatively keep code directory in .NEW?)
REPLACE=1
if [ "$1" == "--noreplace" ]; then REPLACE=0; shift; fi
echo REPLACE=${REPLACE}

# Make sure that python3 is installed
if ! python3 --version >& /dev/null; then echo "ERROR! python3 is not installed"; exit 1; fi

# Make sure that $MG5AMC_HOME exists
if [ "$MG5AMC_HOME" == "" ]; then echo "ERROR! MG5AMC_HOME is not defined"; exit 1; fi
echo -e "\nUsing MG5AMC_HOME=$MG5AMC_HOME on $(hostname)\n"
if [ ! -d $MG5AMC_HOME ]; then echo "ERROR! Directory $MG5AMC_HOME does not exist"; exit 1; fi

# Remove fragments from previous runs
rm -rf ${MG5AMC_HOME}/py.py

# Remove and recreate MG5AMC_HOME/PLUGIN
rm -rf ${MG5AMC_HOME}/PLUGIN
mkdir ${MG5AMC_HOME}/PLUGIN
touch ${MG5AMC_HOME}/PLUGIN/__init__.py
cp -dpr ${SCRDIR}/PLUGIN/${OUTBCK^^}_SA_OUTPUT ${MG5AMC_HOME}/PLUGIN/
ls -lR $MG5AMC_HOME/PLUGIN

# Determine the list of processes to generate
procs=$*
###procs="ee_mumu gg_tt gg_ttg gg_ttgg"
if [ "$procs" == "" ] ; then procs=$(cd $OUTDIR; find . -mindepth 1 -maxdepth 1 -type d); fi

# Iterate through the list of processes to generate
for proc in $procs; do
  proc=$(basename $proc)
  codeGenAndDiff $proc
done
