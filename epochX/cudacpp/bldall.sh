#!/bin/bash

if [ "$1" == "" ] || [ "$2" != "" ]; then
  echo "Usage: $0 <process_directory>"
  exit 1
fi
spdir=${1}/SubProcesses

if [ ! -d ${spdir} ]; then
  echo "ERROR! Directory not found: ${spdir}"
  exit 1
fi
cd ${spdir}

START=$(date +%s)

TMP=${START}
for pdir in P*; do
  echo "--------------------------------------------------------------------------------"
  cd ${pdir}
  pwd
  sleep 1
  make -j bldall
  cd - > /dev/null
  END=$(date +%s)
  echo "ELAPSED: $((END-TMP)) seconds"
  TMP=${END}
done
echo "================================================================================"

END=$(date +%s)
echo "ELAPSED: $((END-START)) seconds"
