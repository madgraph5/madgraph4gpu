#/bin/bash

eemumu=0
ggtt=0
ggttgg=0

function usage()
{
  echo "Usage: $0 <processes [-eemumu] [-ggtt] [-ggttgg]>"
  exit 1
}

while [ "$1" != "" ]; do  
  if [ "$1" == "-eemumu" ]; then
    eemumu=1
    shift
  elif [ "$1" == "-ggtt" ]; then
    ggtt=1
    shift
  elif [ "$1" == "-ggttgg" ]; then
    ggttgg=1
    shift
  else
    usage
  fi
done

# Check that at least one process has been selected
processes=
if [ "${ggttgg}" == "1" ]; then processes="gg_ttgg $processes"; fi
if [ "${ggtt}" == "1" ]; then processes="gg_tt $processes"; fi
if [ "${eemumu}" == "1" ]; then processes="ee_mumu $processes"; fi
if [ "${processes}" == "" ]; then usage; fi

echo "processes: ${processes}"

cd $(dirname $0)
for proc in ${processes}; do
  echo "------------------------------------------------------------------"
  cmdfile=$(mktemp)
  ./diffCode.sh ../${proc}.auto ../${proc} -r | grep diff | awk '{print "cp", $(NF-1), $NF}' > ${cmdfile}
  cat ${cmdfile}
  source ${cmdfile}
  # FIXME: the script is not handling files which only exist in one of the two directories
  echo -e "\nPENDING DIFFERENCES:"
  echo "./diffCode.sh ../${proc}.auto ../${proc} -r"
  ./diffCode.sh ../${proc}.auto ../${proc} -r
  if [ "$?" == "0" ]; then echo "(no differences)"; fi
done
