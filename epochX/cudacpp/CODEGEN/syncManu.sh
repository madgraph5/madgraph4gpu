#/bin/bash

eemumu=0
ggtt=0
ggttg=0
ggttgg=0
ggttggg=0

function usage()
{
  echo "Usage: $0 <processes [-eemumu] [-ggtt] [-ggttg] [-ggttgg] [-ggttggg]>"
  exit 1
}

while [ "$1" != "" ]; do  
  if [ "$1" == "-eemumu" ]; then
    eemumu=1
    shift
  elif [ "$1" == "-ggtt" ]; then
    ggtt=1
    shift
  elif [ "$1" == "-ggttg" ]; then
    ggttg=1
    shift
  elif [ "$1" == "-ggttgg" ]; then
    ggttgg=1
    shift
  elif [ "$1" == "-ggttggg" ]; then
    ggttggg=1
    shift
  else
    usage
  fi
done

# Select processes
processes=
if [ "${ggtt}" == "1" ]; then processes="gg_tt $processes"; fi
if [ "${ggttg}" == "1" ]; then processes="gg_ttg $processes"; fi
if [ "${ggttgg}" == "1" ]; then processes="gg_ttgg $processes"; fi
if [ "${ggttggg}" == "1" ]; then processes="gg_ttggg $processes"; fi
if [ "${eemumu}" == "1" ]; then processes="ee_mumu $processes"; fi

# Optional hack to hardcode additional processes
###processes="heft_gg_h $processes"

# Check that at least one process has been selected
if [ "${processes}" == "" ]; then usage; fi
echo "processes: ${processes}"

cd $(dirname $0)
for proc in ${processes}; do
  cmdfile=$(mktemp)
  echo "------------------------------------------------------------------"
  ./diffCode.sh ../${proc}.auto ../${proc} -r | grep ^diff | awk '{print "cp -dp", $(NF-1), $NF}' > ${cmdfile}
  cat ${cmdfile}
  source ${cmdfile}
  echo -e "\nPENDING DIFFERENCES (before copying to manual any new files added in auto):"
  echo "./diffCode.sh ../${proc}.auto ../${proc} -r"
  ./diffCode.sh ../${proc}.auto ../${proc} -r
  if [ "$?" == "0" ]; then echo "(no differences)"; continue; fi
  echo "------------------------------------------------------------------"
  ./diffCode.sh ../${proc}.auto ../${proc} -r | awk '{if(index($3,".auto")>0){print "cp -dp", $3$4, gensub(".auto","","g",$3$4)}}' | tr ':' '/' > ${cmdfile}
  cat ${cmdfile}
  source ${cmdfile}
  # FIXME: the script is not handling files which only exist in manual (i.e. have been removed in auto)
  echo -e "\nPENDING DIFFERENCES (after copying to manual any new files added in auto):"
  echo "./diffCode.sh ../${proc}.auto ../${proc} -r"
  ./diffCode.sh ../${proc}.auto ../${proc} -r
  if [ "$?" == "0" ]; then echo "(no differences)"; continue; fi
done
