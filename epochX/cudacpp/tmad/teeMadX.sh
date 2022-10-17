#!/bin/bash

scrdir=$(cd $(dirname $0); pwd)
bckend=$(basename $(cd $scrdir; cd ..; pwd)) # cudacpp or alpaka
cd $scrdir

function usage()
{
  echo "Usage: $0 <processes [-eemumu][-ggtt][-ggttg][-ggttgg][-ggttggg]> [-flt|-fltonly|-mix] [-makeonly] [-makeclean] [-rmrdat] [+10x] [-checkonly]" > /dev/stderr
  exit 1
}

procs=
eemumu=
ggtt=
ggttg=
ggttgg=
ggttggg=

suffs="mad"
fptypes="d"
helinls="0"
hrdcods="0"

steps="make test"

rmrdat=

###deb=
deb=" -d" # optional debug mode

makej=
dlp=

add10x=

checkonly=

for arg in $*; do
  if [ "$arg" == "-eemumu" ]; then
    if [ "$eemumu" == "" ]; then procs+=${procs:+ }${arg}; fi
    eemumu=$arg
  elif [ "$arg" == "-ggtt" ]; then
    if [ "$ggtt" == "" ]; then procs+=${procs:+ }${arg}; fi
    ggtt=$arg
  elif [ "$arg" == "-ggttg" ]; then
    if [ "$ggttg" == "" ]; then procs+=${procs:+ }${arg}; fi
    ggttg=$arg
  elif [ "$arg" == "-ggttgg" ]; then
    if [ "$ggttgg" == "" ]; then procs+=${procs:+ }${arg}; fi
    ggttgg=$arg
  elif [ "$arg" == "-ggttggg" ]; then
    if [ "$ggttggg" == "" ]; then procs+=${procs:+ }${arg}; fi
    ggttggg=$arg
  elif [ "$arg" == "-flt" ]; then
    if [ "${fptypes}" != "d" ] && [ "${fptypes}" != "d f" ]; then echo "ERROR! Options -flt, -fltonly and -mix are incompatible"; usage; fi
    fptypes="d f"
  elif [ "$arg" == "-fltonly" ]; then
    if [ "${fptypes}" != "d" ] && [ "${fptypes}" != "f" ]; then echo "ERROR! Options -flt, -fltonly and -mix are incompatible"; usage; fi
    fptypes="f"
  elif [ "$arg" == "-mix" ]; then
    if [ "${fptypes}" != "d" ] && [ "${fptypes}" != "d f m" ]; then echo "ERROR! Options -flt, -fltonly and -mix are incompatible"; usage; fi
    fptypes="d f m"
  #elif [ "$arg" == "-inl" ]; then
  #  if [ "${helinls}" == "1" ]; then echo "ERROR! Options -inl and -inlonly are incompatible"; usage; fi
  #  helinls="0 1"
  #elif [ "$arg" == "-inlonly" ]; then
  #  if [ "${helinls}" == "0 1" ]; then echo "ERROR! Options -inl and -inlonly are incompatible"; usage; fi
  #  helinls="1"
  #elif [ "$arg" == "-hrd" ]; then
  #  if [ "${hrdcods}" == "1" ]; then echo "ERROR! Options -hrd and -hrdonly are incompatible"; usage; fi
  #  hrdcods="0 1"
  #elif [ "$arg" == "-hrdonly" ]; then
  #  if [ "${hrdcods}" == "0 1" ]; then echo "ERROR! Options -hrd and -hrdonly are incompatible"; usage; fi
  #  hrdcods="1"
  elif [ "$arg" == "-makeonly" ]; then
    if [ "${steps}" == "make test" ]; then
      steps="make"
    elif [ "${steps}" == "makeclean make test" ]; then
      steps="makeclean make"
    fi
  elif [ "$arg" == "-makeclean" ]; then
    if [ "${steps}" == "make test" ]; then
      steps="makeclean make test"
    elif [ "${steps}" == "make" ]; then
      steps="makeclean make"
    fi
  #elif [ "$arg" == "-makej" ]; then
  #  makej=-makej
  elif [ "$arg" == "-rmrdat" ]; then
    rmrdat=" $arg"
  elif [ "$arg" == "+10x" ]; then
    add10x="$add10x $arg"
  elif [ "$arg" == "-checkonly" ]; then
    checkonly=" $arg"
  else
    echo "ERROR! Invalid option '$arg'"; usage
  fi  
done

# Check that at least one process has been selected
if [ "${procs}" == "" ]; then usage; fi

status=0
started="STARTED AT $(date)"

for step in $steps; do
  for proc in $procs; do
    for suff in $suffs; do
      for fptype in $fptypes; do
        flt=; if [ "${fptype}" == "f" ]; then flt=" -fltonly"; elif [ "${fptype}" == "m" ]; then flt=" -mixonly"; fi
        for helinl in $helinls; do
          inl=; if [ "${helinl}" == "1" ]; then inl=" -inlonly"; fi
          for hrdcod in $hrdcods; do
            hrd=; if [ "${hrdcod}" == "1" ]; then hrd=" -hrdonly"; fi
            args="${proc}${flt}${inl}${hrd}${deb}${rmrdat}${add10x}${checkonly} ${dlp}"
            ###args="${args} -avxall" # avx, fptype, helinl and hrdcod are now supported for all processes
            if [ "${step}" == "makeclean" ]; then
              printf "\n%80s\n" |tr " " "*"
              printf "*** ./madX.sh -makecleanonly $args"
              printf "\n%80s\n" |tr " " "*"
              if ! ./madX.sh -makecleanonly $args; then exit 1; fi
            elif [ "${step}" == "make" ]; then
              printf "\n%80s\n" |tr " " "*"
              printf "*** ./madX.sh -makeonly ${makej} $args"
              printf "\n%80s\n" |tr " " "*"
              if ! ./madX.sh -makeonly ${makej} $args; then exit 1; fi
            else
              logfile=logs_${proc#-}_${suff}/log_${proc#-}_${suff}_${fptype}_inl${helinl}_hrd${hrdcod}.txt
              if [ "${rndgen}" != "" ]; then logfile=${logfile%.txt}_${rndgen#-}.txt; fi
              if [ "${rmbsmp}" != "" ]; then logfile=${logfile%.txt}_${rmbsmp#-}.txt; fi
              if [ "${checkonly}" != "" ]; then
                logfileold=${logfile}
                logfile=${logfile}.TMP
              fi
              printf "\n%80s\n" |tr " " "*"
              printf "*** ./madX.sh $args | tee $logfile"
              printf "\n%80s\n" |tr " " "*"
              mkdir -p $(dirname $logfile)
              if ! ./madX.sh $args | tee $logfile; then status=2; fi
              if [ "${checkonly}" != "" ]; then
                ./checkOnlyMerge.sh ${logfileold}
                \rm -f ${logfile}
                echo "Results merged into ${logfileold}"
              fi
            fi
          done
        done
      done
    done
  done
  printf "\n%80s\n" |tr " " "#"
done

ended="ENDED   AT $(date)"
echo
echo "$started"
echo "$ended"
exit $status
