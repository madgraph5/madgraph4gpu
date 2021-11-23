#!/bin/bash

cd $(dirname $0)

function usage()
{
  echo "Usage: $0 <procs (-eemumu|-ggtt|-ggttgg)> [-auto|-autoonly] [-flt|-fltonly] [-inl|-inlonly] [-makeonly] [-makeclean]"
  exit 1
}

procs=
eemumu=
ggtt=
ggttgg=
suffs="manu"
fptypes="d"
helinls="0"
steps="make test"
for arg in $*; do
  if [ "$arg" == "-eemumu" ]; then
    if [ "$eemumu" == "" ]; then procs+=${procs:+ }${arg}; fi
    eemumu=$arg
  elif [ "$arg" == "-ggtt" ]; then
    if [ "$ggtt" == "" ]; then procs+=${procs:+ }${arg}; fi
    ggtt=$arg
  elif [ "$arg" == "-ggttgg" ]; then
    if [ "$ggttgg" == "" ]; then procs+=${procs:+ }${arg}; fi
    ggttgg=$arg
  elif [ "$arg" == "-auto" ]; then
    if [ "${suffs}" == "auto" ]; then echo "ERROR! Options -auto and -autoonly are incompatible"; usage; fi
    suffs="manu auto"
  elif [ "$arg" == "-autoonly" ]; then
    if [ "${suffs}" == "manu auto" ]; then echo "ERROR! Options -auto and -autoonly are incompatible"; usage; fi
    suffs="auto"
  elif [ "$arg" == "-flt" ]; then
    if [ "${fptypes}" == "f" ]; then echo "ERROR! Options -flt and -fltonly are incompatible"; usage; fi
    fptypes="d f"
  elif [ "$arg" == "-fltonly" ]; then
    if [ "${fptypes}" == "d f" ]; then echo "ERROR! Options -flt and -fltonly are incompatible"; usage; fi
    fptypes="f"
  elif [ "$arg" == "-inl" ]; then
    if [ "${helinls}" == "1" ]; then echo "ERROR! Options -inl and -inlonly are incompatible"; usage; fi
    helinls="0 1"
  elif [ "$arg" == "-inlonly" ]; then
    if [ "${helinls}" == "0 1" ]; then echo "ERROR! Options -inl and -inlonly are incompatible"; usage; fi
    helinls="1"
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
  else
    echo "ERROR! Invalid option '$arg'"; usage
  fi  
done

#echo "procs=$procs"
#echo "suffs=$suffs"
#echo "df=$df"
#echo "inl=$inl"
#echo "steps=$steps"
###exit 0

started="STARTED AT $(date)"

for step in $steps; do
  for proc in $procs; do
    for suff in $suffs; do
      auto=; if [ "${suff}" == "auto" ]; then auto=" -autoonly"; fi
      if [ "${proc}" == "-ggtt" ] && [ "${suff}" == "manu" ]; then
        ###printf "\n%80s\n" |tr " " "*"
        ###printf "*** WARNING! ${proc#-}_${suff} does not exist"
        ###printf "\n%80s\n" |tr " " "*"
        continue
      fi
      for fptype in $fptypes; do
        flt=; if [ "${fptype}" == "f" ]; then flt=" -fltonly"; fi
        for helinl in $helinls; do
          inl=; if [ "${helinl}" == "1" ]; then inl=" -inlonly"; fi
          args="${proc}${auto}${flt}${inl}"
          args="${args} -avxall" # avx, fptype and helinl are now supported for all processes
          if [ "${step}" == "makeclean" ]; then
            printf "\n%80s\n" |tr " " "*"
            printf "*** ./throughputX.sh -makecleanonly $args"
            printf "\n%80s\n" |tr " " "*"
            if ! ./throughputX.sh -makecleanonly $args; then exit 1; fi
          elif [ "${step}" == "make" ]; then
            printf "\n%80s\n" |tr " " "*"
            printf "*** ./throughputX.sh -makeonly $args"
            printf "\n%80s\n" |tr " " "*"
            if ! ./throughputX.sh -makeonly $args; then exit 1; fi
          else
            logfile=logs_${proc#-}_${suff}/log_${proc#-}_${suff}_${fptype}_inl${helinl}.txt
            printf "\n%80s\n" |tr " " "*"
            printf "*** ./throughputX.sh $args | tee $logfile"
            printf "\n%80s\n" |tr " " "*"
            mkdir -p $(dirname $logfile)
            ./throughputX.sh $args -gtest | tee $logfile
          fi
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
