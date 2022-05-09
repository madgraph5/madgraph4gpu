#!/bin/bash

scrdir=$(cd $(dirname $0); pwd)
bckend=$(basename $(cd $scrdir; cd ..; pwd)) # cudacpp or alpaka
cd $scrdir

function usage()
{
  echo "Usage: $0 <processes [-eemumu][-ggtt][-ggttg][-ggttgg][-ggttggg][-heftggh]> [-auto|-autoonly] [-noalpaka] [-flt|-fltonly] [-inl|-inlonly] [-hrd|-hrdonly] [-common|-curhst] [-rmbhst|-bridge] [-makeonly] [-makeclean] [-makej] [-dlp <dyld_library_path>]"
  exit 1
}

procs=
eemumu=
ggtt=
ggttg=
ggttgg=
ggttggg=
heftggh=
suffs="manu"
alpaka=
fptypes="d"
helinls="0"
hrdcods="0"
rndgen=
rmbsmp=
steps="make test"
makej=
dlp=
dlpset=0

for arg in $*; do
  if [ "${dlpset}" == "1" ]; then
    dlpset=2
    dlp="-dlp $arg"
  elif [ "$arg" == "-dlp" ] && [ "${dlpset}" == "0" ]; then
    dlpset=1
  elif [ "$arg" == "-eemumu" ]; then
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
  elif [ "$arg" == "-heftggh" ]; then
    if [ "$heftggh" == "" ]; then procs+=${procs:+ }${arg}; fi
    heftggh=$arg
  elif [ "$arg" == "-auto" ]; then
    if [ "${suffs}" == "auto" ]; then echo "ERROR! Options -auto and -autoonly are incompatible"; usage; fi
    suffs="manu auto"
  elif [ "$arg" == "-autoonly" ]; then
    if [ "${suffs}" == "manu auto" ]; then echo "ERROR! Options -auto and -autoonly are incompatible"; usage; fi
    suffs="auto"
  elif [ "$arg" == "-noalpaka" ]; then
    alpaka=$arg
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
  elif [ "$arg" == "-hrd" ]; then
    if [ "${hrdcods}" == "1" ]; then echo "ERROR! Options -hrd and -hrdonly are incompatible"; usage; fi
    hrdcods="0 1"
  elif [ "$arg" == "-hrdonly" ]; then
    if [ "${hrdcods}" == "0 1" ]; then echo "ERROR! Options -hrd and -hrdonly are incompatible"; usage; fi
    hrdcods="1"
  elif [ "$arg" == "-common" ]; then
    rndgen=$arg
  elif [ "$arg" == "-curhst" ]; then
    rndgen=$arg
  elif [ "$arg" == "-rmbhst" ]; then
    rmbsmp=$arg
  elif [ "$arg" == "-bridge" ]; then
    rmbsmp=$arg
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
  elif [ "$arg" == "-makej" ]; then
    makej=-makej
  else
    echo "ERROR! Invalid option '$arg'"; usage
  fi  
done

# Workaround for MacOS SIP (SystemIntegrity Protection): set DYLD_LIBRARY_PATH In subprocesses
if [ "${dlpset}" == "1" ]; then usage; fi

# Use only the .auto process directories in the alpaka directory
if [ "$bckend" == "alpaka" ]; then
  echo "WARNING! alpaka directory: using .auto process directories only"
  suffs="auto"
fi

#echo "procs=$procs"
#echo "suffs=$suffs"
#echo "df=$df"
#echo "inl=$inl"
#echo "hrd=$hrd"
#echo "steps=$steps"
###exit 0

started="STARTED AT $(date)"

for step in $steps; do
  for proc in $procs; do
    for suff in $suffs; do
      auto=; if [ "${suff}" == "auto" ]; then auto=" -autoonly"; fi
      ###if [ "${proc}" == "-ggtt" ] && [ "${suff}" == "manu" ]; then
      ###  ###printf "\n%80s\n" |tr " " "*"
      ###  ###printf "*** WARNING! ${proc#-}_${suff} does not exist"
      ###  ###printf "\n%80s\n" |tr " " "*"
      ###  continue
      ###fi
      for fptype in $fptypes; do
        flt=; if [ "${fptype}" == "f" ]; then flt=" -fltonly"; fi
        for helinl in $helinls; do
          inl=; if [ "${helinl}" == "1" ]; then inl=" -inlonly"; fi
          for hrdcod in $hrdcods; do
            hrd=; if [ "${hrdcod}" == "1" ]; then hrd=" -hrdonly"; fi
            args="${proc}${auto}${flt}${inl}${hrd} ${dlp}"
            args="${args} ${alpaka}" # optionally disable alpaka tests
            args="${args} ${rndgen}" # optionally use common random numbers or curand on host
            args="${args} ${rmbsmp}" # optionally use rambo or bridge on host
            args="${args} -avxall" # avx, fptype, helinl and hrdcod are now supported for all processes
            if [ "${step}" == "makeclean" ]; then
              printf "\n%80s\n" |tr " " "*"
              printf "*** ./throughputX.sh -makecleanonly $args"
              printf "\n%80s\n" |tr " " "*"
              if ! ./throughputX.sh -makecleanonly $args; then exit 1; fi
            elif [ "${step}" == "make" ]; then
              printf "\n%80s\n" |tr " " "*"
              printf "*** ./throughputX.sh -makeonly ${makej} $args"
              printf "\n%80s\n" |tr " " "*"
              if ! ./throughputX.sh -makeonly ${makej} $args; then exit 1; fi
            else
              logfile=logs_${proc#-}_${suff}/log_${proc#-}_${suff}_${fptype}_inl${helinl}_hrd${hrdcod}.txt
              if [ "${rndgen}" != "" ]; then logfile=${logfile%.txt}_${rndgen#-}.txt; fi
              if [ "${rmbsmp}" != "" ]; then logfile=${logfile%.txt}_${rmbsmp#-}.txt; fi
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
  done
  printf "\n%80s\n" |tr " " "#"
done

ended="ENDED   AT $(date)"
echo
echo "$started"
echo "$ended"
