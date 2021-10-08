#!/bin/bash

cd $(dirname $0)

if [ "$*" == "" ]; then
  args=""
  eemumu="-eemumu"
  ggttgg="-ggttgg"
  auto="-auto"
else 
  args=
  eemumu=
  ggttgg=
  auto=
  for arg in $*; do
    if [ "$arg" == "-auto" ] || [ "$arg" == "-autoonly" ]; then
      if [ "$auto" != "" ] && [ "$auto" != "$arg" ]; then
        echo "ERROR! Incompatible options -auto and -autoonly"
        exit 1
      fi
      auto=$arg
    elif [ "$arg" == "-eemumu" ]; then
      eemumu=$arg
    elif [ "$arg" == "-ggttgg" ]; then
      ggttgg=$arg
    else
      args="$args $arg"
    fi  
  done
fi

echo "args=$args"
echo "eemumu=$eemumu"
echo "ggttgg=$ggttgg"
echo "auto=$auto"

if ! ./throughputX.sh -makeonly $args $eemumu $ggttgg $auto; then exit 1; fi

if [ "$auto" != "-autoonly" ]; then
  if [ "$eemumu" != "" ]; then ./throughputX.sh $eemumu $args | tee throughputX_log_eemumu.txt; fi
  if [ "$ggttgg" != "" ]; then ./throughputX.sh $ggttgg $args | tee throughputX_log_ggttgg.txt; fi
fi 

if [ "$auto" == "-auto" ] || [ "$auto" == "-autoonly" ]; then
  if [ "$eemumu" != "" ]; then ./throughputX.sh $eemumu $args -autoonly| tee throughputX_log_eemumu_auto.txt; fi
  if [ "$ggttgg" != "" ]; then ./throughputX.sh $ggttgg $args -autoonly| tee throughputX_log_ggttgg_auto.txt; fi
fi
