#!/bin/bash

cd $(dirname $0)

if [ "$*" == "" ]; then
  args="-eemumu -ggttgg"
  auto="-auto"
else 
  args=
  auto=
  for arg in $*; do
    if [ "$arg" == "-auto" ] || [ "$arg" == "-autoonly" ]; then
      if [ "$auto" != "" ] && [ "$auto" != "$arg" ]; then
        echo "ERROR! Incompatible options -auto and -autoonly"
        exit 1
      fi
      auto=$arg
    else
      args="$args $arg"
    fi  
  done
fi
echo "args=$args"
echo "auto=$auto"

if ! ./throughputX.sh -makeonly $args $auto; then exit 1; fi

if [ "$auto" != "-autoonly" ]; then
  ./throughputX.sh $args | tee throughputX_log.txt
fi 

if [ "$auto" == "-auto" ] || [ "$auto" == "-autoonly" ]; then
  ./throughputX.sh $args -autoonly | tee throughputX_log_auto.txt
fi
