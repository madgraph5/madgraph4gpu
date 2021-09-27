#!/bin/bash

cd $(dirname $0)
if [ "$*" == "" ]; then args="-eemumu -ggttgg -auto"; else args="$*"; fi
if ! ./throughputX.sh -makeonly $args; then exit 1; fi
./throughputX.sh $args | tee throughputX_log.txt
