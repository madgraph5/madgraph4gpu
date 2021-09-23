#!/bin/bash

cd $(dirname $0)
./throughputX.sh $* | tee throughputX_log.txt
