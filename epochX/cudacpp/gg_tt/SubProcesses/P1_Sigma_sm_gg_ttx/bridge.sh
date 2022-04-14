#!/bin/bash

###for args in "-p 32 32 1" "-p 32 32 2" "-p 32 256 2"; do
###for args in "-p 32 256 2"; do
for args in "-p 1 32 1" "-p 1 32 2" "-p 1 256 2" "-p 16 256 2" "-p 32 256 2"; do
  for brd in "" "--bridge"; do
    for exe in "./check.exe " "./gcheck.exe"; do
      cmd="$exe $brd $args"; printf "%-36s ==> " "$cmd"
      ###cmd="$exe $brd $args --curhst"; printf "%-46s ==> " "$cmd"
      ###cmd="$exe $brd $args --common"; printf "%-46s ==> " "$cmd"
      eval $cmd | grep MeanM
    done
  done
  echo "------------------------------"
done

