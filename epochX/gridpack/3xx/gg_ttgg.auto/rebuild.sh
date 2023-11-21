#!/bin/bash

cd madevent/Source/
make clean
make

cd ../SubProcesses
for dir in P1_*; do
  cd $dir
  make clean
  make
  cd -
done
