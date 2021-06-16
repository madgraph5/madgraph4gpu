#!/bin/bash

CUDAVERSION=11.1

export PATH=`echo $PATH | sed s/\\\\/cvmfs[a-zA-Z0-9./_-]\\\\+\://g`
export PATH=`echo $PATH | sed s/cuda-[0-9]*.[0-9]/cuda-${CUDAVERSION}/g`
#/eos/user/a/areepsch/Madgraph4gpu

