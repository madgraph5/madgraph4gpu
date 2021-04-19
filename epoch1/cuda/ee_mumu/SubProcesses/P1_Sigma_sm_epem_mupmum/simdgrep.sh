#!/bin/bash

#----------------------------------------------------------------------------------------------
# See http://sponce.web.cern.ch/sponce/CSC/slides/PracticalVectorization.booklet.pdf
# See https://software.intel.com/sites/landingpage/IntrinsicsGuide
# [NB: the trailing 'd' is for double; a trailing 's' would be for single precision]
#----------------------------------------------------------------------------------------------
# https://software.intel.com/sites/landingpage/IntrinsicsGuide/#text=addsd&expand=154
# (x1) : addsd xmm
# https://software.intel.com/sites/landingpage/IntrinsicsGuide/#text=addpd&expand=121,124,127
# (x2) : addpd xmm
# (x4) : vaddpd ymm
# (x8) : vaddpd zmm
#----------------------------------------------------------------------------------------------
# https://software.intel.com/sites/landingpage/IntrinsicsGuide/#text=mulsd&expand=3949
# (x1) : mulsd xmm
# https://software.intel.com/sites/landingpage/IntrinsicsGuide/#text=mulpd&expand=3919,3922,3925
# (x2) : mulpd xmm
# (x4) : vmulpd ymm
# (x8) : vmulpd zmm
#----------------------------------------------------------------------------------------------
# Below NNN is 132, 213, 231
# https://software.intel.com/sites/landingpage/IntrinsicsGuide/#text=vfmadd&expand=2537,2541,2545
# (x2) : vfmaddNNNpd xmm
# (x4) : vfmaddNNNpd ymm
# (x8) : vfmaddNNNpd zmm
#----------------------------------------------------------------------------------------------

# Auxiliary function of mainCountSyms
# Env: $dump must point to the objdump file
function countSyms() {
  if [ "$dump" == "" ] || [ ! -e $dump ]; then echo "ERROR! File '$dump' not found"; exit 1; fi 
  for sym in $*; do
    printf "(%16s : %3d) " $sym $(cat $dump | egrep $sym | wc -l)
  done; printf "\n"
}

###dump=$1 && countSyms # for debugging

#----------------------------------------------------------------------------------------------

function mainCountSyms() {
  # Command line arguments: select file
  if [ "$1" == "" ] || [ "$2" != "" ]; then
    echo "Usage:   $0 <filename>"
    echo "Example: $0 ./check.exe"
    exit 1
  fi
  file=$1
  # Disassemble selected file
  dump=${file}.objdump
  objdump -d -C $file > ${dump}
  # Count symbols in selected file
  single=0
  double=1
  if [ "$single" == "1" ]; then
    countSyms 'addss.*xmm' 'addps.*xmm' 'vaddps.*ymm' 'vaddps.*zmm'
    countSyms 'mulss.*xmm' 'mulps.*xmm' 'vmulps.*ymm' 'vmulps.*zmm'
    countSyms '==========' 'vfmadd...ps.*xmm' 'vfmadd...ps.*ymm' 'vfmadd...ps.*zmm'
  fi
  if [ "$double" == "1" ]; then
    countSyms 'addsd.*xmm' 'addpd.*xmm' 'vaddpd.*ymm' 'vaddpd.*zmm'
    countSyms 'mulsd.*xmm' 'mulpd.*xmm' 'vmulpd.*ymm' 'vmulpd.*zmm'
    countSyms 'subsd.*xmm' 'subpd.*xmm' 'vsubpd.*ymm' 'vsubpd.*zmm'
    countSyms 'movasd.*xmm' 'movapd.*xmm' 'vmovapd.*ymm' 'vmovapd.*zmm'
    countSyms '==========' 'vfmadd...pd.*xmm' 'vfmadd...pd.*ymm' 'vfmadd...pd.*zmm'
  fi
}

#mainCountSyms $*

#----------------------------------------------------------------------------------------------

# Auxiliary function of mainListSyms
# Env: $dump must point to the objdump file
function listSyms() {
  if [ "$dump" == "" ] || [ ! -e $dump ]; then echo "ERROR! File '$dump' not found"; exit 1; fi 
  # Use cut -f3- to print only the assembly code after two leading fields separated by tabs
  cat $dump | awk '/^ +[[:xdigit:]]+:\t/' | cut -f3- | egrep '(x|y|z)mm' | sed -r 's/ .*%(x|y|z)mm.*/ %\1mm/g' | sort | uniq -c | awk '{printf "%15s %4s %5d\n",$2,$3,$1}' | sort
}

###dump=$1 && listSyms # for debugging

#----------------------------------------------------------------------------------------------

function mainListSyms() {
  # Command line arguments: select file
  if [ "$1" == "" ] || [ "$2" != "" ]; then
    echo "Usage:   $0 <filename>"
    echo "Example: $0 ./check.exe"
    exit 1
  fi
  file=$1
  # Disassemble selected file
  dump=${file}.objdump
  objdump -d -C $file > ${dump}
  # List symbols in selected file
  listSyms
}

#mainListSyms $*

#----------------------------------------------------------------------------------------------

function mainCompareSyms() {
  allFileSyms=""
  for avx in none sse4 avx2 512y 512z; do
    file=./build.$avx/CPPProcess.o
    fileSymsRaw=$file.symlist.raw
    mainListSyms $file > $fileSymsRaw
    ls -l $fileSymsRaw
    allFileSyms="$allFileSyms $fileSymsRaw"
  done
  allSymList=all.symlist
  cat $allFileSyms | awk '{print $1, $2}' | sort -u > $allSymList
  for avx in none sse4 avx2 512y 512z; do
    file=./build.$avx/CPPProcess.o
    fileSymsRaw=$file.symlist.raw
    fileSyms=$file.symlist
    cat $fileSymsRaw | awk -vall=$allSymList -vfmt="%15s %4s %5d\n" -vtot=0 '{f1=$1; f2=$2; f3=$3; tot+=f3; while(f3>0){getline < all; if($1==f1 && $2==f2){printf fmt,$1,$2,f3; f3=0} else {printf fmt,$1,$2,0}}} END{printf fmt,"TOTAL","",tot}' > $fileSyms
    ls -l $fileSyms
  done
}

mainCompareSyms

#----------------------------------------------------------------------------------------------

