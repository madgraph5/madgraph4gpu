#!/bin/bash

#----------------------------------------------------------------------------------------------

# Auxiliary function of mainSummarizeSyms
# Env: $dumptmp must point to the stripped objdump file
function stripSyms() {
  if [ "$dumptmp" == "" ] || [ ! -e $dumptmp ]; then echo "ERROR! File '$dumptmp' not found"; exit 1; fi 
  cnt=0
  for sym in "$@"; do
    # Use cut -f3- to print only the assembly code after two leading fields separated by tabs
    (( cnt += $(cat $dumptmp | egrep "$sym" | wc -l) ))
    cat $dumptmp | egrep -v "$sym" > ${dumptmp}.new
    \mv ${dumptmp}.new $dumptmp
  done; echo $cnt
}

#----------------------------------------------------------------------------------------------

function mainSummarizeSyms() {

  # Command line arguments: HelAmps functions only?
  helamps=0
  if [ "$1" == "-helamps" ]; then
    helamps=1
    shift
  fi

  # Command line arguments: strip dir when showing the file name?
  stripdir=0
  if [ "$1" == "-stripdir" ]; then
    stripdir=1
    shift
  fi

  # Command line arguments: select file
  if [ "$1" == "" ] || [ "$2" != "" ]; then
    echo "Usage:   $0 [-helamps] [-stripdir] <filename>"
    echo "Example: $0 ./check.exe"
    exit 1
  fi
  file=$1

  # Disassemble selected file
  # Use cut -f3- to print only the assembly code after two leading fields separated by tabs
  objdump -d -C $file > ${file}.objdump # unnecessary but useful for debugging
  dumptmp=${file}.objdump.tmp
  if [ "$helamps" == "0" ]; then
    objdump -d -C $file | awk '/^ +[[:xdigit:]]+:\t/' | cut -f3- > ${dumptmp}
  else
    objdump -d -C $file | awk -v RS= '/^[[:xdigit:]]+ <MG5_sm::.*/' | awk '/^ +[[:xdigit:]]+:\t/' | cut -f3- > ${dumptmp}
  fi
  ###ls -l $dumptmp

  unamep=$(uname -p)
  #--- ARM ---
  if [ "${unamep}" == "arm" ]; then 

    # FIXME: classifying objdump symbols for ARM has not been done yet
    return

  #--- PPC ---
  # See https://cdn.openpowerfoundation.org/wp-content/uploads/resources/Intrinsics-Reference_final/Intrinsics-Reference-20200811.pdf
  elif [ "${unamep}" == "ppc64le" ]; then 

    # Exclude all instructions not involving "vs" registers
    cat $dumptmp | grep " vs" > ${dumptmp}.new
    \mv ${dumptmp}.new $dumptmp

    # Count and strip "xv" symbols
    cntXV=$(stripSyms '^xv.* vs' '^lxv.* vs' '^stxv.* vs')
    ###echo $cntXV; ls -l $dumptmp

    # Count and strip "xx" symbols
    cntXX=$(stripSyms '^xx.* vs')
    ###echo $cntXX; ls -l $dumptmp

    # Count and strip "xs" symbols
    cntXS=$(stripSyms '^xs.* vs')
    ###echo $cntXS; ls -l $dumptmp

    # Count and strip "mtv" symbols
    cntMTV=$(stripSyms '^mtv.* vs')
    ###echo $cntMTV; ls -l $dumptmp

    # Is there anything else?...
    cntelse=$(stripSyms '.* vs')
    ###echo $cntelse; ls -l $dumptmp
    if [ "$cntelse" != "0" ]; then echo "ERROR! cntelse='$cntelse'"; exit 1; fi 

    # Final report
    filename=$file
    if [ "$stripdir" == "1" ]; then filename=$(basename $file); fi
    printf "=Symbols (vs) in $filename= (^mtv:%5d) (^xs:%5d) (^xx:%5d) (^xv:%5d)\n" $cntMTV $cntXS $cntXX $cntXV

  #--- x86 ---
  else

    # Count and strip AVX512 zmm symbols
    cnt512z=$(stripSyms '^v.*zmm')
    ###echo $cnt512z; ls -l $dumptmp

    # Count and strip AVX512 ymm/xmm symbols
    # [NB: these are AVX512VL symbols, i.e. implementing AVX512 on xmm/ymm registers]
    cnt512y=$(stripSyms '^v.*dqa(32|64).*(x|y)mm' '^v.*(32|64)x2.*(x|y)mm' '^vpcmpneqq.*(x|y)mm' '^vpermi2.*(x|y)mm' '^vblendm.*(x|y)mm' '^vpmovd.*(x|y)mm' '^vrnd.*(x|y)mm')
    ###echo $cnt512y; ls -l $dumptmp

    # Count and strip AVX2 symbols
    cntavx2=$(stripSyms '^v.*(x|y)mm')
    ###echo $cntavx2; ls -l $dumptmp

    # Count and strip ~SSE4 symbols
    ###cntsse4=$(stripSyms '^[^v].*xmm') # 'too many'(*): includes many symbols from "none" build
    cntsse4=$(stripSyms '^[^v].*p(d|s).*xmm') # okish, representative enough of what sse4/nehalem adds
    ###cntsse4=$(stripSyms '^(mul|add|sub)p(d|s).*xmm') # too few: should also include moves
    ###echo $cntsse4; ls -l $dumptmp

    # Is there anything else?...
    ###cntelse=$(stripSyms '.*(x|y|z)mm') # only makes sense when counting 'too many'(*) as sse4
    ###echo $cntelse; ls -l $dumptmp
    ###if [ "$cntelse" != "0" ]; then echo "ERROR! cntelse='$cntelse'"; exit 1; fi 

    # Final report
    filename=$file
    if [ "$stripdir" == "1" ]; then filename=$(basename $file); fi
    printf "=Symbols in $filename= (~sse4:%5d) (avx2:%5d) (512y:%5d) (512z:%5d)\n" $cntsse4 $cntavx2 $cnt512y $cnt512z

  fi
}

mainSummarizeSyms $*

#----------------------------------------------------------------------------------------------
