#!/bin/sh

if [ "$1" == "" ] || [ "$2" != "" ]; then
  echo "Usage:   $0 <proc.mad>"
  echo "Example: $0 pp_dy012j.mad"
  exit 1
fi

procdir=$1
if [ ! -d $procdir ]; then
  echo "ERROR! Directory not found $procdir"
  exit 1
fi

msgs=""
msgs="${msgs} PROGRAM_TOTAL"
msgs="${msgs} Fortran_Other"
msgs="${msgs} Fortran_Initialise(I/O)"
msgs="${msgs} Fortran_Random2Momenta"
msgs="${msgs} Fortran_PDFs"
msgs="${msgs} Fortran_UpdateScaleCouplings"
msgs="${msgs} Fortran_Reweight"
msgs="${msgs} Fortran_Unweight(LHE-I/O)"
msgs="${msgs} Fortran_SamplePutPoint"
msgs="${msgs} Fortran_MEs"
msgs="${msgs} CudaCpp_Initialise"
msgs="${msgs} CudaCpp_Finalise"
msgs="${msgs} CudaCpp_MEs"
msgs="${msgs} OVERALL_NON-MEs"
msgs="${msgs} OVERALL_MEs"

for backend in fortran cppnone cppsse4 cppavx2 cpp512y cpp512z cuda hip; do
  outfile=$procdir/${backend}/output.txt
  echo $outfile
  if [ ! -f $outfile ]; then
    echo "File not found: SKIP backend ${backend}"
  else
    cat $outfile | grep "__CUDACPP_DEBUG: GridPackCmd.launch finished" \
      | sed 's/__CUDACPP_DEBUG: GridPackCmd.launch finished in/[GridPackCmd.launch] GRIDPCK TOTAL   /'
    for msg0 in ${msgs}; do
      msg=${msg0/_/ }
      cat $outfile | grep "\[COUNTERS\]" | grep "${msg}" | sed 's/s for.*//' | sed 's/s$//' \
        | awk -vmsg="${msg}" -vttot=0 '{jtot=$NF; ttot += jtot}; END{if ( ttot!=0 ) print "[madevent COUNTERS] ", msg, ttot}'
    done
  fi
  echo "--------------------------------------------------------------------------------"
done
