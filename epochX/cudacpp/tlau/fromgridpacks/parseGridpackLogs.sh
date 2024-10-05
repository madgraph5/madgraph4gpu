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
msgs="${msgs} Fortran_PhaseSpaceSampling"
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

teefile=$procdir/summary.txt
\rm -f $teefile
touch $teefile
for backend in fortran cppnone cppsse4 cppavx2 cpp512y cpp512z cuda hip; do
  outfile=$procdir/${backend}/output.txt
  echo $outfile | tee -a $teefile
  if [ ! -f $outfile ]; then
    echo "File not found: SKIP backend ${backend}" | tee -a $teefile
  else
    cat $outfile | grep "__CUDACPP_DEBUG: GridPackCmd.launch finished" \
      | sed 's/__CUDACPP_DEBUG: GridPackCmd.launch finished in/[GridPackCmd.launch] GRIDPCK TOTAL                 /' \
      | awk '{printf "%s %-30s %10.4f\n", $1, $2" "$3, $4}' \
      | tee -a $teefile
    for msg0 in ${msgs}; do
      msg=${msg0/_/ }
      cat $outfile | grep "\[COUNTERS\]" \
        | sed -r 's/ (PROGRAM|Fortran|CudaCpp|OVERALL) / \1_/' \
        | grep "${msg0}" | sed -r 's/(\( *[0-9]* *\))//' \
        | sed 's/s for/ for/' | sed 's/s$//' \
        | awk -vmsg="${msg}" -vstot=0 -vetot=0 \
              '{jstot=$4; jetot=$6; stot += jstot; etot += jetot}; \
               END{if ( stot!=0 ) { printf "[madevent COUNTERS]  %-30s %10.4fs", msg, stot;\
                   if ( etot!=0 ) printf " for %10d events\n", etot; else printf "\n" } }' \
	| tee -a $teefile
    done
  fi
  echo "------------------------------------------------------------------------------------------" | tee -a $teefile
done
