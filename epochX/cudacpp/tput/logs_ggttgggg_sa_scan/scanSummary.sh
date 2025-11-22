#/bin/bash -f

cd $(dirname $0)

me0=""
for dpg in dpg1dpf100 dpg10dpf100 dpg100dpf100 dpg200dpf200 dpg1000dpf1000 dpg10000dpf10000; do
  echo $dpg
  for avx in none sse4 avx2 512y 512z; do
    out=${dpg}/log_ggttgggg_m_${avx}.txt
    time=$(cat $out | awk '/SigmaKin/{print $4}')
    tput=$(cat $out | awk '/EvtsPerSec\[MECalcOnly\]/{print $5}')
    if [ "${avx}" == "none" ]; then tput0=$tput; fi
    spup=$(python3 -c "spup=$tput/$tput0; print('%5.2fx'%spup)")
    me=$(cat $out | awk '/MeanMatrixElemValue/{print $4}')
    if [ "${me0}" == "" ]; then me0=$me; fi
    dme=$(python3 -c "print($me-$me0)")
    echo $avx $time $dme $tput $spup
  done
  echo
done
