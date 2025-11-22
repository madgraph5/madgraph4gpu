#/bin/bash -f

topdir=$(cd $(dirname $0)/../..; pwd)

outmref=
for dpg in dpg1000dpf1000 dpg200dpf200 dpg100dpf100 dpg10dpf100 dpg1dpf100; do
  cd ${topdir}/gg_ttgggg.${dpg}.sa/SubProcesses/P1_Sigma_sm_gg_ttxgggg; pwd
  for dcd in dcd0 dcd1; do
    exe=build.cuda_m_inl0_hrd0_${dcd}/check_cuda.exe
    if [ -f ${exe} ]; then
      outs=../../../tput/logs_ggttgggg_dpgs/logs_ggttgggg_${dpg}_sa_m_inl0_hrd0_${dcd}.scaling
      outm=../../../tput/logs_ggttgggg_dpgs/logs_ggttgggg_${dpg}_sa_m_inl0_hrd0_${dcd}.mes
      \rm -f $outs; touch $outs
      \rm -f $outm; touch $outm
      for b in 1 2 4 8 16 32 64 128; do
        tmp=../../../tput/logs_ggttgggg_dpgs/logs_ggttgggg_${dpg}_sa_m_inl0_hrd0_${dcd}_${b}.tmp
        ###echo "CUDACPP_RUNTIME_GOODHELICITIES=ALL ${exe} -p $b 32 1 |& tee ${tmp}"
        ###CUDACPP_RUNTIME_GOODHELICITIES=ALL ${exe} -p $b 32 1 |& tee ${tmp}
        cat $tmp | awk -vb=$b '/assertGpu/{printf "%s %4d %3d\n", "check_cuda.exe: Assertion failed", b, 32}' >> $outs
        cat $tmp | awk -vb=$b '/EvtsPerSec\[MECalcOnly\]/{printf "%s %4d %3d\n", $5, b, 32}' >> $outs
        cat $tmp | awk -vb=$b '/assertGpu/{printf "%s %4d %3d\n", "check_cuda.exe: Assertion failed", b, 32}' >> $outm
        cat $tmp | awk -vb=$b '/MeanMatrixElemValue/{printf "%s %4d %3d\n", $4, b, 32}' >> $outm
      done
      if [ ! -s $outs ]; then \rm -f $outs; fi # empty file (all tests failed)
      if [ ! -s $outm ]; then \rm -f $outm; continue; fi # empty file (all tests failed)
      if [ "$outmref" == "" ]; then outmref=$outm; fi
      echo "diff $outm $outmref"
      diff $outm $outmref
    fi
  done
done
