#/bin/bash -f

topdir=$(cd $(dirname $0)/../..; pwd)

START=$(date)
for dpg in dpg1000dpf1000 dpg200dpf200 dpg100dpf100 dpg10dpf100 dpg1dpf100; do
  cd ${topdir}/gg_ttgggg.${dpg}.sa/SubProcesses/P1_Sigma_sm_gg_ttxgggg; pwd
  for dcd in dcd0 dcd1; do
    exe=build.cuda_m_inl0_hrd0_${dcd}/check_cuda.exe
    if [ -f ${exe} ]; then
      for b in 1 2 4 8 16 32 64 128; do
        tmp=../../../tput/logs_ggttgggg_dpgs/logs_ggttgggg_${dpg}_sa_m_inl0_hrd0_${dcd}_${b}.tmp
        echo "CUDACPP_RUNTIME_GOODHELICITIES=ALL ${exe} -p $b 32 1 |& tee ${tmp}"
        CUDACPP_RUNTIME_GOODHELICITIES=ALL ${exe} -p $b 32 1 |& tee ${tmp}
      done
    fi
  done
done
echo $START
echo $(date)
