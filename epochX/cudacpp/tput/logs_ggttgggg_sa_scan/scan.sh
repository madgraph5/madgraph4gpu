#/bin/bash -f

topdir=$(cd $(dirname $0)/../..; pwd)

START=$(date)
for dpg in dpg1dpf100 dpg10dpf100 dpg100dpf100 dpg200dpf200 dpg1000dpf1000 dpg10000dpf10000; do
  cd ${topdir}/gg_ttgggg.${dpg}.sa/SubProcesses/P1_Sigma_sm_gg_ttxgggg; pwd
  make -j -f cudacpp.mk bldall
  dir=../../../tput/logs_ggttgggg_sa_scan/${dpg}; mkdir -p ${dir}
  for avx in 512z 512y avx2 sse4 none; do
    echo "CUDACPP_RUNTIME_GOODHELICITIES=ALL ./build.${avx}_m_inl0_hrd0/check_cpp.exe -p 1 16 1"
    CUDACPP_RUNTIME_GOODHELICITIES=ALL ./build.${avx}_m_inl0_hrd0/check_cpp.exe -p 1 16 1 | tee ${dir}/log_ggttgggg_m_${avx}.txt
  done
done
echo $START
echo $(date)
