#/bin/bash -f

for avx in 512z 512y avx2 sse4 none; do
  echo "CUDACPP_RUNTIME_GOODHELICITIES=ALL ./build.${avx}_m_inl0_hrd0/check_cpp.exe -p 1 16 1"
  CUDACPP_RUNTIME_GOODHELICITIES=ALL ./build.${avx}_m_inl0_hrd0/check_cpp.exe -p 1 16 1 | tee log_ggttggggg_m_${avx}.txt
done
