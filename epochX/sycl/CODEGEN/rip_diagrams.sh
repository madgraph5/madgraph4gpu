#!/bin/bash
process_file="CPPProcess.cc"
process_out="diagrams.h"
num_diagrams=123

rm -f ${process_out}
echo "#ifndef DIAGRAMS_H " >> ${process_out}
echo "#define DIAGRAMS_H" >> ${process_out}
echo "" >> ${process_out}
echo "#include \"mgOnGpuConfig.h\"" >> ${process_out}
echo "#include \"mgOnGpuTypes.h\"" >> ${process_out}
echo "#include \"mgOnGpuVectors.h\"" >> ${process_out}
echo "" >> ${process_out}
echo "#include \"HelAmps_sm.h\"" >> ${process_out}
echo "" >> ${process_out}
echo "using namespace MG5_sm;" >> ${process_out}
echo "static constexpr size_t np4 = mgOnGpu::np4;" >> ${process_out}
echo "static constexpr size_t neppM = mgOnGpu::neppM;" >> ${process_out}
echo "" >> ${process_out}

for i in $(seq 1 $(echo "${num_diagrams} - 1" | bc))
do
    j=$(echo "${i} + 1" | bc)
    line_start=$(cat ${process_file} | awk "/\<DIAGRAM ${i}\>/{print NR}")
    line_end=$(cat ${process_file} | awk "/\<DIAGRAM ${j}\>/{print NR}")
    line_end=$(echo "${line_end} - 1" | bc)
    diagram=$(cat ${process_file} | awk "/\<DIAGRAM ${i}\>/")
    echo "SYCL_EXTERNAL" >> ${process_out} 
    echo "void diagram_call_${i}_of_${num_diagrams}(" >> ${process_out}
    echo "       const fptype_sv* __restrict__ allmomenta," >> ${process_out}
    echo "       #ifdef MGONGPU_SUPPORTS_MULTICHANNEL" >> ${process_out}
    echo "           fptype* __restrict__ allNumerators," >> ${process_out}
    echo "           fptype* __restrict__ allDenominators," >> ${process_out}
    echo "           const size_t channelId," >> ${process_out}
    echo "       #endif" >> ${process_out}
    echo "       const short*  __restrict__ cHel," >> ${process_out}
    echo "       const cxtype* __restrict__ COUPs," >> ${process_out}
    echo "       const fptype* __restrict__ cIPD," >> ${process_out}
    echo "       cxtype_sv* __restrict__ w_sv," >> ${process_out}
    echo "       cxtype_sv* __restrict__ amp_sv," >> ${process_out}
    echo "       cxtype_sv* __restrict__ jamp_sv" >> ${process_out}
    echo "       ) {" >> ${process_out}
    echo "${diagram}" >> ${process_out} 
    tac ${process_file} | awk "/\<DIAGRAM ${i}\>/{p=0;exit}p;/\<DIAGRAM ${j}\>/{p=1}" | tac >> ${process_out}
    echo "}" >> ${process_out}
    echo "" >> ${process_out}
    diagram_call="    diagram_call_${i}_of_${num_diagrams}(\n            allmomenta,\n            #ifdef MGONGPU_SUPPORTS_MULTICHANNEL\n                allNumerators,\n                allDenominators,\n            channelId,\n            #endif\n            cHel,\n            COUPs,\n            cIPD,\n            w_sv,\n            amp_sv,\n            jamp_sv\n        );\n"
    sed -i "${line_start},${line_end}d" ${process_file}
    sed -i "${line_start}s/^/${diagram_call}/g" ${process_file}
done

line_start=$(cat ${process_file} | awk "/\<DIAGRAM ${num_diagrams}\>/{print NR}")
line_end=$(cat ${process_file} | awk "/\<COLOR ALGEBRA\>/{print NR}")
line_end=$(echo "${line_end} - 1" | bc)
diagram=$(cat ${process_file} | awk "/\<DIAGRAM ${num_diagrams}\>/")

echo "SYCL_EXTERNAL" >> ${process_out} 
echo "void diagram_call_${num_diagrams}_of_${num_diagrams}(" >> ${process_out}
echo "       const fptype* __restrict__ allmomenta," >> ${process_out}
echo "       #ifdef MGONGPU_SUPPORTS_MULTICHANNEL" >> ${process_out}
echo "           fptype* __restrict__ allNumerators," >> ${process_out}
echo "           fptype* __restrict__ allDenominators," >> ${process_out}
echo "           const size_t channelId," >> ${process_out}
echo "       #endif" >> ${process_out}
echo "       const short*  __restrict__ cHel," >> ${process_out}
echo "       const cxtype* __restrict__ COUPs," >> ${process_out}
echo "       const fptype* __restrict__ cIPD," >> ${process_out}
echo "       cxtype_sv* __restrict__ w_sv," >> ${process_out}
echo "       cxtype_sv* __restrict__ amp_sv," >> ${process_out}
echo "       cxtype_sv* __restrict__ jamp_sv" >> ${process_out}
echo "       ) {" >> ${process_out}
echo "${diagram}" >> ${process_out} 
tac ${process_file} | awk "/\<DIAGRAM ${num_diagrams}\>/{p=0;exit}p;/\<COLOR ALGEBRA\>/{p=1}" >> ${process_out}
echo "}" >> ${process_out}
echo "" >> ${process_out}
echo "#endif" >> ${process_out}

diagram_call="    diagram_call_${num_diagrams}_of_${num_diagrams}(\n            allmomenta,\n            #ifdef MGONGPU_SUPPORTS_MULTICHANNEL\n                allNumerators,\n                allDenominators,\n            channelId,\n            #endif\n            cHel,\n            COUPs,\n            cIPD,\n            w_sv,\n            amp_sv,\n            jamp_sv\n        );\n"
sed -i "${line_start},${line_end}d" ${process_file}
sed -i "${line_start}s/^/${diagram_call}/g" ${process_file}

