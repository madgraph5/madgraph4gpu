#!/bin/sh

cd $(dirname $0)/..

procs="eemumu ggtt ggttgg"

for proc in $procs; do

  file=tput/logs_${proc}_manu/log_${proc}_manu_d_inl0_hrd0.txt
  ###cat $file | awk '/Workflow/{print $0}; /MECalcOnly/{print $0}'
  cat $file | awk \
    '/^Process/{proc=""; split($3,a,"_"); proc=a[3]"_"a[4]};\
     /Workflow/{tag=""; split($4,a,":"); tag=a[1]; split($4,a,"+"); split(a[4],b,"/"); tag=tag"/"b[2]};\
     /MECalcOnly/{tput=sprintf("%.3g", $5); tputs[tag]=tput};\
     END{n=split("CUD/none CPP/none CPP/sse4 CPP/avx2 CPP/512y CPP/512z",tags); for(i=1;i<=n;i++){print proc, tags[i], tputs[tags[i]]}}'

done
