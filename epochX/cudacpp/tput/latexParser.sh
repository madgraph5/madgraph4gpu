#!/bin/sh

cd $(dirname $0)/..

file=tput/logs_eemumu_manu/log_eemumu_manu_d_inl0_hrd0.txt
###cat $file | awk '/Workflow/{print $0}; /MECalcOnly/{print $0}'
cat $file | awk \
  '/Workflow/{tag=""; split($4,a,":"); tag=a[1]; split($4,a,"+"); split(a[4],b,"/"); tag=tag"/"b[2]};\
   /MECalcOnly/{tput=sprintf("%.3g", $5); tputs[tag]=tput};\
   END{n=split("CUD/none CPP/none CPP/sse4 CPP/avx2 CPP/512y CPP/512z",tags); for(i=1;i<=n;i++){print tags[i], tputs[tags[i]]}}'
