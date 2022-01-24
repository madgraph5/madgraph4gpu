#!/bin/sh

cd $(dirname $0)/..

procs="eemumu ggtt ggttgg"
files=""
for proc in $procs; do
  file=tput/logs_${proc}_manu/log_${proc}_manu_d_inl0_hrd0.txt
  files="$files $file"
done

###cat $files | awk '/^Process/{print $0}; /Workflow/{print $0}; /MECalcOnly/{print $0}'
cat $files | awk \
    '/^Process/{proc=""; split($3,a,"_"); proc=a[3]"_"a[4]};\
     /Workflow/{tag=""; split($4,a,":"); tag=a[1]; split($4,a,"+"); split(a[4],b,"/"); tag=tag"/"b[2]};\
     /MECalcOnly/{tput=sprintf("%.3g", $5); tput_proc_tag[proc,tag]=tput}; \
     END{ntag=split("CUD/none CPP/none CPP/sse4 CPP/avx2 CPP/512y CPP/512z",tags);\
         nproc=split("EPEM_MUPMUM GG_TTX GG_TTXGG",procs);\
         for(iproc=1;iproc<=nproc;iproc++)\
         {proc=procs[iproc]; for(itag=1;itag<=ntag;itag++){tag=tags[itag]; print proc, tag, tput_proc_tag[proc,tag]}}}'
