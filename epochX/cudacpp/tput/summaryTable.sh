#!/bin/sh

cd $(dirname $0)/..

out=tput/summaryTable.txt
\rm -f $out
touch $out

# Select revisions
revs=""
revs="$revs c2e67b4" # cuda116/gcc102 BASELINE (25 Jan 2022) eemumu/ggtt/ggttgg x f/d x inl0/inl1 + ggttg/ggttggg x f/d
revs="$revs 4f3229d" # cuda116/icx2021 (25 Jan 2022)  eemumu/ggtt/ggttgg x f/d x inl0/inl1 + ggttg/ggttggg x f/d

# Select processes
procs="eemumu ggtt ggttg ggttgg ggttggg"

# Select fptype
###fpts="d"
fpts="d f"

# Select helinl
###inls="inl0"
inls="inl0 inl1"

# Iterate through log files
for fpt in $fpts; do
  echo -e "*** FPTYPE=$fpt ******************************************************************\n" >> $out
  for rev in $revs; do
    for inl in $inls; do
      ###echo "*** REVISION $rev ***" >> $out
      files=""
      for proc in $procs; do
        file=tput/logs_${proc}_manu/log_${proc}_manu_${fpt}_${inl}_hrd0.txt
        if [ -f $file ]; then files="$files $file"; fi
      done
      ###echo "*** FILES $files ***" >> $out
      if [ "$files" == "" ]; then continue; fi
      git checkout $rev $files >& /dev/null
      ###cat $files | awk '/^Process/{print $0}; /Workflow/{print $0}; /MECalcOnly/{print $0}'; exit 0
      cat $files | awk -vrev=$rev -vcomplast=none -vinllast=none -vfptlast=none\
        '/^Process(.)*nvcc/{split($0,a,"["); comp="["a[2]; if ( comp != complast ){print "Revision", rev, comp; complast=comp}};\
         /^Process/{split($0,a,"]"); split(a[2],b,"="); inl=b[2]; if ( inl != inllast ){print "HELINL="inl; inllast=inl}}\
         /^Process/{proc=""; split($3,a,"_"); proc=a[3]"_"a[4]};\
         /^FP precision/{fpt=$4; /*if ( fpt != fptlast ){print "FPTYPE="fpt; fptlast=fpt}*/}\
         /Workflow/{tag=""; split($4,a,":"); tag=a[1]; split($4,a,"+"); split(a[4],b,"/"); tag=tag"/"b[2]};\
         /MECalcOnly/{tput=sprintf("%.2e", $5); tput_proc_tag[proc,tag]=tput}; \
         END{ntag=split("CUD/none CPP/none CPP/sse4 CPP/avx2 CPP/512y CPP/512z",tags);\
             nproc=split("EPEM_MUPMUM GG_TTX GG_TTXG GG_TTXGG GG_TTXGGG",procs);\
             procs_txt["EPEM_MUPMUM"]="eemumu";\
             procs_txt["GG_TTX"]="ggtt";\
             procs_txt["GG_TTXG"]="ggttg";\
             procs_txt["GG_TTXGG"]="ggttgg";\
             procs_txt["GG_TTXGGG"]="ggttggg";\
             printf "%12s", "";\
             for(iproc=1;iproc<=nproc;iproc++){proc=procs[iproc]; printf "%-12s", procs_txt[proc]}; printf "\n";\
             for(itag=1;itag<=ntag;itag++)\
             {tag=tags[itag]; printf "%-12s", tag;\
              for(iproc=1;iproc<=nproc;iproc++){proc=procs[iproc]; printf "%-12s", tput_proc_tag[proc,tag]}; printf "\n"}}' >> $out
      echo "" >> $out
    done
  done
done
git checkout HEAD tput/logs* >& /dev/null
cat $out
