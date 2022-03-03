#!/bin/sh

table=
if [ "$1" == "-hrdcod" ]; then
  table="hrdcod"; shift
elif [ "$1" == "-juwels" ]; then
  table="juwels"; shift
fi
if [ "$1" != "" ]; then echo "Usage: $0 [-hrdcod|-juwels]"; exit 1; fi  

cd $(dirname $0)/..

# Output file
if [ "$table" == "" ]; then
  out=tput/summaryTable.txt
else
  out=tput/summaryTable_${table}.txt
fi
\rm -f $out
touch $out

# Select revisions
revs=""
if [ "$table" == "" ]; then
  ###revs="$revs c2e67b4" # cuda116/gcc102  (25 Jan 2022) BASELINE eemumu/ggtt/ggttgg x f/d x inl0/inl1 + ggttg/ggttggg x f/d x inl0
  ###revs="$revs 4f3229d" # cuda116/icx2021 (25 Jan 2022) ICX TEST eemumu/ggtt/ggttgg x f/d x inl0/inl1 + ggttg/ggttggg x f/d x inl0
  revs="$revs 09e482e" # cuda116/gcc102  (03 Mar 2022) BASELINE eemumu/ggtt/ggttgg x f/d x hrd0/hrd1 x inl0/inl1 + ggttg/ggttggg x f/d x hrd0/hrd1 x inl0
  revs="$revs 88dc717" # cuda116/icx2022 (03 Mar 2022) ICX TEST eemumu/ggtt/ggttgg x f/d x hrd0/hrd1 x inl0/inl1 + ggttg/ggttggg x f/d x hrd0/hrd1 x inl0
elif [ "$table" == "hrdcod" ]; then
  ###revs="$revs 2938acb" # cuda116/gcc102  (27 Jan 2022) HRD TEST eemumu/ggtt(g(g(g))) d inl0 x hrd0/hrd1
  revs="$revs 09e482e" # cuda116/gcc102  (03 Mar 2022) BASELINE eemumu/ggtt/ggttgg x f/d x hrd0/hrd1 x inl0/inl1 + ggttg/ggttggg x f/d x hrd0/hrd1 x inl0
  revs="$revs 88dc717" # cuda116/icx2022 (03 Mar 2022) ICX TEST eemumu/ggtt/ggttgg x f/d x hrd0/hrd1 x inl0/inl1 + ggttg/ggttggg x f/d x hrd0/hrd1 x inl0
elif [ "$table" == "juwels" ]; then
  revs="$revs c2e67b4" # cuda116/gcc102  (25 Jan 2022) BASELINE eemumu/ggtt/ggttgg x f/d x inl0/inl1 + ggttg/ggttggg x f/d
  revs="$revs 65730b2" # cuda115/gcc112  (18 Feb 2022) JUWELSCL eemumu/ggtt/ggttgg x f/d x inl0/inl1 + ggttg/ggttggg x f/d
  revs="$revs df441ad" # cuda115/gcc112  (18 Feb 2022) JUWELSBO eemumu/ggtt/ggttgg x f/d x inl0/inl1 + ggttg/ggttggg x f/d
fi

# Select processes
procs="eemumu ggtt ggttg ggttgg ggttggg"

# Select fptype, helinl, hrdcod
if [ "$table" == "" ]; then
  fpts="d f"
  inls="inl0 inl1"
  hrds="hrd0"
elif [ "$table" == "hrdcod" ]; then
  fpts="d"
  inls="inl0"
  hrds="hrd0 hrd1"
elif [ "$table" == "juwels" ]; then
  fpts="d f"
  inls="inl0 inl1"
  hrds="hrd0"
fi

# Iterate through log files
for fpt in $fpts; do
  echo -e "*** FPTYPE=$fpt ******************************************************************\n" >> $out
  for rev in $revs; do
    echo -e "+++ REVISION $rev +++" >> $out
    nodelast=
    for inl in $inls; do
      for hrd in $hrds; do
        files=""
        for proc in $procs; do
          file=tput/logs_${proc}_manu/log_${proc}_manu_${fpt}_${inl}_${hrd}.txt
          if [ -f $file ]; then files="$files $file"; fi
        done
        ###echo "*** FILES $files ***" >> $out
        if [ "$files" == "" ]; then continue; fi
        git checkout $rev $files >& /dev/null
        node=$(cat $files | grep ^On | sort -u)
	if [ "$nodelast" != "$node" ]; then echo -e "$node\n" >> $out; nodelast=$node; fi
        ###cat $files | awk '/^Process/{print $0}; /Workflow/{print $0}; /MECalcOnly/{print $0}'; exit 0
        cat $files | awk -vrev=$rev -vcomplast=none -vinllast=none -vhrdlast=none -vfptlast=none\
          '/^Process(.)/{split($0,a,"["); comp="["a[2]; if ( complast == "none" ){print comp; complast=comp}};\
           /^Process/{split($0,a,"]"); split(a[2],b,"="); inl=b[2]; if ( inl != inllast ){printf "HELINL="inl; inllast=inl}}\
           /^Process/{split($0,a,"]"); split(a[3],b,"="); hrd=b[2]; if ( hrd != hrdlast ){print " HRDCOD="hrd; hrdlast=hrd}}\
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
                for(iproc=1;iproc<=nproc;iproc++){proc=procs[iproc]; tput=tput_proc_tag[proc,tag]; if(tput==""){tput="--------"}; printf "%-12s", tput}; printf "\n"}}' >> $out
        echo "" >> $out
      done
    done
  done
done
git checkout HEAD tput/logs* >& /dev/null
cat $out
