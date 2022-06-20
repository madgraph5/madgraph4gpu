#!/bin/sh

# Kernel function
function oneTable()
{
  taglist="FORTRAN CPP/none CPP/sse4 CPP/avx2 CPP/512y CPP/512z CUDA/32 CUDA/8192"
  parlist="(1) (2-none) (2-sse4) (2-avx2) (2-512y) (2-512z) (3) (3bis)"
  faclist="1 10 100"
  echo "" > $out
  for proc in $procs; do
    file=tmad/logs_${proc}_${suff}/log_${proc}_${suff}_${fpt}_${inl}_${hrd}.txt
    if [ ! -f $file ]; then continue; fi
    ###echo "*** FILE $file ***"
    cat $file | awk -vproc=$proc -vtaglist="$taglist" -vparlist="$parlist" -vfaclist="$faclist" '
      BEGIN{ntag=split(taglist,tags); npar=split(parlist,pars);
            if(ntag!=npar){print "ERROR! ntag!=npar", ntag, npar; exit(1)};
            for(i=1;i<=ntag;i++){tag1[pars[i]]=tags[i]}}
      BEGIN{nfac=split(faclist,facs)}
      BEGIN{lsep=sprintf("%094d",0); gsub("0","-",lsep)}
      ###/create events.lhe/{print $0}
      /create events.lhe/{par=$2; tag=tag1[par]} # current tag
      /create events.lhe/{fac=substr($5,2)} # current fac
      /\[COUNTERS\]/{if($3=="TOTAL") typ=1; else if($3=="Overhead") typ=2; else if($3=="MEs") typ=3; else{print "ERROR! Unknown type $3"; exit 1};
                     if($4==":") sec=$5; else sec=$8; sec=substr(sec,1,length(sec)-1); # current sec
                     sec3[tag,fac,typ]=sec}
      END{print lsep;
          printf "| %-9s | %-78s |\n", proc, "[sec] Total = Overhead (FORTRAN) + MEs (FORTRAN, CPP or CUDA)";
          print lsep;
          for (itag=1; itag<=ntag; itag++)
          {tag=tags[itag]; printf "| %-9s |", tag;
           for(ifac=1; ifac<=nfac; ifac++)
           {fac=facs[ifac]; printf " %6.2f = %6.2f + %6.2f |", sec3[tag,fac,1], sec3[tag,fac,2], sec3[tag,fac,3]};
           printf "\n"};
          print lsep;
          print "";
         }' >> $out
  done
}

cd $(dirname $0)/..
echo PWD=$(pwd)

procs="eemumu ggtt ggttg ggttgg ggttggg"
#procs="ggttggg"
suff=mad
fpt=d
inl=inl0
hrd=hrd0

out=tmad/summaryTable_default.txt
oneTable
cat $out
