#!/bin/bash

# $1: *_dist_seqlevel.tsv file (for all wells) resulted by calc_distance.py
# $2: parental sequence tsv (sequence ID and its sequence)

dist_tsv_seqlevel=$1
parent_seq=$2

cat $dist_tsv_seqlevel | awk 'BEGIN{OFS="\t"; print "seq","parent_well","parent_seq","category"} {split($1,ar,"_");well=ar[4];cat="";if(well==$2 && $4==1){cat="PASS"}else if(well!=$2 && $4==1){cat="DROP_unexpected"}else if($4!=1){cat="DROP_redundant"}} {print $1,$2,$3,cat}' > detail.tsv

tail -n +2 detail.tsv | cut -f1,4 | sort | uniq > filter.summary

#category count
cut -f2 filter.summary | sort | uniq -c | awk 'BEGIN{print "cat\tcount"} {print $2"\t"$1}' > categorycount.tsv
redundant_cnt=$(grep DROP_redundant categorycount.tsv | cut -f 2)

#unexpected hit count
awk '$4=="DROP_unexpected"{print $3}' detail.tsv | sort | uniq -c | awk '{print $2"\t"$1}' > drop.unexpected.intermediate

#fill 0 for the no hit counts
awk 'FILENAME==ARGV[1]{ar[$1]=$1; print $0} FILENAME==ARGV[2]{if(ar[$1]==""){print $1"\t"0}}' drop.unexpected.intermediate $parent_seq > drop.unexpected

#unexpected well count
sort -k1 drop.unexpected | awk '{split($1,ar,"_"); if(NR==1){tmp=ar[1]; cnt=0}; if(tmp!=ar[1]){print tmp"\t"cnt; tmp=ar[1]; cnt=0; cnt+=$2}else{cnt+=$2}} END{print tmp"\t"cnt}' > drop.unexpected.wellcount

#passed well count
awk 'FILENAME==ARGV[1]{ar[$1]=$1; print $0} FILENAME==ARGV[2]{if(ar[$1]==""){print $1"\t"0}}' <(awk '$4=="PASS"{print $3}' detail.tsv | sort | uniq -c | awk '{print $2"\t"$1}') $parent_seq | sort -k1 | awk '{split($1,ar,"_"); if(NR==1){tmp=ar[1]; cnt=0}; if(tmp!=ar[1]){print tmp"\t"cnt; tmp=ar[1]; cnt=0; cnt+=$2}else{cnt+=$2}} END{print tmp"\t"cnt}' > pass.wellcount

#summarize
cat <(awk '{print $0"\tPASS"}' pass.wellcount) <(awk '{print $0"\tDrop (unexpected parental sequence)"}' drop.unexpected.wellcount) <(echo -e "NULL\t"$redundant_cnt"\tDrop (redundant best-hit parental sequence)") > contatmination.summary
