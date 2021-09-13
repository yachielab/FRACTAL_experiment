#!/bin/bash

in_file=$1
outdir=$2
outname=$3

cat $in_file | \
awk 'BEGIN{OFS=";"} NR>1{print $1"\t"$2,$3,$4}' | \
awk 'NR==1{cell=$1; mut=$2} NR>1&&cell==$1{mut=mut";"$2} NR>1&&cell!=$1{print cell"\t"mut; cell=$1; mut=$2} END{print cell"\t"mut}' | \
awk '{mut=""; prev=""; split($2,ar,";"); for(i=1;i<=length(ar);i++){if(ar[i]!="None"&&ar[i]!=prev){mut=mut";"ar[i]; prev=ar[i]}}; sub("^;","",mut); print $1"\t"mut}' | \
gzip -c > ${outdir}/${outname}.back.MUTATIONTABLE.tsv.gz
