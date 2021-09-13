#!/bin/bash


input_dir=$1
chunk=$2
amp_len=$3
outdir=$4

for i in ${input_dir}/*chunk${chunk}.filtered.MUTATIONTABLE.tsv.gz; do zcat $i | awk -v L=$amp_len '{if($2==""){print $1"\tD_0_"L-1}else{print $0}}' | gzip -c > ${outdir}/$(echo ${i##*/} | cut -d "." -f1).filldel.MUTATIONTABLE.tsv.gz; done
