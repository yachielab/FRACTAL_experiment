#!/bin/bash
#$ -S /bin/bash

part=$1
target_dir=$2
maxchunk=$3
outdir=$4
outname=$5

starting_file=$(ls ${target_dir}/* | grep ${part}_delsub.SEQUENCE.tsv.gz | head -n1)
command="paste <(gzip -cd $starting_file | cut -f1)"
for i in $(seq 0 $(expr $maxchunk - 1)); do
    command="$command <(gzip -cd ${target_dir}/*.${i}.${part}_delsub.SEQUENCE.tsv.gz | cut -f2)"
done

eval $command | \
awk 'BEGIN{OFS="\t"} {col2=""; for(i=1;i<=NF;i++){if(i==1){col1=$1}else{col2=col2""$i}}; print col1,col2}' | \
seqkit tab2fx -o ${outdir}/${outname}.fa.gz
