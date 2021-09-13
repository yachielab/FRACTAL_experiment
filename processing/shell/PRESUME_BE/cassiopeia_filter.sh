#!/bin/bash


input_dir=$1
filekey=$2
cell_list=$3
outdir=$4
outname=$5

zcat $input_dir/*${filekey}*.MUTATIONTABLE.tsv.gz > ${outdir}/${outname}.cat.tmp

awk 'FILENAME==ARGV[1]{ar[$0]=$0} FILENAME==ARGV[2]{if(ar[$1]!=""){print}}' $cell_list ${outdir}/${outname}.cat.tmp | gzip -c > ${outdir}/${outname}.filtered.MUTATIONTABLE.tsv.gz
rm ${outdir}/${outname}.cat.tmp
