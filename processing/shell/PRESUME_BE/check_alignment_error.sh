#!/bin/bash

. ~/.bashrc
. ~/.bash_profile
getpath

in_seq=$1
in_tbl=$2
out=$3
original_fa=$4

python /path/to/processing/py/PRESUME_BE/seq_reverse_fasta.py $in_seq $in_tbl ${out}.filecheck.tsv
seqkit tab2fx ${out}.filecheck.tsv -o ${out}.filecheck.fa

cut -f1 ${out}.filecheck.tsv > ${out}.celllist
rm ${out}.filecheck.tsv

awk 'FILENAME==ARGV[1]{ar[$0]=$0} FILENAME==ARGV[2]{if(ar[$1]!=""){print $0}}' ${out}.celllist <(seqkit fx2tab $original_fa) | seqkit tab2fx -o ${out}.original.fa
rm ${out}.celllist
