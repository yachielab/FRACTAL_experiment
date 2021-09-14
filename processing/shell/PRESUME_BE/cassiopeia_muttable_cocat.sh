#!/bin/bash
#$ -S /bin/bash

set -eu

prefix=$1
target_dir=$2
maxchunk=$3
amp_length=$4
outdir=$5

starting_file=$(ls ${target_dir}/*.tsv.gz | grep ${prefix}_chunk1.filldel.MUTATIONTABLE.tsv.gz | head -n1)
command="paste <(gzip -cd $starting_file | cut -f1)"

for i in $(seq 1 $maxchunk); do
    command="$command <(gzip -cd ${target_dir}/${prefix}_chunk${i}.filldel.MUTATIONTABLE.tsv.gz | cut -f2)"
done

echo $command | tr " " "\n" | grep tsv.gz > ${outdir}/${prefix}_${maxchunk}chunks.usedFiles.log.txt

eval $command | \
awk \
-F "\t" \
-v amp_length=$amp_length \
'BEGIN{OFS="\t";PROCINFO["sorted_in"] = "@ind_num_asc"} \
{col1=$1; \
col2=$2; \
for(i=3;i<=NF;i++){\
    chunk=i-2;\
    split($i,ar,";");\
    for(x in ar){\
        split(ar[x],loc,"_");\
        if(loc[1]=="I" || loc[1]=="S"){\
            loc[2]=chunk*amp_length+loc[2]\
        };\
        if(loc[1]=="D"){\
            loc[2]=chunk*amp_length+loc[2];\
            loc[3]=chunk*amp_length+loc[3];\
        };\
        col2=col2";"loc[1]"_"loc[2]"_"loc[3]
    };\
}; \
print col1,col2}' | \
sed -e 's/\t;/\t/g' | \
gzip -c > ${outdir}/${prefix}_${maxchunk}chunks.MUTATIONTABLE.tsv.gz
