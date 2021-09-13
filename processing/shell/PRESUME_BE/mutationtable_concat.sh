#!/bin/bash
#$ -S /bin/bash

part=$1
target_dir=$2
amp_length=$3
maxchunk=$4
outdir=$5
outname=$6

starting_file=$(ls ${target_dir}/* | grep ${part}.MUTATIONTABLE.tsv.gz | head -n1)
command="paste <(gzip -cd $starting_file | cut -f1)"
for i in $(seq 0 $(expr $maxchunk - 1)); do
    command="$command <(gzip -cd ${target_dir}/*.${i}.${part}.MUTATIONTABLE.tsv.gz | cut -f2)"
done

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
gzip -c > ${outdir}/${outname}.tsv.gz
