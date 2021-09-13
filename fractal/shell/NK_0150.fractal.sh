#!/bin/bash
#$ -S /bin/bash
#$ -N NK_0150
#$ -cwd
#$ -o /home/naoki-konno/data/NK_0150/out
#$ -e /home/naoki-konno/data/NK_0150/err

CWD="/home/naoki-konno/data"

i=1
for placement in ML MP; do
    for ali in unaligned aligned; do 
        for target in bcat bet002; do
            for method in rapidnjNJ raxmlMP fasttreeML; do
                if [ $i -eq ${SGE_TASK_ID} ]; then


                    INPUT="/home/naoki-konno/data/NK_0150/fractalin/${ali}/fractalin.${target}.fa.gz"
                    cd $CWD/NK_0150/fractal/${target}
                    mkdir TASK${i}
                    cd    TASK${i}
                    echo "TASK,${i},${ali},${target},${method},${INPUT}," > info.txt


                    if [ "$ali" = "unaligned" ]; then
                        option="-u"
                    fi

                    while true; do date | tr '\n' '\t'; date +%s | tr '\n' '\t'; qstat | grep _f${i}_ | grep ' r ' | wc -l | tr '\n' '\t'; qstat | grep _f${i}_ | wc -l;sleep 60; done > node_count${SGE_TASK_ID}.txt &
                    FRACTAL ${option} -p ${placement} -i ${INPUT} -m ${method} -O "-l s_vmem=16G -l mem_req=16G -l short" -A "-l s_vmem=128G -l mem_req=128G -l short" -k 10000 -z 1000 -t 30000 -e -d 100 -j _f${i}_ -g -l 1000000 > fractal.out 2> fractal.err
                fi
                i=$(expr $i + 1)
            done
        done
    done
done