#!/bin/bash
#$ -S /bin/bash
#$ -N NK_0142
#$ -cwd
#$ -o /home/ha5653/backup/YachieLab/ikaken/NK_0142/RunTimeSim/simulation/out
#$ -e /home/ha5653/backup/YachieLab/ikaken/NK_0142/RunTimeSim/simulation/err

CODEDIR=/home/ha5653/backup/YachieLab/code
CWD="/home/ha5653/backup/YachieLab/ikaken/NK_0142/RunTimeSim/simulation"
OUTPUT_DIR=${CWD}/TASK_ID${SGE_TASK_ID}

SUBSAMPLE=10000
THRESHOLD=10000
a=0
b=0.0034101310256636724

i=1
for tip_num in 1048576 	2097152 	4194304 	8388608 	16777216 	33554432 	67108864 	134217728 	268435456 	536870912 	1073741824 	2147483648 	4294967296 	8589934592 	17179869184;do
    for node_num in 1 2 5 10 20 30 40 50 60 70 80 90 100 150 200 300 400 600 800 1200 1600 2400 3200 4800 6400 10000 15000 20000;do
        for model in 1 2; do
            if [ $i -eq ${SGE_TASK_ID} ];then
                mkdir ${OUTPUT_DIR}
                cd ${OUTPUT_DIR}
                echo -n "${SGE_TASK_ID},${tip_num},${node_num},${model}," > result.out
                python3 ${CODEDIR}/FRACTAL_exp/MyCodes/python/RunTime_sim.py $node_num $tip_num $SUBSAMPLE $THRESHOLD $model $a $b >> result.out
            fi
            i=`expr $i + 1`
        done
    done
done