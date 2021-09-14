#!/bin/bash
#$ -S /bin/bash
#$ -N NK_0135
#$ -cwd
#$ -o /home/ha5653/backup/YachieLab/ikaken/NK_0135/out
#$ -e /home/ha5653/backup/YachieLab/ikaken/NK_0135/err

PWD=/home/ha5653/backup/YachieLab/ikaken/NK_0135/
PRESUME='/home/ha5653/backup/YachieLab/code/nkPRESUME/PRESUME.py'

N=370

i=0

# 6
DIR=sigma_r_gamma2
mkdir ${PWD}${DIR}
for sigma in $(seq 0 0.0001 0.03); do 
    i=$(expr $i + 1)
    if [ $i -eq ${SGE_TASK_ID} ]; then
        mkdir ${PWD}${DIR}/${DIR}_${sigma}
        timeout 21600 python3 ${PRESUME} --constant 0 -L 1 -n $N -s $sigma --output ${PWD}${DIR}/${DIR}_${sigma} --dist gamma2 > ${PWD}${DIR}/${DIR}_${sigma}/presume.out 2> ${PWD}${DIR}/${DIR}_${sigma}$(echo $cv | tr -d '-')/dM_d.tsv
    fi
done
