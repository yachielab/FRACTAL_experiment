#!/bin/bash
#$ -S /bin/bash
#$ -N NK_0144
#$ -cwd
#$ -o /home/ha5653/backup/YachieLab/ikaken/NK_0144/out
#$ -e /home/ha5653/backup/YachieLab/ikaken/NK_0144/err

NK_0144_dir="/home/ha5653/backup/YachieLab/ikaken/NK_0144"
simulation_dir=${NK_0144_dir}/simulation
PRESUME="/home/ha5653/backup/YachieLab/code/nkPRESUME/PRESUME.py"
MUTATION_RATE=0.02

mkdir ${simulation_dir}/mutation/

i=1
for SEQNUM in 1024 16384 262144 4194304; do
    mkdir ${simulation_dir}/mutation/SEQNUM_${SEQNUM}/
    for ID in `seq 1 10`; do
        if [ $i -eq ${SGE_TASK_ID} ];then
            mkdir ${simulation_dir}/mutation/SEQNUM_${SEQNUM}/ID_${ID}
            cd    ${simulation_dir}/mutation/SEQNUM_${SEQNUM}/ID_${ID}
            # create topology
            python3 ${PRESUME} --constant $MUTATION_RATE -L 1000 --tree ${simulation_dir}/topology/SEQNUM_${SEQNUM}/PRESUMEout/PRESUMEout_combined.nwk --seed `expr ${SEQNUM} + ${ID}` --save
        fi
        i=$(expr $i + 1)
    done
done  