#!/bin/bash
#$ -S /bin/bash
#$ -N NK_0137
#$ -cwd
#$ -o /home/ha5653/backup/YachieLab/ikaken/NK_0137/out
#$ -e /home/ha5653/backup/YachieLab/ikaken/NK_0137/err

# qsub -l s_vmem=16G -l mem_req=16G -t 1-330 /home/ha5653/backup/YachieLab/code/FRACTAL_exp/MyCodes/shell/NK_0137.sh

source ${HOME}/.bash_profile

CODEDIR=/home/ha5653/backup/YachieLab/code

RESULTDIR=/home/ha5653/backup/YachieLab/ikaken/NK_0137/

SEQNUM=16000000

python3 ${CODEDIR}/nkPRESUME/PRESUME.py -n $SEQNUM -L 1 --constant 0 -s 0.05 --qsub --dist gamma2 --seed 0