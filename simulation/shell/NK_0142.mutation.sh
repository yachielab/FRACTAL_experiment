#!/bin/bash
#$ -S /bin/bash
#$ -N NK_0142
#$ -cwd
#$ -o /home/ha5653/backup/YachieLab/ikaken/NK_0142/out
#$ -e /home/ha5653/backup/YachieLab/ikaken/NK_0142/err

# qsub -l s_vmem=512G -l mem_req=512G -l lmem /home/ha5653/backup/YachieLab/code/FRACTAL_exp/MyCodes/shell/NK_0142.sh

source ${HOME}/.bash_profile

CODEDIR=/home/ha5653/backup/YachieLab/code

RESULTDIR=/home/ha5653/backup/YachieLab/ikaken/NK_0142/

cd $RESULTDIR

SEQNUM=200000000

#python3 ${CODEDIR}/nkPRESUME/PRESUME.py -n $SEQNUM -L 1 --constant 0 -s 0.05 --qsub --dist gamma2 --seed 0

#mv PRESUMEout PRESUMEout_topology

python3 ${CODEDIR}/nkPRESUME/PRESUME.py --tree PRESUMEout_topology/PRESUMEout_combined.nwk -L 1000 --gtrgamma GTR{0.03333/0.03333/0.03333/0.03333/0.03333/0.03333}+FU{0.25/0.25/0.25/0.25}+G4{20} -m 0.50 --qsub --seed 0 --debug