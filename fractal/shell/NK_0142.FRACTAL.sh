#!/bin/bash
#$ -S /bin/bash
#$ -N NK_0142
#$ -cwd
#$ -o /home/ha5653/backup/YachieLab/ikaken/NK_0142/out
#$ -e /home/ha5653/backup/YachieLab/ikaken/NK_0142/err

# qsub -l s_vmem=16G -l mem_req=16G -l ljob /home/ha5653/backup/YachieLab/code/FRACTAL_exp/MyCodes/shell/NK_0142.FRACTAL.sh

source ${HOME}/.bash_profile

CODEDIR=/home/ha5653/backup/YachieLab/code

RESULTDIR=/home/ha5653/backup/YachieLab/ikaken/NK_0142/FRACTAL

FRACTAL_IN_FILE=${RESULTDIR}/FRACTALin.fa.gz

JOBNAME=fig2fr

cd $RESULTDIR

cat /home/ha5653/backup/YachieLab/ikaken/NK_0142/PRESUMEout/root.fa.gz /home/ha5653/backup/YachieLab/ikaken/NK_0142/PRESUMEout/PRESUMEout.fa.gz > ${FRACTAL_IN_FILE}

while true; do date | tr '\n' '\t'; date +%s | tr '\n' '\t'; qstat | grep ${JOBNAME} | grep ' r ' | wc -l | tr '\n' '\t'; qstat | grep ${JOBNAME} | wc -l; sleep 60; done > ${RESULTDIR}/node_count.txt &

FRACTAL -i $FRACTAL_IN_FILE -m rapidnjNJ -k 10000 -z 1000 -t 10000 -p ML -e -I "-l s_vmem=128G -l mem_req=128G -l ljob" -O "-l s_vmem=16G -l mem_req=16G" -A "-l s_vmem=512G -l mem_req=512G -l lmem" -d 300 -j ${JOBNAME} -g -l 1000000 > fractal.out 2> fractal.err