#!/bin/bash
#$ -S /bin/bash
#$ -N NK_0148
#$ -cwd
#$ -o /home/ha5653/backup/YachieLab/ikaken/NK_0148/mu_alpha/out
#$ -e /home/ha5653/backup/YachieLab/ikaken/NK_0148/mu_alpha/err

# parameters

# usage : qsub -l s_vmem=16G -l mem_req=16G -t 1-600 /home/ha5653/backup/YachieLab/code/FRACTAL_exp/MyCodes/shell/NK_0148.mu_alpha.simulation.sh

source ${HOME}/.bash_profile

CWD="/home/ha5653/backup/YachieLab/ikaken/NK_0148/mu_alpha"
OUTPUT_DIR=${CWD}/TASK_ID${SGE_TASK_ID}
TREEDISTR="/home/ha5653/backup/YachieLab/code/FRACTAL_exp/MyCodes/R/NormRFdist.R"
PARTIAL="/home/ha5653/backup/YachieLab/code/FRACTAL_exp/MyCodes/python/PartialNRFdist.py"
PRESUME="/home/ha5653/backup/YachieLab/code/nkPRESUME/PRESUME.py"

ESUNUM=32768
SUBSAMPLE=300
THRESHOLD=1000

TREE_TOPOLOGY="/home/ha5653/backup/YachieLab/ikaken/NK_0148/tree_topology/PRESUMEout/PRESUMEout.nwk"

i=1
for m in `seq 0.2 0.2 5`; do
    for g in 0.1 0.125892541 0.158489319 0.199526231 0.251188643 0.316227766 0.398107171 0.501187234 0.630957344 0.794328235 1 1.258925412 1.584893192 1.995262315 2.511886432 3.16227766 3.981071706 5.011872336 6.309573445 7.943282347 10 12.58925412 15.84893192 19.95262315; do
        if [ $i -eq ${SGE_TASK_ID} ];then
            mkdir $OUTPUT_DIR
            cd ${OUTPUT_DIR}
            python3 $PRESUME --gtrgamma GTR{0.03333/0.03333/0.03333/0.03333/0.03333/0.03333}+FU{0.25/0.25/0.25/0.25}+G4{${g}} -m $m --tree ${TREE_TOPOLOGY} --save; wait
            cat ${OUTPUT_DIR}/PRESUMEout/root.fa.gz ${OUTPUT_DIR}/PRESUMEout/PRESUMEout.fa.gz > ${OUTPUT_DIR}/FRACTALin.fa.gz
            echo "mu: $m, alpha: $g" > ${OUTPUT_DIR}/simulation.log
        fi
        i=`expr $i + 1`
    done
done