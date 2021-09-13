#!/bin/bash
#$ -S /bin/bash
#$ -N NK_0148
#$ -cwd
#$ -o /home/ha5653/backup/YachieLab/ikaken/NK_0148/L/out
#$ -e /home/ha5653/backup/YachieLab/ikaken/NK_0148/L/err

# parameters

# usage : qsub -l s_vmem=16G -l mem_req=16G -t 1-600 /home/ha5653/backup/YachieLab/code/FRACTAL_exp/MyCodes/shell/NK_0148.L.simulation.sh

source ${HOME}/.bash_profile

CWD="/home/ha5653/backup/YachieLab/ikaken/NK_0148/L"
OUTPUT_DIR=${CWD}/TASK_ID${SGE_TASK_ID}
TREEDISTR="/home/ha5653/backup/YachieLab/code/FRACTAL_exp/MyCodes/R/NormRFdist.R"
PARTIAL="/home/ha5653/backup/YachieLab/code/FRACTAL_exp/MyCodes/python/PartialNRFdist.py"
PRESUME="/home/ha5653/backup/YachieLab/code/nkPRESUME/PRESUME.py"

ESUNUM=32768
SUBSAMPLE=300
THRESHOLD=1000

TREE_TOPOLOGY="/home/ha5653/backup/YachieLab/ikaken/NK_0148/tree_topology/PRESUMEout/PRESUMEout.nwk"

m=2.6
g=20

i=1
for L in `seq 5 5 1500`; do
    if [ $i -eq ${SGE_TASK_ID} ];then
        mkdir $OUTPUT_DIR
        cd ${OUTPUT_DIR}
        python3 $PRESUME -L $L --gtrgamma GTR{0.03333/0.03333/0.03333/0.03333/0.03333/0.03333}+FU{0.25/0.25/0.25/0.25}+G4{${g}} -m $m --tree ${TREE_TOPOLOGY} --save; wait
        cat ${OUTPUT_DIR}/PRESUMEout/root.fa.gz ${OUTPUT_DIR}/PRESUMEout/PRESUMEout.fa.gz > ${OUTPUT_DIR}/FRACTALin.fa.gz
        echo "L: $L" > ${OUTPUT_DIR}/simulation.log
    fi
    i=`expr $i + 1`
done