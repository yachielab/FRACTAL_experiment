#!/bin/bash
#$ -S /bin/bash
#$ -N NK_0148
#$ -cwd
#$ -o /home/ha5653/backup/YachieLab/ikaken/NK_0148/L/out
#$ -e /home/ha5653/backup/YachieLab/ikaken/NK_0148/L/err

# parameters

# usage : qsub -l s_vmem=16G -l mem_req=16G -t 1-1800 /home/ha5653/backup/YachieLab/code/FRACTAL_exp/MyCodes/shell/NK_0148.L.simulation.sh

source ${HOME}/.bash_profile

CWD="/home/ha5653/backup/YachieLab/ikaken/NK_0148/L"
OUTPUT_DIR=${CWD}/TASK_ID${SGE_TASK_ID}
TREEDISTR="/home/ha5653/backup/YachieLab/code/FRACTAL_exp/MyCodes/R/NormRFdist.R"
PARTIAL="/home/ha5653/backup/YachieLab/code/FRACTAL_exp/MyCodes/python/PartialNRFdist.py"
PRESUME="/home/ha5653/backup/YachieLab/code/nkPRESUME/PRESUME.py"

ESUNUM=32768
SUBSAMPLE=300
THRESHOLD=1000

### default ###
N=32768
m=1
g=20
L=1000
sigma=0
###############

TREE_TOPOLOGY="/home/ha5653/backup/YachieLab/ikaken/NK_0148/tree_topology/PRESUMEout/PRESUMEout.nwk"

dir_id=1
i=1
for L in `seq 5 5 1500`; do
    for METHOD in rapidnjNJ raxmlMP fasttreeML; do 
        for mode in fractal original; do
            if [ $i -eq ${SGE_TASK_ID} ];then
                OUTPUT_DIR=${CWD}/TASK_ID${dir_id}
                cd ${OUTPUT_DIR}
                mkdir ${METHOD}_${mode}; cd ${METHOD}_${mode}
                cp $TREE_TOPOLOGY PRESUMEout.nwk
                echo -n "${SGE_TASK_ID},${dir_id},${METHOD},${mode},${N},${m},${g},${L},${sigma}," > result.out
                if   [ $mode = "original" ]; then
                    /usr/bin/time -f "%M,KB,%e,sec," timeout 10800 FRACTAL -i ${OUTPUT_DIR}/FRACTALin.fa.gz -m $METHOD -k 100          -t 10000000     -e 2>&1 1> /dev/null | tr -d "\n" >> result.out
                    python3 $PARTIAL PRESUMEout.nwk FRACTALout.nwk $TREEDISTR  >> result.out                    
                elif [ $mode = "fractal" ]; then
                    /usr/bin/time -f "%M,KB,%e,sec," timeout 10800 FRACTAL -i ${OUTPUT_DIR}/FRACTALin.fa.gz -m $METHOD -k ${SUBSAMPLE} -t ${THRESHOLD} -e 2>&1 1> /dev/null | tr -d "\n" >> result.out
                    python3 $PARTIAL PRESUMEout.nwk FRACTALout.nwk $TREEDISTR  >> result.out                
                fi
                rm -r FRACTALout
            fi
            i=$(expr $i + 1)      
        done
    done
    dir_id=$(expr $dir_id + 1)
done