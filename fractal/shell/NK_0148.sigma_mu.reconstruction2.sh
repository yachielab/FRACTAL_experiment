#!/bin/bash
#$ -S /bin/bash
#$ -N NK_0148
#$ -o /home/ha5653/backup/YachieLab/ikaken/NK_0148/sigma_mu/out
#$ -e /home/ha5653/backup/YachieLab/ikaken/NK_0148/sigma_mu/err

# parameters

# usage : "qsub -l s_vmem=16G -l mem_req=16G -t 1000 /home/ha5653/backup/YachieLab/code/FRACTAL_exp/MyCodes/shell/NK_0148.sigma_mu.reconstruction.sh

source ${HOME}/.bash_profile

CWD="/home/ha5653/backup/YachieLab/ikaken/NK_0148/sigma_mu"
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

rep=2

dir_id=1
i=1
for m in 0.01	0.012589254	0.015848932	0.019952623	0.025118864	0.031622777	0.039810717	0.050118723	0.063095734	0.079432823	0.1	0.125892541	0.158489319	0.199526231	0.251188643	0.316227766	0.398107171	0.501187234	0.630957344	0.794328235	1; do
    for sigma in `seq 0 0.1 3.0`; do
        for METHOD in rapidnjNJ raxmlMP fasttreeML; do 
            for mode in fractal original; do
                if [ $i -eq ${SGE_TASK_ID} ];then
                    OUTPUT_DIR=${CWD}/TASK_ID${dir_id}_${rep}
                    mkdir ${OUTPUT_DIR}; cd ${OUTPUT_DIR}
                    mkdir ${METHOD}_${mode}; cd ${METHOD}_${mode}

                    TREE_TOPOLOGY=${CWD}/TASK_ID${dir_id}/PRESUMEout/PRESUMEout.nwk
                    cp $TREE_TOPOLOGY PRESUMEout.nwk

                    N=$(cat ${CWD}/TASK_ID${dir_id}/PRESUMEout/PRESUMEout.fa.gz | gunzip | grep '>' | wc -l)

                    echo -n "${SGE_TASK_ID},${dir_id},${METHOD},${mode},${N},${m},${g},${L},${sigma}," > result.out
                    if   [ $mode = "original" ]; then
                        /usr/bin/time -f "%M,KB,%e,sec," timeout 86400 FRACTAL -i ${CWD}/TASK_ID${dir_id}/FRACTALin.fa.gz -m $METHOD -k 100          -t 10000000     -r ${rep} -e 2>&1 1> /dev/null | tr -d "\n" >> result.out
                        python3 $PARTIAL PRESUMEout.nwk FRACTALout.nwk $TREEDISTR  >> result.out                    
                    elif [ $mode = "fractal" ]; then
                        /usr/bin/time -f "%M,KB,%e,sec," timeout 86400 FRACTAL -i ${CWD}/TASK_ID${dir_id}/FRACTALin.fa.gz -m $METHOD -k ${SUBSAMPLE} -t ${THRESHOLD} -r ${rep} -e 2>&1 1> /dev/null | tr -d "\n" >> result.out
                        python3 $PARTIAL PRESUMEout.nwk FRACTALout.nwk $TREEDISTR  >> result.out                
                    fi
                    rm -r FRACTALout
                fi
                i=$(expr $i + 1)      
            done
        done
        dir_id=$(expr $dir_id + 1)
    done
done