#!/bin/bash
#$ -S /bin/bash
#$ -N NK_0148
#$ -cwd
#$ -o /home/ha5653/backup/YachieLab/ikaken/NK_0148/mu_alpha_fractal_param/out
#$ -e /home/ha5653/backup/YachieLab/ikaken/NK_0148/mu_alpha_fractal_param/err

# parameters

# usage : qsub -l s_vmem=16G -l mem_req=16G -t 1-1800 /home/ha5653/backup/YachieLab/code/FRACTAL_exp/MyCodes/shell/NK_0148.mu_alpha.simulation.sh

source ${HOME}/.bash_profile

CWD="/home/ha5653/backup/YachieLab/ikaken/NK_0148/mu_alpha_fractal_param"
OUTPUT_DIR=${CWD}/TASK_ID${SGE_TASK_ID}
TREEDISTR="/home/ha5653/backup/YachieLab/code/FRACTAL_exp/MyCodes/R/NormRFdist.R"
PARTIAL="/home/ha5653/backup/YachieLab/code/FRACTAL_exp/MyCodes/python/PartialNRFdist.py"
PRESUME="/home/ha5653/backup/YachieLab/code/nkPRESUME/PRESUME.py"

ESUNUM=32768

TREE_TOPOLOGY="/home/ha5653/backup/YachieLab/ikaken/NK_0148/tree_topology/PRESUMEout/PRESUMEout.nwk"

### default ###
N=32768
m=1
g=20
L=1000
sigma=0
###############

i=1

for dir_id in 312 600; do
    for SUBSAMPLE in 50	59	71	84	100	119	141	168	200	238	283	336	400	476	566	673	800	951	1131	1345	1600; do
        for THRESHOLD in 50	71	100	141	200	283	400	566	800	1131	1600	2263	3200	4525	6400	9051	12800	18102	25600	36204; do
            if [ ${SUBSAMPLE} -le ${THRESHOLD} ];then
                for METHOD in rapidnjNJ raxmlMP fasttreeML; do 

                    if [ $i -eq ${SGE_TASK_ID} ];then
                        OUTPUT_DIR=${CWD}/TASK_ID${i}
                        mkdir ${OUTPUT_DIR}; cd ${OUTPUT_DIR}
                        mkdir ${METHOD}; cd ${METHOD}
                        cp $TREE_TOPOLOGY PRESUMEout.nwk

                        echo -n "${SGE_TASK_ID},${dir_id},${METHOD},${mode},${N},${SUBSAMPLE},${THRESHOLD}," > result.out

                        /usr/bin/time -f "%M,KB,%e,sec," timeout 10800 FRACTAL -i /home/ha5653/backup/YachieLab/ikaken/NK_0148/mu_alpha/TASK_ID${dir_id}/FRACTALin.fa.gz -m $METHOD -k ${SUBSAMPLE} -t ${THRESHOLD} -e 2>&1 1> /dev/null | tr -d "\n" >> result.out
                        python3 $PARTIAL PRESUMEout.nwk FRACTALout.nwk $TREEDISTR  >> result.out
                            
                        rm -r FRACTALout
                    fi
                    i=$(expr $i + 1)

                done
            fi
        done
    done
done