#!/bin/bash
#$ -S /bin/bash
#$ -N NK_0148
#$ -cwd
#$ -o /home/ha5653/backup/YachieLab/ikaken/NK_0148/sigma_mu/out
#$ -e /home/ha5653/backup/YachieLab/ikaken/NK_0148/sigma_mu/err

# parameters

# usage : "qsub -l s_vmem=16G -l mem_req=16G -t 1000 /home/ha5653/backup/YachieLab/code/FRACTAL_exp/MyCodes/shell/NK_0148.sigma_mu.simulation.sh
source ${HOME}/.bash_profile

CWD="/home/ha5653/backup/YachieLab/ikaken/NK_0148/sigma_mu"
OUTPUT_DIR=${CWD}/TASK_ID${SGE_TASK_ID}
TREEDISTR="/home/ha5653/backup/YachieLab/code/FRACTAL_exp/MyCodes/R/NormRFdist.R"
PARTIAL="/home/ha5653/backup/YachieLab/code/FRACTAL_exp/MyCodes/python/PartialNRFdist.py"
PRESUME="/home/ha5653/backup/YachieLab/code/nkPRESUME/PRESUME.py"

ESUNUM=32768
SUBSAMPLE=300
THRESHOLD=1000

i=1
for m in 0.01	0.012589254	0.015848932	0.019952623	0.025118864	0.031622777	0.039810717	0.050118723	0.063095734	0.079432823	0.1	0.125892541	0.158489319	0.199526231	0.251188643	0.316227766	0.398107171	0.501187234	0.630957344	0.794328235	1; do
    for sigma in `seq 0 0.1 3.0`; do
        if [ $i -eq ${SGE_TASK_ID} ];then
            mkdir $OUTPUT_DIR
            cd ${OUTPUT_DIR}

            while [ ! -e PRESUMEout/PRESUMEout.nwk ]; do
                python3 $PRESUME --gtrgamma GTR{0.03333/0.03333/0.03333/0.03333/0.03333/0.03333}+FU{0.25/0.25/0.25/0.25}+G4{20} -u 1000000000 -m ${m} -s ${sigma} -n ${ESUNUM} -L 1                                 --dist gamma2 >  presume.out
            done
            python3     $PRESUME --gtrgamma GTR{0.03333/0.03333/0.03333/0.03333/0.03333/0.03333}+FU{0.25/0.25/0.25/0.25}+G4{20} -u 1000000000 -m ${m} -s ${sigma} -n ${ESUNUM} -L 1000 --tree PRESUMEout/PRESUMEout.nwk --dist gamma2 >> presume.out
            esu=`cat presume.out| grep alive | tail -n1 | cut -d" " -f4`; wait
            cat PRESUMEout/root.fa.gz PRESUMEout/PRESUMEout.fa.gz > FRACTALin.fa.gz
        fi
        i=`expr $i + 1`
    done
done