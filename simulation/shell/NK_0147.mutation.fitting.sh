#!/bin/bash
#$ -S /bin/bash
#$ -N NK_0147
#$ -cwd
#$ -o /home/ha5653/backup/YachieLab/ikaken/NK_0147/mutation_rate_fitting/out
#$ -e /home/ha5653/backup/YachieLab/ikaken/NK_0147/mutation_rate_fitting/err

PWD=/home/ha5653/backup/YachieLab/ikaken/NK_0147/mutation_rate_fitting
CODEDIR=/home/ha5653/backup/YachieLab/code

i=1
for mu in $(seq 0.01 0.01 1.00); do 

    if [ $i -eq ${SGE_TASK_ID} ];then

        OUTDIR=${PWD}/simulation/ID_${i}_${mu}

        mkdir ${OUTDIR}; cd ${OUTDIR}

        python3 ${CODEDIR}/nkPRESUME/PRESUME.py --tree /home/ha5653/backup/YachieLab/ikaken/NK_0147/mutation_rate_fitting/PRESUMEout_combined_length_tip_or_internal.nwk -L 1000 --gtrgamma GTR{0.821236/2.026286/1.547882/0.805315/3.662521/1.000000}+FU{0.263424/0.222548/0.315115/0.198914}+G4{0.572157} -m ${mu} --seed 0 --debug
    
    fi

    i=$(expr $i + 1)

done