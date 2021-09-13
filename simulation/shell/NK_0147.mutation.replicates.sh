#!/bin/bash
#$ -S /bin/bash
#$ -N NK_0147
#$ -cwd
#$ -o /home/ha5653/backup/YachieLab/ikaken/NK_0147/replicates/out
#$ -e /home/ha5653/backup/YachieLab/ikaken/NK_0147/replicates/err

PWD=/home/ha5653/backup/YachieLab/ikaken/NK_0147/replicates
CODEDIR=/home/ha5653/backup/YachieLab/code

i=1
for mu in 0.077; do 

    for rep in 1 2 3 4; do 

        if [ $i -eq ${SGE_TASK_ID} ];then

            OUTDIR=${PWD}/simulation/ID_${i}_${mu}_${rep}

            mkdir ${OUTDIR}; cd ${OUTDIR}

            python3 ${CODEDIR}/nkPRESUME/PRESUME.py --tree /home/ha5653/backup/YachieLab/ikaken/NK_0147/mutation_rate_fitting/PRESUMEout_combined_length_tip_or_internal.nwk -L 1000 --gtrgamma GTR{0.821236/2.026286/1.547882/0.805315/3.662521/1.000000}+FU{0.263424/0.222548/0.315115/0.198914}+G4{0.572157} -m ${mu} --seed ${rep} --debug
        
        fi

    i=$(expr $i + 1)

    done

done