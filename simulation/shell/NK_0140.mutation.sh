#!/bin/bash
#$ -S /bin/bash
#$ -N NK_0140
#$ -cwd
#$ -o /home/ha5653/backup/YachieLab/ikaken/NK_0140/out
#$ -e /home/ha5653/backup/YachieLab/ikaken/NK_0140/err


source ${HOME}/.bash_profile

CODEDIR=/home/ha5653/backup/YachieLab/code
RESULTDIR=/home/ha5653/backup/YachieLab/ikaken/NK_0140

del_pos_count=/home/ha5653/backup/YachieLab/ikaken/NK_0138/param/prob.del.0.34.txt
ins_pos_count=/home/ha5653/backup/YachieLab/ikaken/NK_0138/param/prob.ins.0.015.txt
del_len_count=/home/ha5653/backup/YachieLab/ikaken/NK_0120/count/len_count.del.txt
ins_len_count=/home/ha5653/backup/YachieLab/ikaken/NK_0120/count/len_count.ins.txt
initial_fasta=${RESULTDIR}/initial.fasta

i=1
for scale in 0.001 0.01 0.05 $(seq 0.1 0.1 1) $(seq 2 1 10); do 

    for No in $(seq 1 100); do

        if [ $i -eq ${SGE_TASK_ID} ]; then

            OUTDIR=${RESULTDIR}/scale${scale}_No${No}
            mkdir ${OUTDIR}
            cd ${OUTDIR}

            python3 ${CODEDIR}/FRACTAL_exp/MyCodes/python/NK_0138.py /home/ha5653/backup/YachieLab/ikaken/NK_0140/param/prob.${No}.csv $scale > ${OUTDIR}/pos_editprob.txt

            python3 ${CODEDIR}/nkPRESUME/PRESUME.py --tree /home/ha5653/backup/YachieLab/ikaken/NK_0138/PRESUMEout_topology/PRESUMEout.nwk -L 255 --editprofile ${OUTDIR}/pos_editprob.txt --inprob ${ins_pos_count} --inlength $ins_len_count --delprob ${del_pos_count} --dellength $del_len_count --seed ${SGE_TASK_ID} -f ${initial_fasta} --debug
        
        fi

        i=$(expr $i + 1)

    done

done