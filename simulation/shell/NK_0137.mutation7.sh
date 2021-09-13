#!/bin/bash
#$ -S /bin/bash
#$ -N NK_0137
#$ -cwd
#$ -o /home/ha5653/backup/YachieLab/ikaken/NK_0137/out
#$ -e /home/ha5653/backup/YachieLab/ikaken/NK_0137/err

# qsub -l s_vmem=128G -l mem_req=128G -t 1-5 /home/ha5653/backup/YachieLab/code/FRACTAL_exp/MyCodes/shell/NK_0137.mutation4.sh

source ${HOME}/.bash_profile

CODEDIR=/home/ha5653/backup/YachieLab/code

RESULTDIR=/home/ha5653/backup/YachieLab/ikaken/NK_0137

del_pos_count=${RESULTDIR}/prob/prob.del.0.17.txt
ins_pos_count=${RESULTDIR}/prob/prob.ins.0.0075.txt
del_len_count=/home/ha5653/backup/YachieLab/ikaken/NK_0120/count/len_count.del.txt
ins_len_count=/home/ha5653/backup/YachieLab/ikaken/NK_0120/count/len_count.ins.txt
initial_fasta=/home/ha5653/backup/YachieLab/ikaken/NK_0137/initial.fasta

#####################
SEQNUM=16000000
#####################

i=1
for scale in 0; do
    if [ $i -eq ${SGE_TASK_ID} ]; then
        mkdir ${RESULTDIR}/SEQ4/scale${scale}
        for ID in $(seq 66 100); do
            mkdir ${RESULTDIR}/SEQ4/scale${scale}/chunk${ID}
            cd    ${RESULTDIR}/SEQ4/scale${scale}/chunk${ID}

            #python3 ${CODEDIR}/FRACTAL_exp/MyCodes/python/NK_0138.py /home/ha5653/backup/YachieLab/ikaken/NK_0140/param/prob.${ID}.csv $scale > ${RESULTDIR}/SEQ4/scale${scale}/chunk${ID}/pos_editprob.txt

            #python3 ${CODEDIR}/nkPRESUME/PRESUME.py --tree ${RESULTDIR}/decomp_tree/ -L 255 --editprofile ${RESULTDIR}/SEQ4/scale${scale}/chunk${ID}/pos_editprob.txt --inprob ${ins_pos_count} --inlength $ins_len_count --delprob ${del_pos_count} --dellength $del_len_count --seed ${ID} -f ${initial_fasta} --debug --qsub

            python3 ${CODEDIR}/nkPRESUME/PRESUME.py --tree ${RESULTDIR}/decomp_tree/ -L 255 --constant 0 --inprob ${ins_pos_count} --inlength $ins_len_count --delprob ${del_pos_count} --dellength $del_len_count --seed ${ID} -f ${initial_fasta} --debug --qsub

            rm -r ${RESULTDIR}/SEQ4/scale${scale}/chunk${ID}/PRESUMEout/intermediate/
        done
    fi
    i=$(expr $i + 1)
done