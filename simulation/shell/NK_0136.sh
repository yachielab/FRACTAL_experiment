#!/bin/bash
#$ -S /bin/bash
#$ -N NK_0136
#$ -cwd
#$ -o /home/ha5653/backup/YachieLab/ikaken/NK_0136/out
#$ -e /home/ha5653/backup/YachieLab/ikaken/NK_0136/err


# qsub -l s_vmem=16G -l mem_req=16G -t 1-330 /home/ha5653/backup/YachieLab/code/FRACTAL_exp/MyCodes/shell/NK_0136.sh

source ${HOME}/.bash_profile

CODEDIR=/home/ha5653/backup/YachieLab/code

RESULTDIR=/home/ha5653/backup/YachieLab/ikaken/NK_0136/

del_pos_count=/home/ha5653/backup/YachieLab/ikaken/NK_0120/count/midpos_count.del.txt
ins_pos_count=/home/ha5653/backup/YachieLab/ikaken/NK_0120/count/pos_count.ins.txt
del_len_count=/home/ha5653/backup/YachieLab/ikaken/NK_0120/count/len_count.del.txt
ins_len_count=/home/ha5653/backup/YachieLab/ikaken/NK_0120/count/len_count.ins.txt

preparation_code=${CODEDIR}/FRACTAL_exp/MyCodes/python/NK_0110.py

SEQNUM=4000
L=255 

cd ${RESULTDIR}

i=1
for prob in  `seq 0.0001 0.0001 0.01` `seq 0.001 0.001 0.1` `seq 0.01 0.01 1`; do

    if [ $i -eq ${SGE_TASK_ID} ];then

        mkdir sim${i}

        cd sim${i}

        python3 $preparation_code $del_pos_count ${prob} ${L} > prob.scaled.txt

        for _ in `seq 1 ${L}`; do echo 0 >> prob.zero.txt; done

        # prep tree topology
        python3 ${CODEDIR}/nkPRESUME/PRESUME.py -n $SEQNUM -L 1 --constant 0 -s 0.05 --seed ${i} --dist gamma2

        mv PRESUMEout PRESUMEout_template

        # insertion only simulation
        python3 ${CODEDIR}/nkPRESUME/PRESUME.py --tree PRESUMEout_template/PRESUMEout.nwk -L ${L} --constant 0 -s 0.05 --inprob prob.scaled.txt --inlength $ins_len_count --delprob prob.zero.txt --dellength $del_len_count --seed ${i}
        mv PRESUMEout PRESUMEout_ins_${prob}

        # deletion only simulation
        python3 ${CODEDIR}/nkPRESUME/PRESUME.py --tree PRESUMEout_template/PRESUMEout.nwk -L ${L} --constant 0 -s 0.05 --inprob prob.zero.txt --inlength $ins_len_count --delprob prob.scaled.txt --dellength $del_len_count --seed ${i}
        mv PRESUMEout PRESUMEout_del_${prob}

    fi
    
    i=`expr $i + 1`

done
