#!/bin/bash
#$ -S /bin/bash
#$ -N NK_0144
#$ -cwd
#$ -o /home/ha5653/backup/YachieLab/ikaken/NK_0144/out
#$ -e /home/ha5653/backup/YachieLab/ikaken/NK_0144/err

CODEDIR=/home/ha5653/backup/YachieLab/code
NK_0144_dir="/home/ha5653/backup/YachieLab/ikaken/NK_0144"
fractal_dir=${NK_0144_dir}/fractal
simulation_dir=${NK_0144_dir}/simulation
TREEDISTR=${CODEDIR}"/FRACTAL_exp/MyCodes/R/NormRFdist.R"
PARTIAL=${CODEDIR}"/FRACTAL_exp/MyCodes/python/PartialNRFdist.py"

mkdir ${fractal_dir}/

SUBSAMPLE=1000

i=1
for SEQNUM in 1024 16384 262144 4194304; do
    mkdir ${fractal_dir}/SEQNUM_${SEQNUM}/
    for ID in `seq 1 10`; do
        mkdir ${fractal_dir}/SEQNUM_${SEQNUM}/ID_${ID}/
        for method in rapidnjNJ raxmlMP fasttreeML; do
            mkdir ${fractal_dir}/SEQNUM_${SEQNUM}/ID_${ID}/method_${method}/
            for threshold in 1000 10000000; do
                if [ $i -eq ${SGE_TASK_ID} ];then

                    mkdir ${fractal_dir}/SEQNUM_${SEQNUM}/ID_${ID}/method_${method}/threshold_${threshold}
                    cd    ${fractal_dir}/SEQNUM_${SEQNUM}/ID_${ID}/method_${method}/threshold_${threshold}
                    echo -n "" >> result.out
 
                    if [ $(cat result.out | tr ',' '\n' | wc -l) -lt 11 ]; then

                        echo -n "${SGE_TASK_ID},$SEQNUM,$ID,$method,$threshold," > result.out
                        
                        # run fractal
                        TIP_SEQ_FILE="${simulation_dir}/mutation/SEQNUM_${SEQNUM}/ID_${ID}/PRESUMEout/PRESUMEout.fa.gz"
                        ROOT_SEQ_FILE="${simulation_dir}/mutation/SEQNUM_${SEQNUM}/ID_${ID}/PRESUMEout/root.fa.gz"
                        FRACTAL_IN_FILE="${fractal_dir}/SEQNUM_${SEQNUM}/ID_${ID}/method_${method}/threshold_${threshold}/PRESUMEout.root.fa.gz"
                        cat $ROOT_SEQ_FILE $TIP_SEQ_FILE | gunzip | gzip> $FRACTAL_IN_FILE

                        rm -r FRACTALout
                        /usr/bin/time -f "%M,KB,%e,sec," FRACTAL -l 10000000 -i $FRACTAL_IN_FILE -m ${method} -k ${SUBSAMPLE} -t ${threshold} -p ML -e -g -A "-l s_vmem=16G -l mem_req=16G -l ljob" 2>&1 1> fractal.out | tr -d "\n" >> result.out

                        # compress FRACTALout
                        tar cvf FRACTALout.tar.gz FRACTALout > /dev/null
                        rm -r FRACTALout
                    fi
                fi
                i=$(expr $i + 1)
            done
        done
    done
done  