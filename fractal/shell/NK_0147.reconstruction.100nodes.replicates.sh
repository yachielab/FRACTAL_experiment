#!/bin/bash
#$ -S /bin/bash
#$ -N NK_0147
#$ -cwd
#$ -o /home/ha5653/backup/YachieLab/ikaken/NK_0147/replicates/reconstruction/out
#$ -e /home/ha5653/backup/YachieLab/ikaken/NK_0147/replicates/reconstruction/err

# parameters

# usage : "qsub -l s_vmem=16G -l mem_req=16G -t 1-1000 ~/NK_0147.sh

source ${HOME}/.bash_profile

TREEDISTR="/home/ha5653/backup/YachieLab/code/FRACTAL_exp/MyCodes/R/NormRFdist.R"
PARTIAL="/home/ha5653/backup/YachieLab/code/FRACTAL_exp/MyCodes/python/PartialNRFdist.py"

SUBSAMPLE=500
THRESHOLD=10000

i=1

CWD="/home/ha5653/backup/YachieLab/ikaken/NK_0147/replicates/reconstruction"

OUTPUT_DIR=${CWD}/TASK_ID${SGE_TASK_ID}

for METHOD in rapidnjNJ raxmlMP fasttreeML; do
    for t in $THRESHOLD; do
        for sim_id in 1 2 3 4; do
            for fractal_id in 1 2 3 4 5; do  
                if [ $i -eq ${SGE_TASK_ID} ];then
                    mkdir ${OUTPUT_DIR}
                    cd ${OUTPUT_DIR}
                    INPUT_FASTA="/home/ha5653/backup/YachieLab/ikaken/NK_0147/replicates/simulation/ID_${sim_id}_0.077_${sim_id}/PRESUMEout/PRESUMEout.fa.gz"
                    TRUE_LINEAGE="/home/ha5653/backup/YachieLab/ikaken/NK_0147/replicates/simulation/ID_${sim_id}_0.077_${sim_id}/PRESUMEout/PRESUMEout.nwk"
                    cat $INPUT_FASTA | gunzip | sed 's/728_249b/root/g' | gzip > FRACTALin.fa.gz
                    cp ${TRUE_LINEAGE} answer.nwk

                    Nseq=`cat ${INPUT_FASTA} | gunzip | grep '>' | wc -l`
                    echo -n "${SGE_TASK_ID},${ID},${Nseq},${METHOD},${t},${sim_id},${fractal_id}" >> ${OUTPUT_DIR}/result.out
                    /usr/bin/time -f "%M,KB,%e,sec," FRACTAL -m $METHOD -i FRACTALin.fa.gz -k $SUBSAMPLE -t $t -e -d 100 -r ${fractal_id} -j f${sim_id}_${fractal_id} 2>&1 1> /dev/null | tr -d "\n" >> ${OUTPUT_DIR}/result.out
                    python3 $PARTIAL FRACTALout.nwk answer.nwk $TREEDISTR >>${OUTPUT_DIR}/result.out
                    #rm -r FRACTALout
                    tar cvf FRACTALout.tar.gz FRACTALout
                fi
                i=`expr $i + 1`
            done
        done
    done
done