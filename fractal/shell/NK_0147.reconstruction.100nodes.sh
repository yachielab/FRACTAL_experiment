#!/bin/bash
#$ -S /bin/bash
#$ -N NK_0147
#$ -cwd
#$ -o /home/ha5653/backup/YachieLab/ikaken/NK_0147/subclades/reconstruction_distributed/out
#$ -e /home/ha5653/backup/YachieLab/ikaken/NK_0147/subclades/reconstruction_distributed/err

# parameters

# usage : "qsub -l s_vmem=16G -l mem_req=16G -t 1-1000 ~/NK_0147.sh

source ${HOME}/.bash_profile

CWD="/home/ha5653/backup/YachieLab/ikaken/NK_0147/subclades/reconstruction_distributed"
OUTPUT_DIR=${CWD}/TASK_ID${SGE_TASK_ID}
TREEDISTR="/home/ha5653/backup/YachieLab/code/FRACTAL_exp/MyCodes/R/NormRFdist.R"
PARTIAL="/home/ha5653/backup/YachieLab/code/FRACTAL_exp/MyCodes/python/PartialNRFdist.py"

SUBSAMPLE=500
THRESHOLD=10000

DUMMY_ROOT="/home/ha5653/backup/YachieLab/ikaken/NK_0147/subclades/dummy_root.fa"

i=1

for METHOD in rapidnjNJ raxmlMP fasttreeML; do
    for t in $THRESHOLD; do
        for ID in 82 83; do
            if [ $i -eq ${SGE_TASK_ID} ];then
                mkdir ${OUTPUT_DIR}
                cd ${OUTPUT_DIR}
                INPUT_FASTA="/home/ha5653/backup/YachieLab/ikaken/NK_0147/subclades/fa/ID_${ID}_*.fa"
                TRUE_LINEAGE="/home/ha5653/backup/YachieLab/ikaken/NK_0147/subclades/nwk/ID_${ID}_*.nwk"
                cat $DUMMY_ROOT $INPUT_FASTA | gzip > FRACTALin.fa.gz
                cp ${TRUE_LINEAGE} answer.nwk

                Nseq=`cat ${INPUT_FASTA} | grep '>' | wc -l`
                echo -n "${SGE_TASK_ID},${ID},${Nseq},${METHOD},${t}," >> ${OUTPUT_DIR}/result.out
                /usr/bin/time -f "%M,KB,%e,sec," FRACTAL -m $METHOD -i FRACTALin.fa.gz -k $SUBSAMPLE -t $t -e -d 100 2>&1 1> /dev/null | tr -d "\n" >> ${OUTPUT_DIR}/result.out
                python3 $PARTIAL FRACTALout.nwk answer.nwk $TREEDISTR >>${OUTPUT_DIR}/result.out
                tar cvf FRACTALout.tar.gz FRACTALout > /dev/null
                rm -r FRACTALout
            fi
            i=`expr $i + 1`
        done
    done
done