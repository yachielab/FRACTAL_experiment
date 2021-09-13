#!/bin/bash
#$ -S /bin/bash
#$ -N NK_0151
#$ -cwd
#$ -o /home/ha5653/backup/YachieLab/ikaken/NK_0151/out
#$ -e /home/ha5653/backup/YachieLab/ikaken/NK_0151/err

# parameters

# usage : "qsub -l s_vmem=16G -l mem_req=16G -t 1-1000 /home/ha5653/backup/YachieLab/code/FRACTAL/MyCodes/shell/NK_0151.sh

source ${HOME}/.bash_profile

CWD="/home/ha5653/backup/YachieLab/ikaken/NK_0151/fractal"
OUTPUT_DIR=${CWD}/TASK_ID${SGE_TASK_ID}
TREEDISTR="/home/ha5653/backup/YachieLab/code/FRACTAL_exp/MyCodes/R/NormRFdist.R"
PARTIAL="/home/ha5653/backup/YachieLab/code/FRACTAL_exp/MyCodes/python/PartialNRFdist.py"

ESUNUM=13897
i=1

SEQ_FILE=/home/ha5653/backup/YachieLab/ikaken/NK_0151/fasta/LTPs132_SSU_aligned.fasta.trim.dereplicated.fractalin.fa
RAND_SEQ_FILE=/home/ha5653/backup/YachieLab/ikaken/NK_0151/fasta/LTPs132_SSU_aligned.fasta.trim.dereplicated.fractalin.shuffled.fa
TREE_FILE=/home/ha5653/backup/YachieLab/ikaken/NK_0151/true_lineage/LTPs132_SSU_tree.dereplicated.newick

for random_ratio in $(seq 0 0.2 5.0); do
    for SUBSAMPLE in 500; do
        for THRESHOLD in 100000 4000; do
            for METHOD in rapidnjNJ raxmlMP fasttreeML; do
                if [ $i -eq ${SGE_TASK_ID} ];then
                    mkdir $OUTPUT_DIR; cd ${OUTPUT_DIR}
                    if [ "$random_ratio" = "0.0" ]; then
                        cat ${SEQ_FILE} > fractalin.mixed.fa
                    else
                        (cat ${SEQ_FILE}; seqkit sample -p ${random_ratio} ${RAND_SEQ_FILE}) > fractalin.mixed.fa
                    fi
                    N=$(cat fractalin.mixed.fa | grep '>' | wc -l)
                    cp ${TREE_FILE} ref_lineage.nwk
                    echo -n "${SGE_TASK_ID},${METHOD},${N},${SUBSAMPLE},${THRESHOLD}," >> result.out
                    
                    /usr/bin/time -f "%M,KB,%e,sec," FRACTAL -m ${METHOD} -i fractalin.mixed.fa -k ${SUBSAMPLE} -t ${THRESHOLD} -e 2>&1 1> fractal.out | tr -d "\n" >> result.out
                    
                    python3 $PARTIAL ref_lineage.nwk FRACTALout.nwk $TREEDISTR         >> result.out
                    tar cvf FRACTALout.tar.gz FRACTALout/
                    rm -r FRACTALout
                fi
                i=`expr $i + 1`
            done
        done
    done
done