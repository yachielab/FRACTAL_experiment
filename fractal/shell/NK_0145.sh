#!/bin/bash
#$ -S /bin/bash
#$ -N NK_0145
#$ -cwd
#$ -o /home/ha5653/backup/YachieLab/ikaken/NK_0145/out
#$ -e /home/ha5653/backup/YachieLab/ikaken/NK_0145/err

# parameters

# usage : "qsub -l s_vmem=16G -l mem_req=16G -t 1-1000 /home/ha5653/backup/YachieLab/code/FRACTAL/MyCodes/shell/NK_0145.sh

source ${HOME}/.bash_profile

CWD="/home/ha5653/backup/YachieLab/ikaken/NK_0145/result"
OUTPUT_DIR=${CWD}/TASK_ID${SGE_TASK_ID}
TREEDISTR="/home/ha5653/backup/YachieLab/code/FRACTAL_exp/MyCodes/R/NormRFdist.R"
PARTIAL="/home/ha5653/backup/YachieLab/code/FRACTAL_exp/MyCodes/python/PartialNRFdist.py"

ESUNUM=13897
i=1

SEQ_FILE=/home/ha5653/backup/YachieLab/ikaken/NK_0145/fractalin/LTPs132_SSU_aligned.fasta.trim.dereplicated.fractalin.fa
TREE_FILE=/home/ha5653/backup/YachieLab/ikaken/NK_0145/dereplicated/LTPs132_SSU_tree.dereplicated.newick

for SUBSAMPLE in 1000	950	900	850	800	750	700	650	600	550	500 450 400 350 300 250 200 150 100 50; do
    for THRESHOLD in 14000	13500	13000	12500	12000	11500	11000	10500	10000	9500	9000	8500	8000	7500	7000	6500	6000	5500	5000	4500	4000	3500	3000	2500	2000; do
        if [ ${SUBSAMPLE} -le ${THRESHOLD} ];then
            for METHOD in rapidnjNJ raxmlMP fasttreeML; do
                if [ $i -eq ${SGE_TASK_ID} ];then
                    mkdir $OUTPUT_DIR; cd ${OUTPUT_DIR}

                    cp ${TREE_FILE} ref_lineage.nwk
                    echo -n "${SUBSAMPLE},${THRESHOLD},${SGE_TASK_ID},${METHOD},${ESUNUM}," >> result.out
                    /usr/bin/time -f "%M,KB,%e,sec," FRACTAL -m ${METHOD} -i ${SEQ_FILE} -k ${SUBSAMPLE} -t ${THRESHOLD} -e 2>&1 1> /dev/null | tr -d "\n" >> result.out
                    python3 $PARTIAL ref_lineage.nwk FRACTALout.nwk $TREEDISTR >> result.out
                    rm -r ${METHOD}/FRACTALout/nodes
                fi
                i=`expr $i + 1`
            done
        fi
    done
done