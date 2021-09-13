#!/bin/bash
#$ -S /bin/bash
#$ -N NK_0149
#$ -cwd
#$ -o /home/ha5653/backup/YachieLab/ikaken/NK_0149/fractal/out
#$ -e /home/ha5653/backup/YachieLab/ikaken/NK_0149/fractal/err


# qsub -l s_vmem=16G -l mem_req=16G -t 1-15 /home/ha5653/backup/YachieLab/code/FRACTAL_exp/MyCodes/shell/NK_0149_FRACTAL.sh

source ${HOME}/.bash_profile

CODEDIR=/home/ha5653/backup/YachieLab/code

RESULTDIR=/home/ha5653/backup/YachieLab/ikaken/NK_0149/fractal

INPUTFA=/home/ha5653/backup/YachieLab/ikaken/NK_0149/simulation/fractalin/FRACTALin.fa.gz

TRUE_LINEAGE=/home/ha5653/backup/YachieLab/ikaken/NK_0149/simulation/lineage/run1.newick


cd ${RESULTDIR}


i=1
for method in rapidnjNJ raxmlMP fasttreeML; do
    
    for r in 1 2 3 4 5; do

        if [ $i -eq ${SGE_TASK_ID} ];then

            OUTPUT_DIR=${RESULTDIR}/TASK_${i}

            mkdir ${OUTPUT_DIR}; cd ${OUTPUT_DIR}
         
            JOBNAME="_${i}_"

            Nseq=$(cat ${INPUTFA} | gunzip | grep '>' | wc -l)

            echo -n "${i},${Nseq},${method},${r}," > ${OUTPUT_DIR}/result.out

            while true; do date | tr '\n' '\t'; date +%s | tr '\n' '\t'; qstat | grep ${JOBNAME} | grep ' r ' | wc -l | tr '\n' '\t'; qstat | grep ${JOBNAME} | wc -l; sleep 60; done > ${OUTPUT_DIR}/node_count.txt &

            /usr/bin/time -f "%M,KB,%e,sec," FRACTAL -i ${INPUTFA} -k 1000 -t 10000 -d 100 -u -e -O "-l s_vmem=16G -l mem_req=16G" -m ${method} -r ${r} -j ${JOBNAME} -f FRACTAL_${method}_${r} -l 1000000 2>&1 1> ${OUTPUT_DIR}/fractal.out | tr -d "\n" >> ${OUTPUT_DIR}/result.out

            cp ${TRUE_LINEAGE} SIMULATION.nwk

            python3 $PARTIAL FRACTALout.nwk SIMULATION.nwk $TREEDISTR >> ${OUTPUT_DIR}/result.out

            tar cvf FRACTALout.tar.gz FRACTALout
            rm -r FRACTALout

        fi

        i=`expr $i + 1`
    
    done

done