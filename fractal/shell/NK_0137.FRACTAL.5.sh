#!/bin/bash
#$ -S /bin/bash
#$ -N NK_0137
#$ -cwd
#$ -o /home/ha5653/backup/YachieLab/ikaken/NK_0137/FRACTAL/out
#$ -e /home/ha5653/backup/YachieLab/ikaken/NK_0137/FRACTAL/err

CODEDIR=/home/ha5653/backup/YachieLab/code
CWD="/home/ha5653/backup/YachieLab/ikaken"
TREEDISTR=${CODEDIR}"/FRACTAL_exp/MyCodes/R/NormRFdist.R"
PARTIAL=${CODEDIR}"/FRACTAL_exp/MyCodes/python/PartialNRFdist.py"
editing_rate=0.025
RESULT_DIR=$CWD/NK_0137/FRACTAL/FRACTAL_${editing_rate}
JOB_NAME=f5_${SGE_TASK_ID}_

mkdir $RESULT_DIR
cd $RESULT_DIR

i=1
for Nchunks in 15 18 20 30 40 50 100 _tail15 _tail18 _tail20 _tail30 _tail40 _tail50; do
    if [ $i -eq ${SGE_TASK_ID} ]; then

        mkdir ${RESULT_DIR}/${Nchunks}chunks
        cd ${RESULT_DIR}/${Nchunks}chunks

        while true; do date | tr '\n' '\t'; date +%s | tr '\n' '\t'; qstat | grep ${JOB_NAME} | grep ' r ' | wc -l | tr '\n' '\t'; qstat | grep ${JOB_NAME} | wc -l;sleep 60; done > node_count.txt &
        
        TIP_SEQ_FILE="/home/ha5653/backup/YachieLab/ikaken/NK_0137/SEQ_from_Yusuke/editing${editing_rate}/fa_concat_merged/chunksize${Nchunks}/seq_concat_chunksize${Nchunks}.fa.gz"
        ROOT_SEQ_FILE="/home/ha5653/backup/YachieLab/ikaken/NK_0137/SEQ_from_Yusuke/editing${editing_rate}/root_concat/root_concat_${Nchunks}chunks.fa.gz"
        FRACTAL_IN_FILE="/home/ha5653/backup/YachieLab/ikaken/NK_0137/SEQ_from_Yusuke/editing${editing_rate}/fa_concat_merged/chunksize${Nchunks}/seq_concat_chunksize${Nchunks}.root.fa.gz"
        
        cat $ROOT_SEQ_FILE $TIP_SEQ_FILE > $FRACTAL_IN_FILE
        FRACTAL -i $FRACTAL_IN_FILE -m rapidnjNJ -k 10000 -z 1000 -t 10000 -p MP -e -I "-l s_vmem=128G -l mem_req=128G -l ljob" -O "-l s_vmem=16G -l mem_req=16G" -A "-l s_vmem=512G -l mem_req=512G -l lmem" -d 100 -j ${JOB_NAME} -g > fractal${JOB_NAME}.out 2> fractal${JOB_NAME}.err
    
        FRACTAL_TREE="${RESULT_DIR}/${Nchunks}chunks/FRACTALout.nwk"
        PRESUME_TREE="/home/ha5653/backup/YachieLab/ikaken/NK_0137/PRESUMEout/PRESUMEout_combined.nwk"
        COPIED_PRESUME_TREE="${RESULT_DIR}/${Nchunks}chunks/PRESUMEout_combined.nwk"

        cp $PRESUME_TREE $COPIED_PRESUME_TREE

        echo -n "$Nchunks," > result.out
        python3 $PARTIAL $FRACTAL_TREE $COPIED_PRESUME_TREE $TREEDISTR >> result.out

    fi
    i=$(expr $i + 1)
done

