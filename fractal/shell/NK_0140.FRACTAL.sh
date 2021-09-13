#!/bin/bash
#$ -S /bin/bash
#$ -N NK_0140
#$ -cwd
#$ -o /home/ha5653/backup/YachieLab/ikaken/NK_0140/out
#$ -e /home/ha5653/backup/YachieLab/ikaken/NK_0140/err


source ${HOME}/.bash_profile

CODEDIR=/home/ha5653/backup/YachieLab/code

RESULTDIR=/home/ha5653/backup/YachieLab/ikaken/NK_0140/result
TREEDISTR=${CODEDIR}"/FRACTAL_exp/MyCodes/R/NormRFdist.R"
PARTIAL=${CODEDIR}"/FRACTAL_exp/MyCodes/python/PartialNRFdist.py"

del_pos_count=/home/ha5653/backup/YachieLab/ikaken/NK_0140/param/prob.del.0.34.txt
ins_pos_count=/home/ha5653/backup/YachieLab/ikaken/NK_0140/param/prob.ins.0.015.txt
del_len_count=/home/ha5653/backup/YachieLab/ikaken/NK_0120/count/len_count.del.txt
ins_len_count=/home/ha5653/backup/YachieLab/ikaken/NK_0120/count/len_count.ins.txt
initial_fasta=${RESULTDIR}/initial.fasta

i=1
for scale in 0.001 0.01 0.05 $(seq 0.1 0.1 1) $(seq 2 1 10); do 

    for No in 5 10 20 40 100; do

        for t in 1000 1000000; do 

            if [ $i -eq ${SGE_TASK_ID} ]; then

                OUTDIR=${RESULTDIR}/scale${scale}_Nchunks${No}_t${t}

                mkdir ${OUTDIR}
                cd ${OUTDIR}

                cp /home/ha5653/backup/YachieLab/ikaken/NK_0138/PRESUMEout_topology/PRESUMEout.nwk ${OUTDIR}/truth.nwk

                echo -n "${scale},${No},${t}," > ${OUTDIR}/result.out

                seqkit concat $(for i in $(seq 1 ${No}); do echo -n "/home/ha5653/backup/YachieLab/ikaken/NK_0140/dataset/scale${scale}_No${i}/PRESUMEout/root.fa.gz ";done)                                | gzip >  ${OUTDIR}/scale${scale}_Nchunks${No}.fa.gz
                seqkit concat $(for i in $(seq 1 ${No}); do echo -n "/home/ha5653/backup/YachieLab/ikaken/NK_0140/dataset/scale${scale}_No${i}/alignment/scale${scale}_No${i}_delsub.SEQUENCE.fa.gz ";done) | gzip >> ${OUTDIR}/scale${scale}_Nchunks${No}.fa.gz

                /usr/bin/time -f "%M,KB,%e,sec," FRACTAL -i ${OUTDIR}/scale${scale}_Nchunks${No}.fa.gz -k 1000 -t ${t} -m rapidnjNJ -p MP -e 2>&1  1> /dev/null | tr -d "\n" >> ${OUTDIR}/result.out

                python3 $PARTIAL ${OUTDIR}/truth.nwk ${OUTDIR}/FRACTALout.nwk $TREEDISTR >>${OUTDIR}/result.out

            fi

            i=$(expr $i + 1)
        
        done

    done

done