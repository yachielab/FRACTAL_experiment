#!/bin/bash
#$ -S /bin/bash
#$ -N NK_0137
#$ -cwd
#$ -o /home/ha5653/backup/YachieLab/ikaken/NK_0137/out
#$ -e /home/ha5653/backup/YachieLab/ikaken/NK_0137/err

CODEDIR=/home/ha5653/backup/YachieLab/code
NK_0137_dir="/home/ha5653/backup/YachieLab/ikaken/NK_0137/VS_CASSIOPEIA"
TREEDISTR=${CODEDIR}"/FRACTAL_exp/MyCodes/R/NormRFdist.R"
PARTIAL=${CODEDIR}"/FRACTAL_exp/MyCodes/python/PartialNRFdist.py"

source ${HOME}/.bash_profile

conda activate NK_0146_Keito

echo $(hostname)

export GRB_LICENSE_FILE=/home/ha5653/software/gurobi902/license/$(hostname)/gurobi.lic

echo ${GRB_LICENSE_FILE}

i=1

for x in 0.25 0.025 0.0; do
    for Nchunks in 5 10 20 40 100; do
        for Ntips in 100 1000 10000; do
            for mode in subclade; do
                for charmat_filepath in $(ls /home/ha5653/backup/YachieLab/ikaken/NK_0137/VS_CASSIOPEIA/from_yusuke/x${x}/cassiopeia_format/${mode}/${Nchunks}chunks/*tip${Ntips}_*character_matrix.txt | head -n 1); do

                    if [ ${mode} = "subclade" ]; then

                        cladename=$(echo ${charmat_filepath} | cut -d "/" -f14 | cut -d"_" -f2| sed 's/tiplabels//g')
                        true_nwk_path="/home/ha5653/backup/YachieLab/ikaken/NK_0137/VS_CASSIOPEIA/from_yusuke/newick/${mode}/tip${Ntips}_BEmodel${cladename}.nwk"
                    
                    elif [ ${mode} = "random" ]; then

                        seedname=$(echo ${charmat_filepath} | cut -d "/" -f14 | cut -d"_" -f1)
                        true_nwk_path="/home/ha5653/backup/YachieLab/ikaken/NK_0137/VS_CASSIOPEIA/from_yusuke/newick/${mode}/${seedname}_tip${Ntips}_${Ntips}_sampled.nwk"

                    fi

                    
                    for method in greedy hybrid ilp; do 

                        if [ ${method} = "hybrid" ]; then

                            option="--num_threads 1 --max_neighborhood_size 6000 --time_limit 5000 --cutoff 200"

                        else

                            option="--num_threads 1 --time_limit 12600"

                        fi
                        
                        if [ $i -eq ${SGE_TASK_ID} ];then

                            OUTDIR=${NK_0137_dir}/cassiopeia/TASK_ID${SGE_TASK_ID}

                            mkdir ${OUTDIR}; cd ${OUTDIR}


                            #remove redundant sequences

                            (head -n1 ${charmat_filepath}; cat ${charmat_filepath}| awk '{$1=$1";";print $0}' | tr ' ' '\t' | tail -n +2| sort -k2 -t";" | uniq -f1 | tr -d ";") > ${OUTDIR}/charmat.txt
                            charmat_filepath=${OUTDIR}/charmat.txt

                            Nseqs=$(cat ${charmat_filepath} | tail -n +2| wc -l)

                            echo -n "${SGE_TASK_ID},${x},${Nchunks},${Nseqs},${charmat_filepath},${true_nwk_path},${method}," >> ${OUTDIR}/result.out

                            /usr/bin/time -f "%M,KB,%e,sec," timeout 86400 reconstruct-lineage --${method} ${option} ${charmat_filepath} ${OUTDIR}/cassiopeia.pre.nwk 2>&1 1> fractal.out | tr -d "\n" >> ${OUTDIR}/result.out
                            cat ${charmat_filepath} | awk '{c=$1"\t";for(i=2;i<=NF;i++) c=c $i"|"; print c}' | sed -e 's/|$//'|awk '{print $2"\t"$1}' | tail -n +2 > ${OUTDIR}/oldname2newname.txt

                            cat ${OUTDIR}/cassiopeia.pre.nwk | tree_rename -s | tree_rename -n ${OUTDIR}/oldname2newname.txt > ${OUTDIR}/cassiopeia.nwk

                            
                            copied_true_nwk_path=${OUTDIR}/copied_true.nwk

                            cp $true_nwk_path $copied_true_nwk_path

                            python3 $PARTIAL ${OUTDIR}/cassiopeia.nwk ${copied_true_nwk_path} $TREEDISTR | tr -d '\n' >> ${OUTDIR}/result.out

                            echo -n "," >> ${OUTDIR}/result.out

                            tree_compare -s ${OUTDIR}/cassiopeia.nwk.ext -t ${copied_true_nwk_path}.ext --seed 0 >> ${OUTDIR}/result.out

                        fi

                        i=$(expr $i + 1)
                    done
                done
            done
        done
    done
done
