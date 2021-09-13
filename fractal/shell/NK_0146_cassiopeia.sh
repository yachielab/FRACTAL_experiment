#!/bin/bash
#$ -S /bin/bash
#$ -N NK_0146
#$ -cwd
#$ -o /home/ha5653/backup/YachieLab/ikaken/NK_0146/out
#$ -e /home/ha5653/backup/YachieLab/ikaken/NK_0146/err

CODEDIR=/home/ha5653/backup/YachieLab/code
NK_0146_dir="/home/ha5653/backup/YachieLab/ikaken/NK_0146"
TREEDISTR=${CODEDIR}"/FRACTAL_exp/MyCodes/R/NormRFdist.R"
PARTIAL=${CODEDIR}"/FRACTAL_exp/MyCodes/python/PartialNRFdist.py"

source ${HOME}/.bash_profile

conda activate NK_0146_Keito

echo $(hostname)

export GRB_LICENSE_FILE=/home/ha5653/software/gurobi902/license/$(hostname)/gurobi.lic

echo ${GRB_LICENSE_FILE}

i=1

for rep in $(seq 0 4); do
    for method in greedy hybrid ilp; do 
        for mut_rate in 5 6 7 75 8 85 88 9 92 93 94 95 96 97 98 99; do
            charmat_filename=true_network_one_minus_mutation_rate_0_${mut_rate}_run_${rep}_pkl.charmat.txt
            charmat_filepath=${NK_0146_dir}/dataset/charmat/${charmat_filename}

        

            if [ ${method} = "hybrid" ]; then

                option="--num_threads 1 --max_neighborhood_size 6000 --time_limit 5000 --cutoff 10"

            else

                option="--num_threads 1 --time_limit 12600"

            fi

            if [ $i -eq ${SGE_TASK_ID} ];then

                OUTDIR=${NK_0146_dir}/cassiopeia/TASK_ID${SGE_TASK_ID}

                mkdir ${OUTDIR}

                echo -n "${SGE_TASK_ID},${mut_rate},${rep},${method}," >> ${OUTDIR}/result.out

                /usr/bin/time -f "%M,KB,%e,sec," reconstruct-lineage --${method} ${option} ${charmat_filepath} ${OUTDIR}/cassiopeia.pre.nwk 2>&1 1> ${OUTDIR}/cassiopeia.outerr | tr -d "\n" >> ${OUTDIR}/result.out
                cat ${charmat_filepath} | awk '{c=$1"\t";for(i=2;i<=NF;i++) c=c $i"|"; print c}' | sed -e 's/|$//'|awk '{print $2"\t"$1}' | tail -n +2 > ${OUTDIR}/oldname2newname.txt

                cat ${OUTDIR}/cassiopeia.pre.nwk | tree_rename -s | tree_rename -n ${OUTDIR}/oldname2newname.txt > ${OUTDIR}/cassiopeia.nwk

                true_nwk_path=${NK_0146_dir}/dataset/newick/$(echo "${charmat_filename}" | cut -f1 -d".").tree.nwk

                copied_true_nwk_path=${OUTDIR}/copied_true.nwk

                cp $true_nwk_path $copied_true_nwk_path

                python3 $PARTIAL ${OUTDIR}/cassiopeia.nwk ${copied_true_nwk_path} $TREEDISTR | tr '\n' ','  >> ${OUTDIR}/result.out
                tree_compare -s ${OUTDIR}/cassiopeia.nwk -t ${copied_true_nwk_path}.ext --seed 0 >> ${OUTDIR}/result.out
            fi
        i=$(expr $i + 1)
        done
    done
done