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

SUBSAMPLE=100
THRESHOLD=100

for rep in $(seq 0 4); do
    for method in greedy hybrid ilp; do 
        for mut_rate in 5 6 7 75 8 85 88 9 92 93 94 95 96 97 98 99; do
            charmat_filename=true_network_one_minus_mutation_rate_0_${mut_rate}_run_${rep}_pkl.charmat.txt
            charmat_filepath=${NK_0146_dir}/dataset/charmat/${charmat_filename}
            editpattern_filepath=${NK_0146_dir}/dataset/editpattern/${charmat_filename}.edit

            if [ $i -eq ${SGE_TASK_ID} ];then

                OUTDIR=${NK_0146_dir}/fractal_cassiopeia/TASK_ID${SGE_TASK_ID}

                mkdir ${OUTDIR}; cd ${OUTDIR}

                echo "CHARMATRIX_path=${charmat_filepath}" >  cassiopeia_in_fractal.sh
                echo "method=${method}"                    >> cassiopeia_in_fractal.sh

                echo -n "${SGE_TASK_ID},${mut_rate},${rep},${method}," >> ${OUTDIR}/result.out
                
                cat ${CODEDIR}/FRACTAL_exp/MyCodes/shell/NK_0146.cassiopeia_in_fractal.sh >> cassiopeia_in_fractal.sh

                export PATH=$(pwd)/:${PATH}

                chmod u+x cassiopeia_in_fractal.sh

                /usr/bin/time -f "%M,KB,%e,sec," FRACTAL -i ${editpattern_filepath} -E -p MP -k ${SUBSAMPLE} -t ${THRESHOLD} -e -s cassiopeia_in_fractal.sh 2>&1 1> fractal.outerr | tr -d "\n" >> result.out

                true_nwk_path=${NK_0146_dir}/dataset/newick/$(echo "${charmat_filename}" | cut -f1 -d".").tree.nwk

                copied_true_nwk_path=${OUTDIR}/copied_true.nwk

                cp $true_nwk_path $copied_true_nwk_path

                python3 $PARTIAL ${OUTDIR}/FRACTALout.nwk ${copied_true_nwk_path} $TREEDISTR | tr '\n' ',' >> ${OUTDIR}/result.out

                tree_compare -s ${OUTDIR}/FRACTALout.nwk.ext -t ${copied_true_nwk_path}.ext --seed 0 >> ${OUTDIR}/result.out
            fi
        i=$(expr $i + 1)
        done
    done
done