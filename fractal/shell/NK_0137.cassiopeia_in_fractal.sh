FASTA_file_path=$1


i=$(basename ${FASTA_file_path} | sed 's/RENAMED_//g' | sed 's/.edit.fa.aligned//g') 

DIR=$(dirname $FASTA_file_path)

if [ ${method} = "hybrid" ]; then

    option="--num_threads 1 --time_limit 60 --cutoff 4"

else

    option="--num_threads 1 --time_limit 60"

fi



if [ -e ${DIR}/SUBSAMPLE.edit ]; then

    if [ ! -e ${DIR}/charmatrix.reform.sample.$(expr ${i} - 1).txt ]; then
        echo "" > ${DIR}/charmatrix.reform.sample.$(expr ${i} - 1).txt
    fi

    paste <(cat ${DIR}/SUBSAMPLE.edit | cut -f1) <(seqkit fx2tab $FASTA_file_path | cut -f1) > ${DIR}/name2renamed.${i}.txt # name > sXXX
    
    for name in $(cat ${DIR}/SUBSAMPLE.edit | cut -f1); do cat $CHARMATRIX_path ${DIR}/charmatrix.reform.sample.$(expr ${i} - 1).txt | awk -v name=${name} {'if($1==name){print $0}'};done > ${DIR}/charmatrix.reform.sample.txt

    cat ${CHARMATRIX_path}| head -n 1              >  ${DIR}/charmatrix.reform.sample.renamed.txt
    csvjoin -t -c1,1 --no-header-row ${DIR}/name2renamed.${i}.txt ${DIR}/charmatrix.reform.sample.txt | tail -n +2 | tr ',' '\t' | cut -f1 --complement >> ${DIR}/charmatrix.reform.sample.renamed.txt
    

    reconstruct-lineage --${method} ${option} ${DIR}/charmatrix.reform.sample.renamed.txt ${DIR}/cassiopeia.nwk 
    cat ${DIR}/charmatrix.reform.sample.renamed.txt | awk '{c=$1"\t";for(i=2;i<=NF;i++) c=c $i"|"; print c}' | sed -e 's/|$//'|awk '{print $2"\t"$1}' | tail -n +2 |sed 's/cell//g' > ${DIR}/oldname2newname.txt

    cat ${DIR}/cassiopeia.nwk | tree_rename -s | tree_rename -n ${DIR}/oldname2newname.txt | python3 /home/ha5653/backup/YachieLab/code/FRACTAL_exp/MyCodes/python/add_outgroup.py s0| tree_bifurcation --seed 0 --rooted > ${FASTA_file_path}.tree

else 
    (cat ${CHARMATRIX_path}| head -n 1; for name in $(seqkit fx2tab $FASTA_file_path | cut -f1); do cat $CHARMATRIX_path | awk -v name=${name} {'if($1==name){print $0}'};done)> ${DIR}/charmatrix.reform.sample.txt

    reconstruct-lineage --${method} ${option} ${DIR}/charmatrix.reform.sample.txt ${DIR}/cassiopeia.nwk 
    cat ${DIR}/charmatrix.reform.sample.txt | awk '{c=$1"\t";for(i=2;i<=NF;i++) c=c $i"|"; print c}' | sed -e 's/|$//'|awk '{print $2"\t"$1}' | tail -n +2 |sed 's/cell//g' > ${DIR}/oldname2newname.txt

    cat ${DIR}/cassiopeia.nwk | tree_rename -s | tree_rename -n ${DIR}/oldname2newname.txt | python3 /home/ha5653/backup/YachieLab/code/FRACTAL_exp/MyCodes/python/add_outgroup.py root| tree_bifurcation --seed 0 --rooted > ${FASTA_file_path}.tree

fi

cp ${FASTA_file_path}.tree ${DIR}/../PARAM/RAxML_result.PARAM_${i}

mv ${DIR}/charmatrix.reform.sample.renamed.txt ${DIR}/charmatrix.reform.sample.${i}.txt