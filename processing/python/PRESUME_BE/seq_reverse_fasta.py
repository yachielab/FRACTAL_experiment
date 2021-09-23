import pandas as pd
import sys
import re

def get_insertion(mutpat):
    mut_list=mutpat.split(";")
    ins_pos_dic=dict()
    for each_mutation in mut_list:
        if "I" in each_mutation:
            pos=each_mutation.split("_")[1]
            sequence=each_mutation.split("_")[2]
            ins_pos_dic[pos]=sequence
    return ins_pos_dic


def insert_insert(seq_mutpat):
    seq=seq_mutpat.split("+")[0]
    mutpat=seq_mutpat.split("+")[1]
    
    if mutpat=="nomutation":
        seq_out=seq
        seq_out=re.sub("-","",seq_out)
        return seq_out

    ins_pos_dic=get_insertion(mutpat)
    seq_out=""
    if "-0.5" in ins_pos_dic:
        seq_out=ins_pos_dic["-0.5"]

    for pos,c in enumerate(seq):
        if str(pos+0.5) in ins_pos_dic:
            seq_out+=c+ins_pos_dic[str(pos+0.5)]
        else:
            seq_out+=c
    seq_out=re.sub("-","",seq_out)
    return seq_out
        

if __name__ == "__main__":    
    delsub=pd.read_csv(sys.argv[1],sep="\t",header=None)
    muttbl=pd.read_csv(sys.argv[2],sep="\t",header=None)
    outfile=sys.argv[3]

    muttbl=muttbl.fillna("nomutation")
    concat_series=delsub[1].str.cat(muttbl[1],sep="+")
    concat_series=concat_series.map(insert_insert)

    out_df=pd.DataFrame({'header':delsub[0],'seq':concat_series})
    out_df.to_csv(outfile,sep="\t",header=False,index=False)