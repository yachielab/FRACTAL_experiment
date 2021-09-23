#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
import argparse
from argparse import ArgumentParser
import re
import glob
import collections
import time
import os
import pysam
import pandas as pd

class parseError(Exception):
    pass

def get_option():
    argparser = ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    argparser.add_argument('-sam', '--sam', type=str,required=True,help="Sam file with MD-tag (by samtools calmd)")
    # argparser.add_argument('-ref', '--needleReference', type=str,required=True,help="Reference fasta file for needle alignment")
    argparser.add_argument('-annot', '--annotation', type=str,required=True,help="Reference annotation which include cut site")
    argparser.add_argument('-idenrate', '--idenrate', type=float,default=0.85,help="Hit ratio in aligned region between query and reference")
    argparser.add_argument('-identity', '--identity', type=int,default=90,help="Minimum number of matches between query and reference")
    argparser.add_argument('-trim5prime', '--trim5prime', type=int,default=19,help="Trim sequence from 5prime end")
    argparser.add_argument('-primer5', '--primer5', type=int,default=17,help="Primer binding site length started from 5' end (after trimmed)")
    argparser.add_argument('-primer5ed', '--primer5ed', type=int,default=2,help="Maximum edit distance allowed at the 5' primer binding site")
    argparser.add_argument('-noRemoveFiles', '--noRemoveFiles', action="store_true",help="Don't remove intermediate files")
    argparser.add_argument('-no_annot', '--no_annot', action="store_true",help="Don't use cut site selection")
    argparser.add_argument('-o', '--outname', type=str,default="gestalt",help='output file name prefix')
    argparser.add_argument('-d', '--outdir', type=str, default=".", help='output directory')
    argparser.add_argument('-extend', '--extend', type=int,default=0,help="Additional N bp upstream of cut site")
    argparser.add_argument('-amplicon_length', '--amplicon_length', type=int,default=255,help="Amplicon size")
    return argparser.parse_args()

def parseSam(sam):
    # outParsedFile=outdir+outname+".mutationpattern"
    samfile=pysam.AlignmentFile(sam,mode="r")
    parsedMutationPattern={}
    for line in samfile:
        mut_temp={}
        query=line.query_name
        cigar=line.cigarstring
        # cigar_tup=line.cigartuples
        # md_tag=line.tags
        seq=line.query_sequence
        
        #parse cigar string
        parsed_cigar=""
        n=""
        for c in cigar:
            if re.search("[0-9]",c):
                n+=c
            else:
                parsed_cigar+=c*int(n)
                n=""
        
        #parse seq
        parsed_seq=""
        d_len=0
        for cnt,c in enumerate(parsed_cigar):
            if c=="D":
                d_len+=1
                parsed_seq+="D"
            else:
                parsed_seq+=seq[cnt-d_len]
        
        parsed_cigar=re.sub("D+$","",parsed_cigar)
        parsed_seq=re.sub("D+$","",parsed_seq)
        
        mut_temp={"cigar":parsed_cigar,"seq":parsed_seq}
        parsedMutationPattern[query]=mut_temp

    return parsedMutationPattern

def alignmentQC(parsedMutationPattern,ref_annot,primer5ed,idenrate,identity):
    cigar_now=parsedMutationPattern.split(";")[0]
    seq_now=parsedMutationPattern.split(";")[1]
    primer5_matching=seq_now[ref_annot["primer_5prime"][0]:ref_annot["primer_5prime"][1]]
    editDist=0
    for c in primer5_matching:
        if not c=="=":
            editDist+=1

    seq_match=0
    seq_ident=0
    for pos,cigar_c in enumerate(cigar_now):
        if cigar_c=="M":
            seq_match+=1
            if seq_now[pos]=="=":
                seq_ident+=1
    seq_idenrate=seq_ident/seq_match
    
    if seq_idenrate<=idenrate or seq_ident<=identity or editDist>primer5ed:
        return None
    else:
        return parsedMutationPattern

def generateMutationTable(parsedMutationPattern):
    detectedPatternList=[]
    delflag=False
    insflag=False
    subflag=False
    n_ins=0
    mut_index=0
    cigar_now=parsedMutationPattern.split(";")[0]
    seq_now=parsedMutationPattern.split(";")[1]
    for cnt,cigstr in enumerate(cigar_now):
        if cigstr=="I":
            n_ins+=1
        pos=cnt-n_ins

        if cigstr=="M":
            delflag=False
            insflag=False
            if not seq_now[cnt]=="=":
                if not subflag:
                    mut_index+=1
                    subflag=True
                    pat_now="_".join(["S",str(pos),seq_now[cnt]])
                    detectedPatternList.append(pat_now)
                else:
                    detectedPatternList[mut_index-1]+=seq_now[cnt]
            else:
                subflag=False

        if cigstr=="I":
            delflag=False
            subflag=False
            if not insflag:
                mut_index+=1
                insflag=True
                pat_now="_".join(["I",str(pos+0.5),seq_now[cnt]])
                detectedPatternList.append(pat_now)
            else:
                detectedPatternList[mut_index-1]+=seq_now[cnt]

        if cigstr=="D":
            insflag=False
            subflag=False
            if not delflag:
                mut_index+=1
                delflag=True
                pat_now="_".join(["D",str(pos),str(pos)])
                detectedPatternList.append(pat_now)
            else:
                detectedPatternList[mut_index-1]=re.sub("[0-9]+$",str(pos),detectedPatternList[mut_index-1])

    for pos,pat in enumerate(detectedPatternList):
        pat_split=pat.split("_")
        if pat_split[0]=="I" or pat_split[0]=="S":
            insertion=pat_split[2]
            if re.search("^N+$",insertion):
                del detectedPatternList[pos]

    return ";".join(detectedPatternList)

def judge_in_cutsite(pos,cutsite):
    for each_cut in cutsite:
        if pos>=each_cut[0] and pos<=each_cut[1]:
            return True
    return False

def judge_include_cutsite(pos1,pos2,cutsite):
    for each_cut in cutsite:
        if pos1<=each_cut[0] and pos2>=each_cut[1]:
            return True
    return False

def return_cutsite_edit(pats,cutsite):
    pats=pats.split(";")
    return_list=[]
    for p in pats:
        p_split=p.split("_")
        if p_split[0]=="I":
            pos=float(p_split[1])
            if judge_in_cutsite(pos,cutsite):
                return_list.append(p)
        elif p_split[0]=="S":
            pos=int(p_split[1])+len(p_split[2])-1
            if judge_in_cutsite(pos,cutsite):
                return_list.append(p)
        elif p_split[0]=="D":
            pos_1=int(p_split[1])
            pos_2=int(p_split[2])
            if judge_in_cutsite(pos_1,cutsite) or judge_in_cutsite(pos_2,cutsite) or judge_include_cutsite(pos_1,pos_2,cutsite):
                return_list.append(p)
    return ";".join(return_list)

def extract_valid_edits(table,cutsite,outfile):
    table_df = pd.read_csv(table,sep="\t",header=None)
    table_df[1]=table_df[1].apply(return_cutsite_edit,cutsite=cutsite)
    table_df.to_csv(outfile,sep="\t",index=False,header=False)
    return table_df

def generate_ins_set(extracted_pattern):
    ins_set=set()
    for row in extracted_pattern:
        pat_split=row.split(";")
        for each_pat in pat_split:
            indelsub=each_pat.split("_")[0]
            if indelsub=="I":
                ins_set.add(each_pat)
    return ins_set

def generate_ins_series(extracted_pattern):
    ins_string=""
    pat_split=extracted_pattern.split(";")
    for each_pat in pat_split:
        indelsub=each_pat.split("_")[0]
        if indelsub=="I":
            ins_string+=";"+each_pat
    ins_string=re.sub("^;","",ins_string)
    return ins_string

def generate_delsub(extracted_pattern,amplicon_length):
    pat_split=extracted_pattern.split(";")
    pointer=-1
    matching_string=""
    for each_pat in pat_split:
        indelsub=each_pat.split("_")[0]
        if indelsub=="I":
            continue
        pos=int(each_pat.split("_")[1])
        additional=each_pat.split("_")[2]

        if indelsub=="D":
            additional=int(additional)
            if additional>=amplicon_length:
                additional=amplicon_length-1
            matching_string+="=;"*(pos-pointer-1)
            matching_string+="D;"*(additional-pos+1)
            pointer=additional
        
        elif indelsub=="S":
            matching_string+="=;"*(pos-pointer-1)
            p=[i for i in additional]
            p=";".join(p)+";"
            matching_string+=p
            pointer=pos+len(additional)-1            
    
    if pointer < amplicon_length-1:
        #if scGST(sequence read length is fixed and end deletion is regarded as failure to reach)
        matching_string+="=;"*((amplicon_length-1)-pointer)

        #if PRESUME(all amplicaon length is covered)
        #matching_string+="D;"*((amplicon_length-1)-pointer)

    matching_string=re.sub(";$","",matching_string)
    return matching_string

def generate_ins_table(insertion_pattern,list_ins):
    insertion_pattern=insertion_pattern.split(";")
    binary_list=[]
    for ref_ins in list_ins:
        if ref_ins in insertion_pattern:
            binary_list.append("1")
        else:
            binary_list.append("0")
    return ";".join(binary_list)

if __name__ == "__main__":
    opt=get_option()
    sam=opt.sam
    annotation=opt.annotation
    idenrate=opt.idenrate
    identity=opt.identity
    trim5prime=opt.trim5prime
    primer5=opt.primer5
    primer5ed=opt.primer5ed
    noremove=opt.noRemoveFiles
    outname=opt.outname
    extend=opt.extend
    amplicon_length=opt.amplicon_length
    no_annot=opt.no_annot
    outdir=re.sub("/$","",opt.outdir)+"/"
    outfile=outdir+"/"+outname+".MUTATIONTABLE.tsv"

    ref_annot={}
    with open(annotation,mode="rt") as r:
        for nrow,line in enumerate(r):
            if nrow>0:
                line=line.split(",")
                ref_annot[line[0]]=list(map(int,line[1:]))


    parsedMutationPattern=parseSam(sam)
    df_parsed_from_sam=pd.DataFrame.from_dict(parsedMutationPattern,orient="index")
    ser_concat_parsedSam=df_parsed_from_sam["cigar"].str.cat(df_parsed_from_sam["seq"],sep=";")
    ser_concat_parsedSam=ser_concat_parsedSam.map(lambda x:alignmentQC(x,ref_annot,primer5ed,idenrate,identity))
    concat_df=pd.DataFrame(ser_concat_parsedSam)
    concat_df=concat_df.dropna(how="any")
    concat_ser=concat_df.iloc[:,0]

    final_df=pd.DataFrame(concat_ser.map(generateMutationTable))
    final_df.to_csv(outfile,sep="\t",header=False)


    #extract reliable mutation
    outfile_reliable=re.sub("/$","",outdir)+"/"+outname+"_reliable_mutation.tsv"
    outfile_ins=re.sub("/$","",outdir)+"/"+outname+"_insertion_per_cell.tsv"
    outfile_delsub=re.sub("/$","",outdir)+"/"+outname+"_delsub_per_cell.tsv"
    if not no_annot:
        cutsite=[]
        with open(annotation,mode="rt") as r:
            for line in r:
                line=line.replace("\n","")
                line_csv=line.split(",")
                if re.search("^cut",line_csv[0]):
                    cutsite.append([int(line_csv[1])-extend,int(line_csv[2])])
        table_extracted_df=extract_valid_edits(outfile,cutsite,outfile_reliable)
    else:
        table_extracted_df=pd.read_csv(outfile,sep="\t",header=None)

    set_ins=generate_ins_set(table_extracted_df[1])
    ser_ins=table_extracted_df[1].apply(generate_ins_series)
    ser_delsub=table_extracted_df[1].apply(generate_delsub,amplicon_length=amplicon_length)
    table_extracted_df[1]=ser_delsub

    len_set=set()
    for i in table_extracted_df[1]:
        len_set.add(len(i.split(";")))
    delsub_df=pd.concat([table_extracted_df[0],table_extracted_df[1].str.split(";",expand=True)],axis=1)

    list_ins=list(set_ins)
    list_ins=sorted(list_ins,key=lambda x:float(x.split("_")[1]))
    ser_ins=ser_ins.apply(generate_ins_table,list_ins=list_ins)
    ins_df=ser_ins.str.split(";",expand=True)
    ins_df=pd.concat([table_extracted_df[0],ins_df],axis=1)

    delsub_df.to_csv(outfile_delsub,sep="\t",index=False,header=False)
    ins_df.to_csv(outfile_ins,sep="\t",index=False,header=False)