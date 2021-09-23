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
import numpy as np
import gzip
import shutil

def get_option():
    argparser = ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    argparser.add_argument('-sam', '--sam', type=str,required=True,help="Sam file with MD-tag (by samtools calmd)")
    argparser.add_argument('-amplicon_length','--amplicon_length',type=int,default=255,help="Amplicon length of each array")
    argparser.add_argument('-root', '--root', type=str,default="",help="Root fasta, ungzipped (if you need sequence file)")
    argparser.add_argument('-trashbox', '--trashbox', type=str,default="",help="Directory path where delsub tsv file is archived if --no_export_delsub")
    argparser.add_argument('-no_export_seq', '--no_export_seq', action="store_true",help='DO NOT export sequence file')
    argparser.add_argument('-no_export_delsub', '--no_export_delsub', action="store_true",help='DO NOT export delsub per cell file')
    argparser.add_argument('-export_binary_mat','--export_binary_mat',action="store_true",help='export pattern binary matrix')
    argparser.add_argument('-export_insertion_mat','--export_insertion_mat',action="store_true",help='export insertion pattern binary matrix')
    argparser.add_argument('-o', '--outname', type=str,default="gestalt",help='output file name prefix')
    argparser.add_argument('-d', '--outdir', type=str, default=".", help='output directory')
    return argparser.parse_args()

def parseSam(sam):
    # outParsedFile=outdir+outname+".mutationpattern"
    samfile=pysam.AlignmentFile(sam,mode="rb")
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
        
        # parsed_cigar=re.sub("D+$","",parsed_cigar)
        # parsed_seq=re.sub("D+$","",parsed_seq)
        
        mut_temp={"cigar":parsed_cigar,"seq":parsed_seq}
        parsedMutationPattern[query]=mut_temp

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
    if type(extracted_pattern) is not str and np.isnan(extracted_pattern):
        matching_string="=;"*amplicon_length
    else:
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
            # matching_string+="D;"*((amplicon_length-1)-pointer)

    matching_string=re.sub(";$","",matching_string)
    return matching_string

def generate_delsub_seq(delsub,root):
    seq=""
    rootseq=""
    with open(root,mode="rt") as r:
        for line in r:
            if not re.search("^>",line):
                rootseq+=line.replace("\n","")
    delsub=delsub.split(";")
    for pos,nuc in enumerate(delsub):
        if nuc=="=":
            seq+=rootseq[pos]
        elif nuc=="D":
            seq+="-"
        else:
            seq+=nuc
    return seq

def generate_ins_table(insertion_pattern,list_ins):
    insertion_pattern=insertion_pattern.split(";")
    binary_list=[]
    for ref_ins in list_ins:
        if ref_ins in insertion_pattern:
            binary_list.append("1")
        else:
            binary_list.append("0")
    return ";".join(binary_list)

def build_binary_mat(pats,list_mutation):
    res_list=[]
    pats_list=pats.split(";")
    for i in list_mutation:
        if i in pats_list:
            res_list.append("1")
        else:
            res_list.append("0")
    return ";".join(res_list)

if __name__ == "__main__":
    opt=get_option()
    sam=opt.sam
    amplicon_length=opt.amplicon_length
    root=opt.root
    trashbox=opt.trashbox
    no_export_seq=opt.no_export_seq
    no_export_delsub=opt.no_export_delsub
    export_binary_mat=opt.export_binary_mat
    export_insertion_mat=opt.export_insertion_mat
    outdir=opt.outdir
    outname=opt.outname
    outfile=outdir+"/"+outname+".MUTATIONTABLE.tsv.gz"

    print("Generating table...",flush=True)
    parsedMutationPattern=parseSam(sam)
    df_parsed_from_sam=pd.DataFrame.from_dict(parsedMutationPattern,orient="index")

    ser_concat_parsedSam=df_parsed_from_sam["cigar"].str.cat(df_parsed_from_sam["seq"],sep=";")
    out_df=pd.DataFrame(ser_concat_parsedSam.map(generateMutationTable))
    out_df.to_csv(outfile,sep="\t",header=False,compression="gzip")

    #extract reliable mutation
    outfile_ins=re.sub("/$","",outdir)+"/"+outname+"_insertion_per_cell.tsv.gz"
    outfile_delsub=re.sub("/$","",outdir)+"/"+outname+"_delsub_per_cell.tsv.gz"
    outfile_delsub_seq=re.sub("/$","",outdir)+"/"+outname+"_delsub.SEQUENCE.tsv.gz"
    outfile_binary=re.sub("/$","",outdir)+"/"+outname+"_delsub.BINARY.tsv.gz"

    print("Start processing Insertion and Del/Sub...",flush=True)
    out_df=pd.read_csv(outfile,sep="\t",header=None)
    print("Generating delsub table...",flush=True)
    out_df[1]=out_df[1].map(lambda x: generate_delsub(x,amplicon_length=amplicon_length)) 
    print("Merging delsub table...",flush=True)
    delsub_df=pd.concat([out_df[0],out_df[1].str.split(";",expand=True)],axis=1) 
    print("Exporting tables...",flush=True)
    delsub_df.to_csv(outfile_delsub,sep="\t",index=False,header=False,compression="gzip")

    if not no_export_seq:
        print("Exporting sequence file...",flush=True)
        delsub_seq=out_df[1].map(lambda x:generate_delsub_seq(x,root=root))
        delsub_seq_df=pd.concat([out_df[0],delsub_seq],axis=1)
        delsub_seq_df.to_csv(outfile_delsub_seq,sep="\t",header=False,index=False,compression="gzip")

    if no_export_delsub:
        print("Move delsub file to archive space to save the disk space...")
        shutil.move(outfile_delsub,trashbox+"/"+outname+"_delsub_per_cell.tsv.gz")

    if export_insertion_mat:
        list_ins=list(generate_ins_set(out_df[1]))
        ser_ins=out_df[1].map(generate_ins_series)
        print("Generating insertion table...",flush=True)
        list_ins=sorted(list_ins,key=lambda x:float(x.split("_")[1]))
        ser_ins=ser_ins.map(lambda x: generate_ins_table(x,list_ins=list_ins))
        print("Merging inserion table...",flush=True)
        ins_df=ser_ins.str.split(";",expand=True)
        ins_df=pd.concat([out_df[0],ins_df],axis=1)
        ins_df.to_csv(outfile_ins,sep="\t",index=False,header=False,compression="gzip")

    if export_binary_mat:
        print("Building binary pattern matrix...",flush=True)
        df=pd.read_csv(outfile,sep="\t",header=None)
        set_mutation=set()
        for i in df[1]:
            i=set(i.split(";"))
            set_mutation|=i
        list_mutation=list(set_mutation)
        df[1]=df[1].map(lambda x:build_binary_mat(x,list_mutation))
        s=df[1].str.split(";",expand=True)
        df_new=pd.concat([df[0],s],axis=1)
        print("Exporting binary pattern matrix...",flush=True)
        df_new.to_csv(outfile_binary,sep="\t",index=False,header=False,compression="gzip")
