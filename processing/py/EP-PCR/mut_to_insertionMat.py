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
    argparser.add_argument('-in_file', '--mutationtable_file', type=str,required=True,help="mutationtable_file")
    argparser.add_argument('-amplicon_length','--amplicon_length',type=int,default=255,help="Amplicon length of each array")
    argparser.add_argument('-o', '--outname', type=str,default="gestalt",help='output file name prefix')
    argparser.add_argument('-d', '--outdir', type=str, default=".", help='output directory')
    return argparser.parse_args()

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
    amplicon_length=opt.amplicon_length
    outdir=opt.outdir
    outname=opt.outname
    outfile_ins = outdir + "/" + outname + "_insertion.tsv.gz"
    mutationtable_file=opt.mutationtable_file

    print("Start processing Insertion and Del/Sub...",flush=True)
    out_df=pd.read_csv(mutationtable_file,sep="\t",header=None)
    # print("Generating delsub table...",flush=True)
    # out_df[1]=out_df[1].map(lambda x: generate_delsub(x,amplicon_length=amplicon_length))
    # print("Merging delsub table...",flush=True)
    # delsub_df=pd.concat([out_df[0],out_df[1].str.split(";",expand=True)],axis=1)
    
    
    list_ins=list(generate_ins_set(out_df[1]))
    ser_ins=out_df[1].map(generate_ins_series)

    print("Generating insertion table...",flush=True)
    list_ins=sorted(list_ins,key=lambda x:float(x.split("_")[1]))
    ser_ins=ser_ins.map(lambda x: generate_ins_table(x,list_ins=list_ins))
    print("Merging inserion table...",flush=True)
    ins_df=ser_ins.str.split(";",expand=True)
    ins_df=pd.concat([out_df[0],ins_df],axis=1)
    ins_df.to_csv(outfile_ins,sep="\t",index=False,header=False,compression="gzip")
    