#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from argparse import ArgumentParser
import re
import time
import os
import pandas as pd
import numpy as np

class MYERROR(Exception):
    pass

def get_option():
    argparser = ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    argparser.add_argument('-i', '--inputTable', type=str,required=True,help="Mutation table (concatenated across chunks)")
    argparser.add_argument('-target', '--target_site_per_unit', type=str,default="0-11,12-38,39-65,66-92,93-119,120-146,147-173,174-200,201-227,228-254",help="gRNA target site per unit. Start with 0!")
    argparser.add_argument('-n_chunk','--num_chunk',type=int,default=1,help="Number of chunk contained in the input mutation table")
    argparser.add_argument('-amp_len','--amplicon_length',type=int,default=255,help="amplicon seq length")
    argparser.add_argument('-o', '--outname', type=str,default="",help='output file name prefix')
    argparser.add_argument('-d', '--outdir', type=str, default=".", help='output directory')
    return argparser.parse_args()

def parse_target_list(target):
    target=target.split(",")
    target=[[int(k) for k in i.split('-')] for i in target]
    return target


def get_target_all(target,n_chunk,amp_len):
    target_list_all=[]
    for chunk in range(n_chunk):
        for t in target:
            target_start=t[0]+chunk*amp_len
            target_end=t[1]+chunk*amp_len
            target_list_all.append([target_start,target_end])
    return target_list_all

def parse_mut(each_mut,mutation_list_per_target,flag_deletion,mutated_position,deletion_end,target_start,target_end,index,mutation_split_by_target):
    if not flag_deletion: #substitution or insertion
        mutation_list_per_target.append(each_mut)
    else: #deletion
        if deletion_end <= target_end: #intra-site deletion
            # if mutated_position==target_start and deletion_end==target_end: #complete deletion of single target site
            #     mutation_list_per_target.append("None")
            # else:
            mutation_list_per_target.append(each_mut)

        else: #deletion_end > target_end, which means inter-site deletion!
            # if mutated_position==target_start:
            #     mutation_list_per_target.append("None")
            # else: #in case of inter-site deletion, deletion from the starting position to the end of gRNA target site is taken as a mutation pattern for the current gRNA target site.
            #     mutation_list_per_target.append("_".join(["D",str(mutated_position),str(target_end)]))
            
            mutation_list_per_target.append(each_mut)

            #Done for the first gRNA target site
            mutation_split_by_target.append(";".join(mutation_list_per_target))

            #Go to the next target site
            mutation_list_per_target=[]
            index+=1
            flag_stay=True
            while flag_stay:                        
                target_start=target_list_all[index][0]
                target_end=target_list_all[index][1]
        
                if deletion_end <= target_end: #end of iteration
                    # mutation_list_per_target.append("_".join(["D",str(target_start),str(deletion_end)]))
                    mutation_list_per_target.append(each_mut)
                    flag_stay=False

                # elif deletion_end == target_end: #end of iteration
                #     mutation_list_per_target.append("None")
                #     # mutation_list_per_target=[]
                #     # index+=1
                #     flag_stay=False
                elif deletion_end > target_end: #The deletion is still going over the current target site. continue the iteration
                    mutation_split_by_target.append(each_mut)
                    mutation_list_per_target=[]
                    index+=1
                    flag_stay=True
                else:
                    print("Something wrong!")
                    print(target_start,target_end,each_mut,mutated_position)
                    raise MYERROR()
    return mutation_split_by_target,mutation_list_per_target,index

def parse_mutation_list_2(mutations,target_list_all): #run for a mutation pattern separated by ";" in each cell
    if pd.isna(mutations):
        mutation_split_by_target=("D_0_"+str(target_list_all[len(target_list_all)-1][1]+1))*len(target_list_all)
        return "|".join(mutation_split_by_target)  

    mutations=mutations.split(";")
    mutation_split_by_target=[] #final mutation list split by target
    mutation_list_per_target=[] #intermediate mutation list in one target
    
    index=0
    mutated_position=0
    deletion_end=0
    
    for each_mut in mutations:
        target_start=target_list_all[index][0] #current target site starting position
        target_end=target_list_all[index][1] #current target site ending position
        if re.search("S",each_mut): #substitution
            flag_deletion=False
            mutated_position=int(each_mut.split("_")[1])

        elif re.search("I",each_mut): #insertion
            flag_deletion=False
            mutated_position=float(each_mut.split("_")[1])+0.5

        else: #deletion
            flag_deletion=True
            mutated_position=int(each_mut.split("_")[1])
            deletion_end=int(each_mut.split("_")[2])

        #Assign all mutations into each gRNA target site
        if mutated_position < target_start: #Error treatment
            print("Something wrong!")
            print(target_start,target_end,each_mut,mutated_position,mutations)
            raise MYERROR("Error!"+each_mut)

        elif mutated_position <= target_end: #mutation starting position is before the current gRNA target site ending position
            """
            For the mutation inside the target site, if it is inter-site deletion, 
            mutation pattern will be iterated in the inter target sites.
            """
            mutation_split_by_target,mutation_list_per_target,index=parse_mut(each_mut,mutation_list_per_target,flag_deletion,mutated_position,deletion_end,target_start,target_end,index,mutation_split_by_target)
        
        elif mutated_position > target_end: #mutation starts after the end of the current gRNA target site. (=curresnt target site is done)
            if mutated_position==target_list_all[len(target_list_all)-1][1]+1: #mutation starts at the end of concatenated arrays (=end of process)
                mutation_list_per_target.append(each_mut)
                continue

            if not mutation_list_per_target: #non-edited target site
                mutation_list_per_target.append("None")

            mutation_split_by_target.append(";".join(mutation_list_per_target)) #finalize the mutation pattern in the current gRNA target sote

            #Go to the next gRNA target site
            mutation_list_per_target=[]
            index+=1
            target_start=target_list_all[index][0]
            target_end=target_list_all[index][1]
            if mutated_position <= target_end: #current mutation is in the current target site
                mutation_split_by_target,mutation_list_per_target,index=parse_mut(each_mut,mutation_list_per_target,flag_deletion,mutated_position,deletion_end,target_start,target_end,index,mutation_split_by_target)
            else:
                while True:
                    mutation_split_by_target.append("None")
                    index+=1
                    target_start=target_list_all[index][0]
                    target_end=target_list_all[index][1]
                    if mutated_position <= target_end:
                        break
                mutation_split_by_target,mutation_list_per_target,index=parse_mut(each_mut,mutation_list_per_target,flag_deletion,mutated_position,deletion_end,target_start,target_end,index,mutation_split_by_target)
    
    mutation_split_by_target.append(";".join(mutation_list_per_target))

    if len(target_list_all) > len(mutation_split_by_target):
        for i in range(len(mutation_split_by_target),len(target_list_all)):
            mutation_split_by_target.append("None")
    return "|".join(mutation_split_by_target)        

if __name__ == "__main__":
    opt=get_option()
    inputTable=opt.inputTable
    target=opt.target_site_per_unit
    n_chunk=opt.num_chunk
    amp_len=opt.amplicon_length
    outdir=opt.outdir
    outname=opt.outname

    #Preparing a list of target sites per single unit (255 bps in the GESTALT case)
    target=parse_target_list(target)
    tbl=pd.read_csv(inputTable,sep="\t",header=None)

    print("parsing start")
    #Preparing a list of target sites throughout whole chunks (concatenated)
    target_list_all=get_target_all(target,n_chunk,amp_len)

    #Make all mutations be separated by each gRNA target
    processed=tbl[1].map(lambda x: parse_mutation_list_2(x,target_list_all))
    processed=processed.str.split("|",expand=True)

    nrow=processed.shape[0]
    ncol=processed.shape[1]
    additional_col=[]

    #Number of target sites must be multiples of 3. if not, additional column is added with "None" patterns
    for i in range(3-ncol%3):
        additional_col.append("add"+str(i))
    if (3-ncol%3)>0 and (3-ncol%3)<3:
        for i in additional_col:
            processed[i]=["None" for k in range(nrow)]
    n_int=int(processed.shape[1]/3)
    int_bcs=["int"+str(i+1) for i in range(n_int)]


    processed_arr=np.array(processed)

    #reshape the table to make 3 columns
    processed_arr=processed_arr.reshape(-1,3)
    processed=pd.DataFrame(processed_arr)
    processed["intBC"]=int_bcs*nrow

    cells=[]
    for i in tbl[0]:
        cells+=[i]*n_int

    processed=pd.concat([pd.Series(cells),processed],axis=1)
    processed.columns=["cellBC","r1","r2","r3","intBC"]

    processed.to_csv(outdir+"/"+outname+"_to_cassiopeia.tsv",sep='\t',header=True,index=False)