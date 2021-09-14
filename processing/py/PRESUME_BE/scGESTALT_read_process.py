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

class CMDERROR(Exception):
    pass

def get_option():
    argparser = ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    argparser.add_argument('-i', '--inputFasta', type=str,required=True,help="fasta file splitted by cell barcode")
    argparser.add_argument('-cb', '--cellbarcode', type=str,required=True,help="cell barcode of this file")
    argparser.add_argument('-header_pat','--header_pat',type=str,default="^.+:(?P<umi>[ATGCN]+):.+$",help="regular expression pattern which indicates cell barcode as 'umi'.\nDefault: ^.+:(?P<umi>[ATGCN]+):.+$ (for scGESTALT)")
    argparser.add_argument('-n_umi','--num_umi',type=int,default=1,help="analyze cell which has at least this number of barcode molecules")
    argparser.add_argument('-dup_umi', '--duplicatedUMI', type=float,default=-1,help="For each UMI, same barocde sequence appeared more than this ratio is reagarded\nas the consensus sequence for the UMI. If < 0, this precedure is not performed\nand CD-hit will be directly carried out.")
    argparser.add_argument('-clst_umi', '--clusterUMI', type=float,default=0.75,help="For each UMI, a cluster which dominates more than this ratio is selected to be\n'true' UMI sequences.")
    argparser.add_argument('-clst_cb', '--clusterCB', type=float,default=0.75,help="For each UMI, same barocde sequence appeared more than this ratio is reagarded\nas the consensus sequence for the UMI: max=1")
    argparser.add_argument('-ref', '--needleReference', type=str,required=True,help="Reference fasta file for needle alignment")
    argparser.add_argument('-annot', '--annotation', type=str,required=True,help="Reference annotation which include cut site")
    argparser.add_argument('-idenrate', '--idenrate', type=float,default=0.85,help="Hit ratio in aligned region between query and reference")
    argparser.add_argument('-identity', '--identity', type=int,default=90,help="Minimum number of matches between query and reference")
    argparser.add_argument('-trim5prime', '--trim5prime', type=int,default=19,help="Trim sequence from 5prime end")
    argparser.add_argument('-primer5', '--primer5', type=int,default=17,help="Primer binding site length started from 5' end (after trimmed)")
    argparser.add_argument('-primer5ed', '--primer5ed', type=int,default=2,help="Maximum edit distance allowed at the 5' primer binding site")
    argparser.add_argument('-noRemoveFiles', '--noRemoveFiles', action="store_true",help="Don't remove intermediate files")
    argparser.add_argument('-o', '--outname', type=str,default="gestalt",help='output file name prefix')
    argparser.add_argument('-d', '--outdir', type=str, default=".", help='output directory')
    return argparser.parse_args()

def is_valid_file(fpath):  
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0

def importFasta(inputFasta):
    fasta={}
    with open(inputFasta,mode="rt") as r:
        for line in r:
            line=line.replace("\n","")
            if re.search("^>",line):
                header=line
            else:
                if header in fasta:
                    fasta[header]+=line
                else:
                    fasta[header]=line
    return fasta

def input_preprocess(input_fa_path,header_pat,n_umi,trim5prime,outname,outdir):
    passed=True
    umi_list=[]
    fasta={}
    with open(input_fa_path,mode="rt") as r:
        for line in r:
            line=line.replace("\n","")
            if re.search("^>",line):
                umi_match=re.search(header_pat,line)
                umi_now=umi_match.groupdict()["umi"]
                if umi_now not in umi_list:
                    fasta[umi_now]=collections.defaultdict(str)
                umi_list.append(umi_now)

                header=line.replace(" ","|")   
                flag=True             
            else:
                if flag:
                    line=line[trim5prime:]
                    flag=False
                fasta[umi_now][header]+=line
    
    if len(set(umi_list))<n_umi:
        log="number of UMI: "+str(len(set(umi_list)))+"\n"+"thresholding for depth: Dropped"+"\n"
        passed=False
        return log,passed
    
    for umi in fasta:
        with open(outdir+outname+".umi."+umi+".fa",mode="wt") as w:
            for header in fasta[umi]:
                w.write(header+"\n")
                w.write(fasta[umi][header]+"\n")
    
    log="number of UMI: "+str(len(set(umi_list)))+"\n"+"thresholding for depth: Passed"+"\n"
    return log,passed

def fastaDedup(inputFasta_paths,outname,outdir):
    fasta=dict()
    for fa_path in inputFasta_paths:
        fasta_tmp=importFasta(fa_path)
        fasta.update(fasta_tmp)
    with open(outdir+outname+".dedup.fa",mode="wt") as w:
        for key in fasta:
            w.write(key+"\n")
            w.write(fasta[key]+"\n")
    dedupFastaName=outdir+outname+".dedup.fa"
    return dedupFastaName

def findConsensusSeq_molecule(inputFasta,duplicatedUMI,umi_now,clusterUMI,outname,outdir):
    outfile_cons=re.sub(outname+"\.umi\.[ATGCN]+\.fa",outname+".consensus.umi."+umi_now+".fa",inputFasta)
    outfile_cdhit=re.sub(outname+"\.umi\.[ATGCN]+\.fa",outname+".cdhit.umi."+umi_now+".fa",inputFasta)
    outfile_cdhit_largestCluster=re.sub(outname+"\.umi\.[ATGCN]+\.fa",outname+".cdhit.largestCluster.umi."+umi_now+".fa",inputFasta)
    outfile_mafft=re.sub(outname+"\.umi\.[ATGCN]+\.fa",outname+".cdhit.umi."+umi_now+".mafft",inputFasta)
    outfile_mafft_pre=re.sub(outname+"\.umi\.[ATGCN]+\.fa",outname+".cdhit.umi."+umi_now+".mafft.pre",inputFasta)
    outfile_cons=re.sub(outname+"\.umi\.[ATGCN]+\.fa",outname+".cdhit.umi."+umi_now+".cons",inputFasta)
    t0=time.time()

    if not is_valid_file(inputFasta):
        log="Empty file\nDropped\n"
        passed=False
        return log,passed

    if duplicatedUMI>=0:
        fasta=importFasta(inputFasta)
        c=collections.Counter(list(fasta.values()))
        count_key=list(c.keys())
        count_num=list(c.values())
        if max(count_num)/sum(count_num) >= duplicatedUMI:
            passed=True
            consensusSeq=count_key[count_num.index(max(count_num))]
            with open(outfile_cons,mode="wt") as w:
                w.write(">consensus_"+umi_now+"\n"+consensusSeq+"\n")
        else:
            passed=False
        log="Dominat sequence ratio: "+str(max(count_num)/sum(count_num))+"\n"
        return log,passed
    else:
        #CD-HIT
        cd_hit_cmd_list=["cd-hit-est","-d","0","-c","0.9","-n","7","-i",inputFasta,"-o",outfile_cdhit]
        cd_hit_cmd=" ".join(cd_hit_cmd_list)
        print("Start running CD-HIT for UMI=",umi_now,"/ command:",cd_hit_cmd,flush=True)
        o=subprocess.run(cd_hit_cmd,shell=True,capture_output=True)
        if not o.returncode==0:
            print(o.stderr)
            print(o.stdout)
            raise CMDERROR("CDHIT command error!")
        t_cdhit=time.time()
        print("CD-HIT completed for UMI=",umi_now,"/ computing time:",round(t_cdhit-t0),"sec",flush=True)
        cdhitCluster=outfile_cdhit+".clstr"
        log,passed=cdhitClusterQualityCheck(cdhitCluster,inputFasta,clusterUMI,outfile_cdhit_largestCluster)

        if not passed:
            return log,passed
        
        #MAFFT
        mafft_cmd_list=["mafft",outfile_cdhit_largestCluster,">",outfile_mafft_pre]
        mafft_cmd=" ".join(mafft_cmd_list)
        print("Start running MAFFT for UMI=",umi_now,"/ command:",mafft_cmd,flush=True)
        o=subprocess.run(mafft_cmd,shell=True,capture_output=True)
        if not o.returncode==0:
            print(o.stderr)
            print(o.stdout)
            raise CMDERROR("MAFFT command error!")
        t_mafft=time.time()
        print("MAFFT completed for UMI=",umi_now,"/ computing time:",round(t_mafft-t_cdhit),"sec",flush=True)
        
        #Treat delitions for alignment result
        fa_mafft=importFasta(outfile_mafft_pre)
        fa_seq=list(fa_mafft.values())
        survive=[]
        for c,tup in enumerate(zip(*fa_seq)):
            if "-" in tup:
                num_bar=0
                for i in tup:
                    if i=="-":
                        num_bar+=1
                if num_bar/len(tup)<=0.5:
                    survive.append(c)
            else:
                survive.append(c)
        for key in fa_mafft:
            fa_mafft[key]="".join([fa_mafft[key][pos] for pos in survive])
        num_record=len(fa_mafft)
        with open(outfile_mafft,mode="wt") as w:
            for k in fa_mafft:
                w.write(k+"\n")
                w.write(fa_mafft[k]+"\n")
        
        #CONS
        if num_record==1:
            print("Number of read for UMI=",umi_now,"is 1 and returned the raw sequence.")
            fa_raw=importFasta(outfile_mafft)
            fa_raw={">consensus_"+umi_now:fa_raw[k] for k in fa_raw}
            with open(outfile_cons,mode="wt") as w:
                for k in fa_raw:
                    w.write(k+"\n")
                    w.write(fa_raw[k]+"\n")
        else:
            cons_cmd_list=["cons",outfile_mafft,outfile_cons,"-name","consensus_"+umi_now]
            cons_cmd=" ".join(cons_cmd_list)
            print("Start running CONS for UMI=",umi_now,"/ command:",cons_cmd,flush=True)
            o=subprocess.run(cons_cmd,shell=True,capture_output=True)
            if not o.returncode==0:
                print(o.stderr)
                print(o.stdout)
                raise CMDERROR("CONS command error!")
            t_cons=time.time()
            print("CONS completed for UMI=",umi_now,"/ computing time:",round(t_cons-t_mafft),"sec",flush=True)
        
        print("Finding consensus sequence for UMI=",umi_now,"has been completed!\n")
        return log,passed

def findConsensusSeq_cell(dedupFasta,clusterCB,outname,outdir):
    t0=time.time()
    outfile_cdhit=outdir+outname+".cdhit.cb.fa"
    outfile_cdhit_largestCluster=outdir+outname+".cdhit.largestCluster.cb.fa"
    outfile_mafft=outdir+outname+".cdhit.cb.mafft"
    outfile_mafft_pre=outdir+outname+".cdhit.cb.pre.mafft"
    outfile_cons=outdir+outname+".RESULT.CONSENSUS.fa"

    if not is_valid_file(dedupFasta):
        log="Empty file\nDropped\n"
        passed=False
        return log,passed

    #CD-HIT
    cd_hit_cmd_list=["cd-hit-est","-d","0","-c","0.9","-n","7","-mask","N","-i",dedupFasta,"-o",outfile_cdhit]
    cd_hit_cmd=" ".join(cd_hit_cmd_list)
    print("Start running CD-HIT / command:",cd_hit_cmd,flush=True)
    o=subprocess.run(cd_hit_cmd,shell=True,capture_output=True)
    if not o.returncode==0:
        print(o.stderr)
        print(o.stdout)
        raise CMDERROR("CDHIT command error!")
    t_cdhit=time.time()
    print("CD-HIT completed / computing time:",round(t_cdhit-t0),"sec",flush=True)
    cdhitCluster=outfile_cdhit+".clstr"
    log,passed=cdhitClusterQualityCheck(cdhitCluster,dedupFasta,clusterCB,outfile_cdhit_largestCluster)

    if not passed:
        return log,passed
    
    #MAFFT
    mafft_cmd_list=["mafft",outfile_cdhit_largestCluster,">",outfile_mafft_pre]
    mafft_cmd=" ".join(mafft_cmd_list)
    print("Start running MAFFT / command:",mafft_cmd,flush=True)
    o=subprocess.run(mafft_cmd,shell=True,capture_output=True)
    if not o.returncode==0:
        print(o.stderr)
        print(o.stdout)
        raise CMDERROR("MAFFT command error!")
    t_mafft=time.time()
    print("MAFFT completed / computing time:",round(t_mafft-t_cdhit),"sec",flush=True)
    
    #Treat delitions for alignment result
    fa_mafft=importFasta(outfile_mafft_pre)
    fa_seq=list(fa_mafft.values())
    survive=[]
    for c,tup in enumerate(zip(*fa_seq)):
        if "-" in tup:
            num_bar=0
            for i in tup:
                if i=="-":
                    num_bar+=1
            if num_bar/len(tup)<=0.5:
                survive.append(c)
        else:
            survive.append(c)
    for key in fa_mafft:
        fa_mafft[key]="".join([fa_mafft[key][pos] for pos in survive])
    num_record=len(fa_mafft)
    with open(outfile_mafft,mode="wt") as w:
        for k in fa_mafft:
            w.write(k+"\n")
            w.write(fa_mafft[k]+"\n")
            
    #CONS
    if num_record==1:
        print("Number of read is 1 and returned the raw sequence.")
        cp_cmd_list=["cp",outfile_mafft,outfile_cons]
        cp_cmd=" ".join(cp_cmd_list)
        o=subprocess.run(cp_cmd,shell=True,capture_output=True)
        if not o.returncode==0:
            print(o.stderr)
            print(o.stdout)
            raise CMDERROR("CP command error!")
    else:
        cons_cmd_list=["cons",outfile_mafft,outfile_cons,"-name","consensus_"+umi_now]
        cons_cmd=" ".join(cons_cmd_list)
        print("Start running CONS / command:",cons_cmd,flush=True)
        o=subprocess.run(cons_cmd,shell=True,capture_output=True)
        if not o.returncode==0:
            print(o.stderr)
            print(o.stdout)
            raise CMDERROR("CONS command error!")
        t_cons=time.time()
        print("CONS completed / computing time:",round(t_cons-t_mafft),"sec",flush=True)
    
    print("Finding consensus sequence has been completed!\n")
    return log,passed


def cdhitClusterQualityCheck(cdhitCluster,inputFasta,thresh,largestClusterOutName):
    fasta=importFasta(inputFasta)
    cdhit_cluster_headerpat=r"^.+nt, (?P<fasta_header>\>.+)\.\.\. .+$"
    cluster_numDict={}
    cluster_headerDict={}
    with open(cdhitCluster,mode="rt") as r:
        for line in r:
            line=line.replace("\n","")
            if re.search("^>",line):
                cluster_now=line
                cluster_numDict[cluster_now]=0
                cluster_headerDict[cluster_now]=[]
            else:
                m=re.search(cdhit_cluster_headerpat,line)
                
                extractedHeader=m.groupdict()["fasta_header"]
                cluster_numDict[cluster_now]+=1
                cluster_headerDict[cluster_now].append(extractedHeader)
    
    cluster_dominant_ratio=max(cluster_numDict.values())/sum(cluster_numDict.values())
    print("The largest cluster dominance ratio:",cluster_dominant_ratio,flush=True)
    if cluster_dominant_ratio>=thresh:
        largest_cluster_position=list(cluster_numDict.values()).index(max(cluster_numDict.values()))
        largest_cluster=list(cluster_numDict.keys())[largest_cluster_position]
        largest_cluster_headers=cluster_headerDict[largest_cluster]
        largest_cluster_fasta={k:fasta[k] for k in fasta if k in largest_cluster_headers}

        with open(largestClusterOutName,mode="wt") as w:
            for key in largest_cluster_fasta:
                w.write(key+"\n"+largest_cluster_fasta[key]+"\n")
        passed=True
        log="The largest cluster dominance ratio: "+str(cluster_dominant_ratio)+"\n"+"Quality filtering: Passed\n"
    else:
        passed=False
        log="The largest cluster dominance ratio: "+str(cluster_dominant_ratio)+"\n"+"Quality filtering: Dropped\n"
    
    return log,passed

def needleAlign(inputFasta,ref,outname,outdir):
    t0=time.time()
    outSamName=outdir+outname+".needle.sam"
    needle_cmd_list=["needleall","-aformat3","sam","-gapextend","0.5","-gapopen","10.0","-awidth3","5000","-asequence",ref,"-bsequence",inputFasta,"-outfile",outSamName]
    needle_cmd=" ".join(needle_cmd_list)
    print("Start running NEEDLE / command:",needle_cmd,flush=True)
    o=subprocess.run(needle_cmd,shell=True,capture_output=True)
    if not o.returncode==0:
        print(o.stderr)
        print(o.stdout)
        raise CMDERROR("NEEDLE command error!")
    t_needle=time.time()
    print("NEEDLE completed / computing time:",round(t_needle-t0),"sec\n",flush=True)
    return outSamName

def fixSam(needle_result,ref,outname,outdir):
    outSamName=outdir+outname+".needle.edit.sam"
    outMDSamName=outdir+outname+".MD.sam"
    sam=[]
    ref_fa=importFasta(ref)
    fa_length=len(list(ref_fa.values())[0])
    with open(needle_result,mode="rt") as r:
        for line in r:
            if re.search("^@PG",line):
                line="@SQ\tSN:reference\tLN:"+str(fa_length)+"\n"+line
            sam.append(line)
    with open(outSamName,mode="wt") as w:
        for line in sam:
            w.write(line)
    samtools_cmd_list=["samtools","calmd","-e",outSamName,ref,">",outMDSamName]
    samtools_cmd=" ".join(samtools_cmd_list)
    print("Start running SAMTOOLS CALMD / command:",samtools_cmd,flush=True)
    o=subprocess.run(samtools_cmd,shell=True,capture_output=True)
    if not o.returncode==0:
        print(o.stderr)
        print(o.stdout)
        raise CMDERROR("SAMTOOLS CALMD command error!")
    return outMDSamName

def parseSam(needle_edit_result,idenrate,identity,outname,outdir):
    # outParsedFile=outdir+outname+".mutationpattern"
    samfile=pysam.AlignmentFile(needle_edit_result,mode="r")
    for line in samfile:
        cigar=line.cigarstring
        cigar_tup=line.cigartuples
        md_tag=line.tags
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
    # print(parsed_cigar)
    # print(parsed_seq)
    parsedMutationPattern={"cigar":parsed_cigar,"seq":parsed_seq}
    return parsedMutationPattern
    
def alignmentQC(parsedMutationPattern,ref_annot,primer5ed,idenrate,identity,outname,outdir):
    passed=True
    primer5_matching=parsedMutationPattern["seq"][ref_annot["primer_5prime"][0]:ref_annot["primer_5prime"][1]]
    editDist=0
    for c in primer5_matching:
        if not c=="=":
            editDist+=1
    seq_match=0
    seq_ident=0
    for pos,cigar_c in enumerate(parsedMutationPattern["cigar"]):
        if cigar_c=="M":
            seq_match+=1
            if parsedMutationPattern["seq"][pos]=="=":
                seq_ident+=1
    seq_idenrate=seq_ident/seq_match
    if seq_idenrate<=idenrate or seq_ident<=identity or editDist>primer5ed:
        passed=False
        log="Number of mismatch at the 5' primer site: "+str(editDist)+"\n"+"Query matching rate: "+str(seq_idenrate)+"\n"+"Number of query match: "+str(seq_ident)+"\n"+"Quality filtering: Dropped\n"
    else:
        log="Number of mismatch at the 5' primer site: "+str(editDist)+"\n"+"Query matching rate: "+str(seq_idenrate)+"\n"+"Number of query match: "+str(seq_ident)+"\n"+"Quality filtering: Passed\n"

    return log,passed

def generateMutationTable(parsedMutationPattern,ref_annot):
    detectedPatternList=[]
    delflag=False
    insflag=False
    subflag=False
    n_ins=0
    mut_index=0
    for cnt,cigstr in enumerate(parsedMutationPattern["cigar"]):
        if cigstr=="I":
            n_ins+=1
        pos=cnt-n_ins

        if cigstr=="M":
            delflag=False
            insflag=False
            if not parsedMutationPattern["seq"][cnt]=="=":
                if not subflag:
                    mut_index+=1
                    subflag=True
                    pat_now="_".join(["S",str(pos),parsedMutationPattern["seq"][cnt]])
                    detectedPatternList.append(pat_now)
                else:
                    detectedPatternList[mut_index-1]+=parsedMutationPattern["seq"][cnt]
            else:
                subflag=False

        if cigstr=="I":
            delflag=False
            subflag=False
            if not insflag:
                mut_index+=1
                insflag=True
                pat_now="_".join(["I",str(pos+0.5),parsedMutationPattern["seq"][cnt]])
                detectedPatternList.append(pat_now)
            else:
                detectedPatternList[mut_index-1]+=parsedMutationPattern["seq"][cnt]

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
        

if __name__ == "__main__":
    opt=get_option()
    inputFilepath=opt.inputFasta
    cellbarcode=opt.cellbarcode
    header_pat=opt.header_pat
    n_umi=opt.num_umi
    duplicatedUMI=opt.duplicatedUMI
    clusterUMI=opt.clusterUMI
    clusterCB=opt.clusterCB
    needleReference=opt.needleReference
    annotation=opt.annotation
    idenrate=opt.idenrate
    identity=opt.identity
    trim5prime=opt.trim5prime
    primer5=opt.primer5
    primer5ed=opt.primer5ed
    noremove=opt.noRemoveFiles
    outname=opt.outname+"_"+opt.cellbarcode
    outdir=re.sub("/$","",opt.outdir)+"/"

    #preprocess input fasta, count number of UMIs in this cell and split fasta into several chunks based on UMI sequence
    print("start preprocessing...")
    log,passed=input_preprocess(inputFilepath,header_pat,n_umi,trim5prime,outname,outdir)
    log="#Depth check\n"+log

    passed_cell=False
    passed_alignmentQC=False

    if passed:
        #Decide the consensus GESTALT barcode sequence for each UMI. 
        umi_splitted_filepath=glob.glob(outdir+outname+".umi.*.fa")
        cnt_highQualUMI=0
        for cnt,fa in enumerate(umi_splitted_filepath):
            umi_extraction_pattern=r"^.+\.umi\.(?P<umi_now>[^\.]+).fa$"
            m=re.search(umi_extraction_pattern,fa)
            umi_now=m.groupdict()["umi_now"]

            log_findCons_umi,passed=findConsensusSeq_molecule(fa,duplicatedUMI,umi_now,clusterUMI,outname,outdir)
            print(log_findCons_umi+"\n")

            if passed:
                cnt_highQualUMI+=1

        log+="\n#Deciding consensus seq for each UMI\nnumber of high quality UMI: "+str(cnt_highQualUMI)+"\nnumber of total UMI: "+str(cnt+1)+"\n"
        
        #Dedup the GESTALT barcode sequences for all UMIs in this cell.
        filepaths_each_umi=glob.glob(outdir+outname+".cdhit.umi.*.cons")
        dedupFastaName=fastaDedup(filepaths_each_umi,outname,outdir)

        #Decide the consensus GESTALT barcode sequence for this cell.
        log_findCons_cell,passed_cell=findConsensusSeq_cell(dedupFastaName,clusterCB,outname,outdir)
        log+="\n#Deciding conseusun seq for each cell\n"+log_findCons_cell

        if not noremove:
            umi_splitted_fastapath=glob.glob(outdir+outname+".umi.*.fa")
            umi_splitted_alignmentpath=glob.glob(outdir+outname+".cdhit*umi.*")
            for f in umi_splitted_fastapath+umi_splitted_alignmentpath:
                os.remove(f)

    if passed_cell:
        #Align barcode sequence to reference file and add MD tag to output sam file
        outfile_cons=outdir+outname+".RESULT.CONSENSUS.fa"
        needle_result=needleAlign(outfile_cons,needleReference,outname,outdir)
        taggedSam=fixSam(needle_result,needleReference,outname,outdir)

        #Parse sam file
        parsedMutationPattern=parseSam(taggedSam,idenrate,identity,outname,outdir)
        outAlignedFile=outdir+outname+".ALIGNED.tsv"
        with open(outAlignedFile,mode="wt") as w:
            w.write("\t".join([cellbarcode,"CIGAR",parsedMutationPattern["cigar"]])+"\n")
            w.write("\t".join([cellbarcode,"QUERY",parsedMutationPattern["seq"]])+"\n")
        
        #Quality filtering
        ref_annot={}
        with open(annotation,mode="rt") as r:
            for nrow,line in enumerate(r):
                if nrow>0:
                    line=line.split(",")
                    ref_annot[line[0]]=list(map(int,line[1:]))
        log_alignmentQC,passed_alignmentQC=alignmentQC(parsedMutationPattern,ref_annot,primer5ed,idenrate,identity,outname,outdir)
        log+="\n#Alignment quality check\n"+log_alignmentQC

    if passed_alignmentQC:
        mutationTable=generateMutationTable(parsedMutationPattern,ref_annot)
        with open(outdir+outname+".MUTATIONTABLE.tsv",mode="wt") as w:
            w.write(cellbarcode+"\t"+mutationTable+"\n")
        
    with open(outdir+outname+".log",mode="wt") as w:
        w.write(log)
        
    print("program completed!")
