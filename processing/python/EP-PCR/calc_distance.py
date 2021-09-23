import sys
import itertools as it
import Levenshtein as lv
import collections as cl
import pandas as pd

seq_tsv=sys.argv[1]
parent_tsv=sys.argv[2]
outpath=sys.argv[3]
ualn_tips=sys.argv[4]
#thresh_len=int(sys.argv[5])

def get_levdist(query,parent_db,label_now):
    d_list=[]
    dseq_list=[]
    res=[]
    res_seq=[]
    for p in parent_db:
        well=p[0].split("_")[0]
        seq=p[1]

        d=lv.distance(query,seq)
        d_list.append([well,d])
        dseq_list.append([p[0],d])
    d_list=sorted(d_list,key=lambda x:x[1],reverse=False)
    dseq_list=sorted(dseq_list,key=lambda x:x[1],reverse=False)

    for n,i in enumerate(d_list):
        if n==0:
            d=i[1]
        elif i[1]>d:
            break
        res.append(i[0])
        res_seq.append(dseq_list[n][0])


    cntr=cl.Counter(res)
    tot=sum(cntr.values())
    cntr={k:cntr[k]/tot for k in cntr}
    out_str=["\t".join([label_now]+list(map(str,list(i)))) for i in cntr.items()]

    cntr_seq=cl.Counter(res_seq)
    tot_seq=sum(cntr_seq.values())
    cntr_seq={k:cntr_seq[k]/tot_seq for k in cntr_seq}
    out_seq_str=["\t".join([label_now,i[0].split("_")[0]]+list(map(str,list(i)))) for i in cntr_seq.items()]
    return out_str,out_seq_str


label_to_seq=dict()
with open(seq_tsv,mode="rt") as r:
    for l in r:
        l=l.replace("\n","")
        l=l.split("\t")
        #if len(l[1])>=thresh_len:
        label_to_seq[l[0]]=l[1]

well_parent=[]
with open(parent_tsv,mode="rt") as r:
    for line in r:
        line=line.replace("\n","")
        line=line.split("\t")
        well_parent.append(line)

ualn_tip_set=set()
with open(ualn_tips,mode="rt") as r:
    for i in r:
        i=i.replace("\n","")
        ualn_tip_set.add(i)

print("starting...")
dist_out=[]
dist_seqlevel_out=[]
for l in label_to_seq:
    if not l in ualn_tip_set:
        continue
    dist_now,dist_seqlevel_now=get_levdist(label_to_seq[l],well_parent,l)
    dist_out.append("\n".join(dist_now))
    dist_seqlevel_out.append("\n".join(dist_seqlevel_now))

dist_out="\n".join(dist_out)+"\n"
dist_seqlevel_out="\n".join(dist_seqlevel_out)+"\n"

with open(outpath+"_dist.tsv",mode="wt") as w:
    w.write(dist_out)

with open(outpath+"_dist_seqlevel.tsv",mode="wt") as w:
    w.write(dist_seqlevel_out)

with open(outpath+"_thresh.tsv",mode="wt") as w:
    for i in label_to_seq:
        w.write(i+"\t"+label_to_seq[i]+"\n")

all_wells=[i[0].split("_")[0] for i in well_parent]
all_wells=set(all_wells)
df = pd.read_csv(outpath+"_dist_seqlevel.tsv",sep="\t",header=None)
df.columns=["tip","nearest","nearest_seq","prop"]
df=df[["nearest","prop"]]
df=df.query('prop==1')

df_g=df.groupby("nearest")
df_sum=df_g.sum()

well_now=list(label_to_seq.keys())[0].split("_")[3]
not_detected_wells=all_wells-set(df["nearest"])
not_detected={k:0 for k in not_detected_wells}
#not_detected=pd.Series(not_detected)
for w in not_detected_wells:
    df_sum.loc[w]=not_detected[w]
#df_sum=pd.concat([df_sum,not_detected])
df_sum=df_sum/sum(df["prop"])
df_sum["well"]=well_now
print(df_sum.head())
df_sum.to_csv(outpath+"_summary.tsv",sep="\t",header=False,index=True)
