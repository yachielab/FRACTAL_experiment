# RunTime_sim.py
# written by Naoki Konno

import matplotlib.pyplot as plt
import seaborn as sns
import copy
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import sys

class node():
    def __init__(self,t): 
        self.At=t

class node_queue():
    def __init__(self, node_num): # idM & mseq means id & sequence of mother ESU, respectively. 
        self.queue=list(node(0) for i in range(node_num))
        self.L=node_num
    def node_pop(self):
        self.L-=1
        no = self.queue.pop(0)
        return no
    def node_push(self, node):
        i=0
        if(self.L==0): self.queue.append(node)
        else:
            while(node.At<self.queue[self.L-1-i].At): 
                i+=1
                if(i==self.L):break
            self.queue.insert(self.L-i,node)
        self.L+=1
        
class job():
    def __init__(self,t,tip_num):
        self.Bt=t
        self.N=tip_num

class job_queue():
    def __init__(self,tip_num):
        self.queue=[job(0,tip_num)]
        self.M=1 # number of jobs
    def job_pop(self):
        self.M-=1
        jo = self.queue.pop(0)
        return jo
    def job_push(self, job):
        i=0
        if(self.M==0): self.queue.append(job)
        else:
            while(job.Bt < self.queue[self.M-1-i].Bt): 
                i+=1
                if(i==self.M):break
            self.queue.insert(self.M-i,job)
        self.M+=1

class FRACTAL_parameters():
    def __init__(self,subsample_size,threshold):
        self.subsample_size=subsample_size
        self.threshold=threshold

def f(N,a,b): # phylogenetic tree reconstruction + phylogenetic placement (divided by number of nodes)
    return a+b*N

def SimulateRunTime(node_num,tip_num,parameters, model, a,b):
    NODES=node_queue(node_num)
    JOBS=job_queue(tip_num)
    while(True):
        if(model==1) :
            no=NODES.node_pop()
            jo=JOBS.job_pop()

            # calculate run time
            RUNTIME=f(jo.N,a,b)

            finish=max(no.At,jo.Bt)+RUNTIME

            no.At=finish; NODES.node_push(no)
        elif(model==2): 
            jo=JOBS.job_pop()
            assigned_node_num=max(1,int(node_num*jo.N/tip_num)) # number of node used for computing assigned node number 
            node_list=[]
            for _ in range(assigned_node_num): node_list.append(NODES.node_pop())

            RUNTIME=f(jo.N,a,b)/assigned_node_num
            finish=max(max([no.At for no in node_list]),jo.Bt)+RUNTIME
            
            for no in node_list: 
                no.At=finish; NODES.node_push(no)

        if(jo.N>parameters.threshold):
            clusters=2 # number of clusters in a cycle
            jo.Bt=finish; jo.N = int(jo.N/clusters)
            for _ in range(clusters): JOBS.job_push(copy.deepcopy(jo))
        #finalize
        if(JOBS.M<=0):
            return max([no.At for no in NODES.queue])

param=FRACTAL_parameters(int(sys.argv[3]),int(sys.argv[4]))
RunTime=SimulateRunTime(int(sys.argv[1]),int(sys.argv[2]),param, int(sys.argv[5]), float(sys.argv[6]),float(sys.argv[7]))
print(str(RunTime)+",sec")