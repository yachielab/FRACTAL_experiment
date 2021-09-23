from Bio import Phylo
import sys
import subprocess
import os
import extraction_old

def PartialNormRFdist(treefile1, treefile2, codefile): #treefile1: complete tree file (newick), trefile2: partial tree file (newick), codefile: .R code which calculate NRFdist
    try:
    #if(True):
        tree2=Phylo.read(treefile2,'newick')
        namelist=list(c.name for c in tree2.get_terminals())
        nameset=set(namelist)
        tree1ext=extraction_old.tree_extraction(treefile1,nameset)
        Phylo.write(tree1ext,treefile1+".ext",'newick')
        namelist=list(c.name for c in tree1ext.get_terminals())
        nameset=set(namelist)
        tree2ext=extraction_old.tree_extraction(treefile2,nameset) # cross extraction to remove directory-named terminal
        Phylo.write(tree2ext,treefile2+".ext",'newick')
        NRFdist=float((subprocess.Popen("Rscript "+codefile+" "+treefile1+".ext"+" "+treefile2+".ext| grep [1]|cut -d\" \" -f2", stdout=subprocess.PIPE,shell=True).communicate()[0]).decode('utf-8'))
        print(str(len(namelist))+","+str(NRFdist))
    except:
        print("Error"+","+"Error")
    #os.remove(treefile1+".ext")
    #os.remove(treefile2+".ext")

if __name__ == "__main__":

    sys.setrecursionlimit(1000000000)

    args=sys.argv
    PartialNormRFdist(args[1],args[2],args[3])