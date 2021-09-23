# command line argument : "R --vanilla --quiet --args xxxx.nwk yyyy.nwk"
# what the program do : calculate Robinson Fould Distance between 2 trees

library("Rcpp")
library("ape")
library("castor")
library("phangorn")

nwk1 = commandArgs(trailingOnly=TRUE)[1] #get 1st argument
nwk2 = commandArgs(trailingOnly=TRUE)[2] #get 2nd argument

tree1 <- read.tree(nwk1)
#tree3 <- multifurcations_to_bifurcations(tree1)$tree
tree2 <- read.tree(nwk2)
#tree4 <- multifurcations_to_bifurcations(tree2)$tree

#tree1
#tree2
#tree3
#tree4
RF.dist(tree1,tree2, rooted=FALSE, normalize=TRUE)
#RF.dist(tree3,tree4, rooted=FALSE, normalize=TRUE)