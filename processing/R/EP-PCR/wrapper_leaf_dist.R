args <- commandArgs(trailingOnly = TRUE)
workdir <- args[1]
tree.nwk <- args[2]
tip_to_parentalseq <- args[3]
n_combi <- as.integer(args[4])
n_per_seq <- as.integer(args[5])
seed <- as.integer(args[6])


source("/path/to/R/errorPCR_toolkit.R")
setwd(workdir)

print("start...")
tree_imported <- read.newick(tree.nwk)
tree_phylo <- as.phylo(tree_imported)
tree_tbl<- as_tibble(tree_imported)

tip_to_par <- read.table(tip_to_parentalseq)
tip_to_parseq <- tip_to_par$V2
names(tip_to_parseq) <- tip_to_par$V1

tree_tbl_tips <- tree_tbl %>% filter(!node %in% parent)
tree_tbl_tips$parseq <- tip_to_parseq[tree_tbl_tips$label]

large.parent <- aggregate(rep(1,nrow(tree_tbl_tips)),list(tree_tbl_tips$parseq),sum) %>% filter(x>=n_per_seq)
large.parent <- large.parent$Group.1

set.seed(seed)
tree_tbl_tips.large <- tree_tbl_tips %>% filter(parseq %in% large.parent) %>% group_by(parseq) %>% sample_n(n_per_seq)
tree_tbl_tips.small <- tree_tbl_tips %>% filter(! parseq %in% large.parent)

if(nrow(tree_tbl_tips.small)>0){
    tree_tbl_tips <- rbind(tree_tbl_tips.large %>% as.data.frame,tree_tbl_tips.small %>% as.data.frame)
}else{
    tree_tbl_tips <- as.data.frame(tree_tbl_tips.large)
}

tree.sub <- keep.tip(tree_phylo,tree_tbl_tips$label)
tree_tbl.sub <- as_tibble(tree.sub)

tree_tbl.tips <- filter(tree_tbl,! node %in% parent)
tree_tbl.sub.tips <- filter(tree_tbl.sub,! node %in% parent)
tree_tbl.orig.sub.merge <- merge(tree_tbl.tips,tree_tbl.sub.tips,by.x=4,by.y=4)
node_orig_to_sub.tips <- tree_tbl.orig.sub.merge$node.y ; names(node_orig_to_sub.tips) <- tree_tbl.orig.sub.merge$node.x
node_orig_to_sub.parent <- tree_tbl.orig.sub.merge$parent.y ; names(node_orig_to_sub.parent) <- tree_tbl.orig.sub.merge$parent.x

tree_tbl_tips$parent <- node_orig_to_sub.parent[tree_tbl_tips$parent %>% as.character]
tree_tbl_tips$node   <- node_orig_to_sub.tips[tree_tbl_tips$node %>% as.character]
tree_tbl_tips <- as_tibble(tree_tbl_tips)

node_used <- tree_tbl_tips$node %>% as.integer
node_list <- lapply(node_used,eq_return)
names(node_list) <- paste0("tip",node_used)



tbl_list_by_parseq <- list()
for(k in unique(tree_tbl_tips$parseq)){
  tbl_list_by_parseq[[k]] <- filter(tree_tbl_tips,parseq==k)
}

use_parents <- tree_tbl_tips$parseq %>% unique %>% as.character()
cnt <- 0
combination_sampled_list <- list()

print(use_parents)

for(i in use_parents){
    l_tmp <- list()
    for(j in use_parents){
        cnt <- cnt+1
        #if(paste0(j,":",i) %in% names(combination_sampled_list)){
        #next
        #}
        grd <- expand.grid(tbl_list_by_parseq[[i]]$node,tbl_list_by_parseq[[j]]$node)
        if(i==j){
        grd <- filter(grd,Var1<Var2)
        }
        if(nrow(grd)>n_combi){
        set.seed(seed+cnt)
        grd <- grd[sample(1:nrow(grd),n_combi),]
        }
        grd <- paste0(grd$Var1,":",grd$Var2)
        l_tmp[[paste0(i,":",j)]]<- grd
        combination_sampled_list[[paste0(i,":",j)]] <- "DONE"
    }
    l_tmp <- stack(l_tmp) %>% as_tibble
    l_tmp$group <- i
    fname <- paste0(i,"_sampled_combination.tsv")
    write.table(l_tmp,fname,col.names=T,row.names=F,quote=F,sep="\t")

    #shell command
    cmd <- paste0("bash /path/to/shell/EP-PCR/exe_caldist.sh ",
                  fname,
                  " subsample_tree.rds ",
                  i)
    print(cmd)
    system(cmd)
    print(paste(cnt,"combinations were processed."))

}
print("DONE!")
