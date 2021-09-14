source("/path/to/R/EP-PCR/errorPCR_toolkit.R")


args <- commandArgs(trailingOnly = TRUE)
combi <- args[1]
tree_RDS <- args[2]
outname <- args[3]

combi_df <- read.table(combi,header=T,sep="\t",stringsAsFactors=F)
tree.sub <- readRDS(tree_RDS)

#get number of tips in each node
node_to_ntips <- fnc_dfs_Ntips_stack(tree.sub,rootnode(tree.sub))

#get path from tip to root
tip_to_node_hash <- path_to_root(tree.sub,rootnode(tree.sub))

combi_df <- filter(combi_df,values!=":")
combi_df$dist <- sapply(combi_df$values,calc_dist_ntips,tip_to_node_hash,node_to_ntips)/Ntip(tree.sub)

combi_median <- aggregate(combi_df$dist,list(combi_df$ind),median)
write.table(combi_median,paste0(outname,"_ntips_median.tsv"),sep="\t",col.names=F,row.names=F,quote=F)
write.table(combi_df,paste0(outname,"_ntips_distance.tsv"),sep="\t",col.names=F,row.names=F,quote=F)
