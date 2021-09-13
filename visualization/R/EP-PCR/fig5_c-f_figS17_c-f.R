source("~/work/yachie/fractal_gestalt/FRACTAL_mutation_call/script/R/errorPCR_toolkit.R")
source("~/work/fractal_naoki/FRACTAL_mutation_call/script/R/DFS_src.R")
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggtree)
library(tidyr)
library(igraph)
library(ape)
library(treeio)
library(hash)
library(aplot)
library(ggforce)
library(brunnermunzel)


##################################
# treepath: lineage tree (newick)
# use_tips: list of tips to use
##################################
ntips <- 15000
use_tips <- args[3]


tree_imported <- read.newick(treepath)
tree_phylo <- as.phylo(tree_imported)
tips_consistent <- read.table(use_tips,stringsAsFactors = F,header=F)$V1
use_tips <- intersect(tree_phylo$tip.label,tips_consistent)
tree_phylo <- keep.tip(tree_phylo,use_tips)
tbl_tree <- as_tibble(tree_phylo)

#Depth first search to get number of tips assigned to each node, stack
node_to_n_tips <- fnc_dfs_Ntips_stack(tree_phylo,rootnode(tree_phylo))

#Generate a parent to children hash
parent_to_child <- gen_parent_to_child(tbl_tree)

#Get upstream
tbl_tree$parent <- tbl_tree$parent %>% as.character
tbl_tree$node <- tbl_tree$node %>% as.character
up_node_list <- get_upstream_tree_dfs(tbl_tree,
                                      as.character(rootnode(tree_phylo)),
                                      parent_to_child,
                                      node_to_n_tips,
                                      ntips)

#processing
upstream_tree_tbl <- filter(tbl_tree,node %in% up_node_list)
upstream_tree_tbl$parent <- upstream_tree_tbl$parent %>% as.integer
upstream_tree_tbl$node <- upstream_tree_tbl$node %>% as.integer
upstream_tree_tbl_filled <- fill_label(upstream_tree_tbl)
tree_list <- reassign_node(upstream_tree_tbl_filled)
clade_vec <- get_original_node_to_clade(tree_list)
upstream_tip <- filter(upstream_tree_tbl_filled,!node %in% parent)$label
upstream_tip_node <- clade_vec[is.element(clade_vec,upstream_tip)] %>% names %>% as.integer()
tbl_tree$parent <- as.integer(tbl_tree$parent)
tbl_tree$node <- as.integer(tbl_tree$node)
offs_table <- offspring(tbl_tree,upstream_tip_node,self_include = T)
offs_table <- tibble(id=names(offs_table),offs_table) %>% unnest(cols = c(offs_table))
all_node_to_clade_vec <- clade_vec[offs_table$id]
names(all_node_to_clade_vec) <- offs_table$node

tbl_tree_expand <- separate(tbl_tree,label,paste0("col",1:5),sep="_")
tbl_tree_expand$Clade <- all_node_to_clade_vec[as.character(tbl_tree_expand$node)]

count_table <- aggregate(tbl_tree_expand$col4,
                         by=list(tbl_tree_expand$Clade),table) %>% as.matrix %>% as.data.frame()
count_table <- data.frame(Clade=count_table[,1],
                          apply(count_table[,2:ncol(count_table)],c(1,2),as.integer)/rowSums(apply(count_table[,2:ncol(count_table)],c(1,2),as.integer)),
                          total_count=rowSums(apply(count_table[,2:ncol(count_table)],c(1,2),as.integer)))
count_long <- pivot_longer(count_table[,-ncol(count_table)],
                           cols = colnames(count_table[,-ncol(count_table)])[2:ncol(count_table[,-ncol(count_table)])])
count_long$name <- sub("x.","",count_long$name)

#Entropy
fnc_entropy <- function(x){
  E <- 0
  for(i in x){if(i!=0){E <- E+i*log2(i)}}
  return(-1*E)
}
clade_entropy <- aggregate(count_long$value,list(count_long$Clade),fnc_entropy) %>% as_tibble


seed <- 1
tbl_tree <- as_tibble(tree_phylo)
rand.color <- readRDS("/path/to/well.COLORCODE.96.bcat.thick.rds") #if CATNNB
rand.color <- readRDS("/path/to/well.COLORCODE.96.bet002.thick.rds") #if CATNNB


#Figure 4c
g1_all <- ggplot(count_long,aes(x=Clade,y=value,fill=name))+
  geom_bar(stat="identity",position = "stack")+
  theme_classic()+
  NoLegend()+
  scale_fill_manual(values = rand.color)+
  scale_x_discrete(position = "top")+
  scale_y_continuous(breaks=c(0,1),labels = c("1","0"))+
  theme(axis.text = element_blank(),
        axis.title= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0,0,0,0),units="mm"),
        legend.position='none')

g3_all <- ggplot(count_table,aes(x=Clade,y=total_count))+
  theme_classic()+
  ylab("N tips")+
  xlab("Clades")+
  geom_bar(stat = "identity",position = "dodge",fill="steelblue",col="steelblue",size=0.1)+
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks=element_blank(),
        plot.margin = unit(c(0,0,0,0),units="mm"),
        axis.line = element_blank())+
  scale_y_log10(breaks=c(1,100,3000))

g4_all <- ggplot(clade_entropy,aes(x=Group.1,y=x))+
  theme_classic()+
  ylab("Entropy")+
  geom_bar(stat = "identity",position = "dodge",fill="red",col="red",size=0.1)+
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks=element_blank(),
        plot.margin = unit(c(0,0,0,0),units="mm"),
        axis.line = element_blank()
  )+
  scale_y_continuous(breaks=c(0,5),limits=c(0,fnc_entropy(rep(1/96,96))))

g2_all <- ggtree(as.phylo(tree_list$tree),size=.2,branch.length = "none") +
  geom_tiplab(size=0)+
  scale_x_reverse()+
  coord_flip()+
  theme(plot.margin = unit(c(0,0,0,0),units="mm"),
        axis.title.y = element_text(angle=-90))

g_merge_all <- g1_all %>% insert_top(g3_all,height = 0.2) %>% insert_top(g2_all,height = 0.4) %>% insert_bottom(g4_all,height = 0.2) 
ggsave("/poth/to/output_dir/upstream_tree.tip15000.jpeg",g_merge_all,width = 11.6,height = 9.2,units = "cm")




# Figure 4d
tree_data <- as.phylo(tree_list$tree); colnames(tree_data$edge) <- c("V1","V2")
node_to_n_tips_again <- fnc_dfs_Ntips_stack(tree_data,rootnode(tree_data))

#Generate a parent to children hash
parent_to_child_again <- gen_parent_to_child(tree_list$tree)

#Get upstream
tbl_tree_again <- tree_list$tree
tbl_tree_again$parent <- tbl_tree_again$parent %>% as.character
tbl_tree_again$node <- tbl_tree_again$node %>% as.character
up_node_list_again <- get_upstream_tree_dfs(tbl_tree_again,
                                            as.character(rootnode(tree_data)),
                                            parent_to_child_again,
                                            node_to_n_tips_again,
                                            50)

#Get offspring table
upstream_tree_tbl_again <- filter(tbl_tree_again,node %in% up_node_list_again)
upstream_tree_tbl_again$parent <- upstream_tree_tbl_again$parent %>% as.integer
upstream_tree_tbl_again$node   <- upstream_tree_tbl_again$node %>% as.integer
upstream_tree_tbl_filled_again <- fill_label(upstream_tree_tbl_again)
tree_list_again <- reassign_node(upstream_tree_tbl_filled_again)
clade_vec_again <- get_original_node_to_clade(tree_list_again)
upstream_tip_again <- filter(upstream_tree_tbl_filled_again,!node %in% parent)$label
upstream_tip_node_again <- clade_vec_again[is.element(clade_vec_again,upstream_tip_again)] %>% names %>% as.integer()
tbl_tree_again$parent <- as.integer(tbl_tree_again$parent)
tbl_tree_again$node <- as.integer(tbl_tree_again$node)
offs_table_again <- offspring(tbl_tree_again,upstream_tip_node_again,self_include = T)
available_nodes <- c()
for(i in names(offs_table_again)){if(nrow(offs_table_again[[i]])>1){available_nodes <- c(available_nodes,i)}}


d = fortify(as.phylo(tree_list$tree))

d = subset(d, isTip)
ordered_tips <- with(d, label[order(y, decreasing=T)])
clade_to_uptree_node <- filter(tree_list$tree,!node %in% parent)$node; names(clade_to_uptree_node) <- filter(tree_list$tree,!node %in% parent)$label
ordered_tips <- clade_to_uptree_node[ordered_tips]

#bcat: 12130
#bet002: 1949

use_sub_clade <- "12130"

selected_nodes <- na.omit(offs_table_again[[use_sub_clade]])$node
selected_clades <- na.omit(offs_table_again[[use_sub_clade]])$label
tip_category <- rep("No",length(ordered_tips))
names(tip_category) <- ordered_tips
tip_category[as.character(selected_nodes)] <- "Yes1"
tip_col_df <- data.frame(node=names(tip_category) %>% as.integer,Selected=tip_category)
tip_color <- c(adjustcolor("white",alpha.f = 0),adjustcolor("white",alpha.f = 0),"orangered")
names(tip_color) <- c("No","Yes1","Yes2")


g1 <- ggplot(count_long %>% filter(Clade %in% selected_clades),aes(x=Clade,y=value,fill=name))+
  geom_bar(stat="identity",position = "stack")+
  theme_classic()+
  NoLegend()+
  scale_fill_manual(values = rand.color)+
  scale_x_discrete(position = "top")+
  scale_y_continuous(breaks=c(0,0.5,1),labels = c("1","0.5","0"))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 0,hjust = 0.95,size=7),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size=0.2),
        plot.margin = unit(c(0,0,0,0),units="mm"),
        legend.position='none')+
  ylab("Well portion per clade")

extracted_clade_tbl <- as.phylo(tree_list$tree) %>% extract.clade(node=use_sub_clade %>% as.integer) %>% as_tibble
clade_to_extracted_node <- na.omit(extracted_clade_tbl)$node; names(clade_to_extracted_node) <- na.omit(extracted_clade_tbl)$label
extracted_node_category <- rep("Yes1",length(clade_to_extracted_node)); names(extracted_node_category) <- na.omit(extracted_clade_tbl)$node

tip_dol_df_extracted <- data.frame(node=as.integer(names(extracted_node_category)),Selected=extracted_node_category,stringsAsFactors = F)
g2 <- ggtree(as.phylo(extracted_clade_tbl),size=0.1,branch.length="none") %<+% tip_dol_df_extracted+
  geom_tiplab(size=0)+
  scale_x_reverse()+
  coord_flip()+
  theme(plot.margin = unit(c(0,0,0,0),units="mm"),
        axis.title = element_blank())

g3 <- ggplot(count_table %>% filter(Clade %in% selected_clades),aes(x=Clade,y=total_count))+
  theme_classic()+
  ylab("N tips")+
  xlab("Clades")+
  geom_bar(stat = "identity",position = "dodge",fill="steelblue",col="steelblue",size=0.1)+
  theme(axis.text.x = element_blank(),
        axis.ticks.y = element_line(size=0.2),
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(0,0,0,0),units="mm"),
        axis.text.y = element_text(angle=0,hjust = 0.95,size=7),
        axis.title = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y=element_blank())+
  scale_y_log10(breaks=c(1,1e+1,1e+2,1e+3,1e+4),limits=c(1,15000))

g4 <- ggplot(clade_entropy %>% filter(Group.1 %in% selected_clades),aes(x=Group.1,y=x))+
  theme_classic()+
  ylab("Entropy")+
  geom_bar(stat = "identity",position = "dodge",fill="red",col="red",size=0.1)+
  theme(axis.text.x = element_blank(),
        plot.margin = unit(c(0,0,0,0),units="mm"),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_line(size=0.2),
        axis.text.y = element_text(angle=0,hjust = 0.95,size=7),
        axis.title.y = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank())+
  scale_y_continuous(breaks=c(0,5),limits=c(0,fnc_entropy(rep(1/96,96))))
g_merge <- g1 %>% insert_top(g3,height = 0.2) %>% insert_top(g2,height = 0.4) %>% insert_bottom(g4,height = 0.2) 
ggsave("/path/to/output_dir/upstream_tree.subtree.jpeg",g_merge,width = 2.52,height = 9.2,units = "cm")


## Figure 4e
count_table.df <- data.frame(Clade=count_table$Clade,total_count=count_table$total_count)
entropy_count_merge <- merge(count_table.df,clade_entropy,by.x=1,by.y=1)
entropy_count_merge <- entropy_count_merge %>% filter(total_count>=1000)

tbl_tree_expand <- filter(tbl_tree_expand,! node %in% parent)
set.seed(seed)
tbl_tree_expand$col4 <- sample(tbl_tree_expand$col4,nrow(tbl_tree_expand))
tbl_tree_expand <- tbl_tree_expand[tbl_tree_expand$Clade %in% entropy_count_merge$Clade,]
scr.count_table <- aggregate(tbl_tree_expand$col4,
                         by=list(tbl_tree_expand$Clade),table) %>% as.matrix %>% as.data.frame()
scr.count_table <- data.frame(Clade=scr.count_table[,1],
                          apply(scr.count_table[,2:ncol(scr.count_table)],c(1,2),as.integer)/rowSums(apply(scr.count_table[,2:ncol(scr.count_table)],c(1,2),as.integer)),
                          total_count=rowSums(apply(scr.count_table[,2:ncol(scr.count_table)],c(1,2),as.integer)))
scr.count_long <- pivot_longer(scr.count_table[,-ncol(scr.count_table)],
                           cols = colnames(scr.count_table[,-ncol(scr.count_table)])[2:ncol(scr.count_table[,-ncol(scr.count_table)])])
scr.count_long$name <- sub("x.","",scr.count_long$name)
scr.clade_entropy <- aggregate(scr.count_long$value,list(scr.count_long$Clade),fnc_entropy) %>% as_tibble

entropy_count_merge.scr <- merge(entropy_count_merge,scr.clade_entropy,by.x=1,by.y=1)
colnames(entropy_count_merge.scr) <- c("Clade","tot","Real","Scramble")
entropy_count_merge.scr.long <- pivot_longer(entropy_count_merge.scr,cols = c("Real","Scramble"))

entropy_col <- c("Real"=adjustcolor("red",alpha.f = 0.5),"Scramble"=adjustcolor("blue",alpha.f = 0.5))
g_hist <- ggplot(entropy_count_merge.scr.long,aes(x=value,fill=name))+
  geom_histogram(bins = 50,size=0.2,position = "identity")+
  theme_classic()+
  scale_fill_manual(values = entropy_col)+
  theme(axis.title=element_blank(),
        axis.line=element_line(size=0.2),
        axis.text = element_text(size=7),
        axis.ticks = element_line(size=0.2),
        plot.background = element_rect(fill = "transparent",color = NA))+
  NoLegend()

ggsave("/path/to/output_dir/entropy_hist.jpeg",g_hist,width=2.7,height=2.1,units = "cm")

t <- brunnermunzel.test(entropy_count_merge.scr$Real,entropy_count_merge.scr$Scramble)
print(t)



## Figure 4f
# CATTNNB: Clade3115
# BET002: Clade654
clade_use <- "Clade3115"
col_ATGCN <- c(G="#F2F059", C="#74B2B7", A="#79E5B7", T="#FF776C", N="grey77")
clade_to_original_node <- as.integer(names(clade_vec)); names(clade_to_original_node) <- clade_vec
node_use <- clade_to_original_node[clade_use]
tree_use <- tree_phylo %>% extract.clade(node=node_use %>% as.integer)
tree_use_tbl <- tree_use %>% as_tibble

##################################
# delsub: Aligned sequences with deletions and substitutions
# ins: Binary matrix of insertion hits for each cell
##################################

del_pat <- read.table(delsub,sep="\t",header = F)
colnames(del_pat) <- c("Cell",as.character(seq(1,ncol(del_pat)-1)))
del_pat <- filter(del_pat,Cell %in% filter(tree_use_tbl,!node %in% parent)$label)
ins_pat <- read.table(ins,sep="\t",header = F)

tree_use_tbl.tips <- filter(tree_use_tbl,! node %in% parent)
tip_to_well <- str_split_fixed(tree_use_tbl.tips$label,"_",5)[,4]
well_tbl <- data.frame(tree_use_tbl.tips$label,tip_to_well)

colnames(ins_pat) <- c("Cell",as.character(seq(1,ncol(ins_pat)-1)))
colnames(well_tbl) <- c("Cell","Well")

del_pat_long <- pivot_longer(del_pat,
                             cols = -Cell,
                             names_to = "pos",
                             values_to = "pattern_raw")
del_pat_long$pos <- as.integer(del_pat_long$pos)
del_pat_long$Cell <- factor(del_pat_long$Cell,levels = well_tbl$Cell)
del_pat_long$Del_Sub_pattern <- factor(del_pat_long$pattern_raw,
                                       levels = c("A","T","G","C","N","D","="))

ins_pat_long <- pivot_longer(ins_pat,
                             cols = -Cell,
                             names_to = "pattern_index",
                             values_to = "pattern_raw")
ins_pat_long$Insertion_pattern <- factor(ins_pat_long$pattern_raw,
                                         levels = c(1,0))

col_ATGCN[c("D","=")] <- c("firebrick","grey90")
colnames(well_tbl)<- c("Cell","Well")
well_tbl$Cell <- factor(well_tbl$Cell,levels=well_tbl$Cell)
well_tbl <- filter(well_tbl,Cell %in% del_pat$Cell)

g_annot <- ggplot(well_tbl,aes(x=1,y=Cell,fill=Well))+
  geom_tile()+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        plot.margin = unit(c(0,0,0,0),units="mm"),
        axis.ticks = element_blank())+
  scale_fill_manual(values = rand.color)+
  NoLegend()

g <- ggplot(del_pat_long,aes(x=pos,y=Cell,fill=Del_Sub_pattern))+
  theme_classic()+
  geom_tile()+
  scale_fill_manual(values = col_ATGCN)+
  theme(axis.title = element_blank(),
        axis.text.y  = element_blank(),
        axis.text.x = element_text(size=7),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(size=0.2),
        axis.line = element_blank(),
        legend.position = "right",
        plot.margin = unit(c(1,1,1,1),units="mm"))+
  NoLegend()

g_ins <- ggplot(ins_pat_long,aes(x=pattern_index,y=Cell,fill=Insertion_pattern))+
  theme_classic()+
  geom_tile()+
  scale_fill_manual(values = c("blue","grey90"))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0.2,0.2,0.2,0.2),units="mm"))+
  NoLegend()

tree_use_tbl$is_root <- ifelse(tree_use_tbl$parent==tree_use_tbl$node,"root","other")
gtree <- ggtree(tree_use,size=.2,branch.length = "none") %<+% tree_use_tbl + 
  scale_color_manual(values = c(adjustcolor("white",alpha.f = 0),"orangered"))+
  theme(plot.margin = unit(c(0,0,0,0),units="mm"))
gmerge_pat <- g %>% insert_right(g_ins, width=.4) %>% insert_left(g_annot,width = 0.15) %>% insert_left(gtree,width = .5)

ggsave("/path/to/output_dir/selectedClade.delsubins.nolength.jpeg",gmerge_pat,width = 8,height = 14.25,units = "cm",dpi=600)

