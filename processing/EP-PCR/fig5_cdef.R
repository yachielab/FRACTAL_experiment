source("~/work/yachie/fractal_gestalt/FRACTAL_mutation_call/script/R/errorPCR_toolkit.R")

seed <- 1
setwd("~/work/yachie/fractal_gestalt/local/PCR_10022020/re_analysis/fractal_tree/")
outpath <- "~/work/yachie/Konno_FRACTAL/Revise round 1//Figures/Fig5 - PCR CTNNB1/"
outrds <- "~/work/yachie/Konno_FRACTAL/Revise round 1/data/Yusuke/epPCR_reanalysis/RDS/bcat_filtered/"
treepath <- "bcat_unaligned/dfs_15000_nolength/bcat_unaligned_nolength_15000tipsNode.filteredTree.rds"
tree_imported <- readRDS(treepath)
# tree_imported <- read.newick(treepath)
tree_phylo <- as.phylo(tree_imported)
tbl_tree <- as_tibble(tree_imported)

count_table <-readRDS("bcat_unaligned/dfs_15000_nolength/bcat_unaligned_nolength_15000tipsNode.count_table.rds")
count_long <- readRDS("bcat_unaligned/dfs_15000_nolength/bcat_unaligned_nolength_15000tipsNode.count_long.rds")
tree_list <-  readRDS("bcat_unaligned/dfs_15000_nolength/bcat_unaligned_nolength_tree_list.15000tipsNode.rds")
clade_entropy <- readRDS("bcat_unaligned/dfs_15000_nolength/bcat_unaligned_nolength_15000tipsNode.clade_entropy.rds")
clade_vec <- readRDS("bcat_unaligned/dfs_15000_nolength/bcat_unaligned_nolength_clade_vec.15000tipsNode.rds")
count_exp <- readRDS("bcat_unaligned/dfs_15000_nolength/bcat_unaligned_nolength_tbl_tree_expand.15000tipsNode.rds")
# 
filter_thick_col <- function(x){
  val <- paste(as.vector(col2rgb(x))) %>% as.integer
  if(sum(val)<600){return(x)}else{return(NA)}
}
#Set color code
# color_all = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
# color_all <- sapply(color_all,filter_thick_col)
# color_all <- na.omit(color_all)
# 
# wells <- expand.grid(c("A","B","C","D","E","F","G","H"),1:12) %>% as.data.frame()
# wells <- paste0(wells$Var1,wells$Var2)
# set.seed(1)
# rand.color <- sample(color_all,length(wells))
# names(rand.color) <- wells
# rand.color["E10"] <- "blue3"
# saveRDS(rand.color,"well.COLORCODE.96.bcat.thick.rds")
rand.color <- readRDS("well.COLORCODE.96.bcat.thick.rds")
# rand.color <- readRDS("well.COLORCODE.96.unfiltered.thick.rds") #for BET002, filtered

#Panel C
g1_all <- ggplot(count_long,aes(x=Clade,y=value,fill=name))+
  geom_bar(stat="identity",position = "stack")+
  theme_classic()+
  NoLegend()+
  scale_fill_manual(values = rand.color)+
  scale_x_discrete(position = "top")+
  scale_y_continuous(breaks=c(0,1),labels = c("1","0"))+
  # ylab("Well portion per clade")+
  theme(axis.text = element_blank(),
        # axis.text.y = element_text(angle = -90,hjust = 0.5),
        axis.title= element_blank(),
        # axis.title.y = element_text(angle = -90),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0,0,0,0),units="mm"),
        legend.position='none')

g3_all <- ggplot(count_table,aes(x=Clade,y=total_count))+
  theme_classic()+
  ylab("N tips")+
  xlab("Clades")+
  geom_bar(stat = "identity",position = "dodge",fill="steelblue",col="steelblue",size=0.1)+
  # scale_y_continuous(trans = trans_reverser('log10'),breaks = c(1,3000))+
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks=element_blank(),
        # axis.ticks.y=element_line(size=0.5),
        plot.margin = unit(c(0,0,0,0),units="mm"),
        # axis.text.y = element_text(angle=-90,hjust = 0.5),
        # axis.title.y = element_text(angle=-90),
        axis.line = element_blank())+
  scale_y_log10(breaks=c(1,100,3000))

g4_all <- ggplot(clade_entropy,aes(x=Group.1,y=x))+
  theme_classic()+
  ylab("Entropy")+
  geom_bar(stat = "identity",position = "dodge",fill="red",col="red",size=0.1)+
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks=element_blank(),
        # axis.ticks.y=element_line(size=0.5),
        plot.margin = unit(c(0,0,0,0),units="mm"),
        # axis.text.y = element_text(angle=-90,hjust = 0.5),
        # axis.title.y = element_text(angle=-90),
        axis.line = element_blank()
        # axis.line = element_blank()
        # axis.text.x = element_blank(),
        # plot.margin = unit(c(0,0,0,0),units="mm"),
        # axis.ticks.x=element_blank(),
        # axis.text.y = element_text(angle=-90,hjust = 0.5),
        # axis.title.y = element_text(angle=-90),
        # axis.line.x = element_blank(),
        # axis.title.x = element_blank()
  )+
  scale_y_continuous(breaks=c(0,5),limits=c(0,fnc_entropy(rep(1/96,96))))
# scale_y_reverse(breaks=c(0,5))

# g2_all <- ggtree(as.phylo(tree_list$tree),size=.2,branch.length = "none") +
#   geom_tiplab(size=0)+
#   scale_x_reverse()+
#   coord_flip()+
#   theme(plot.margin = unit(c(0,0,0,0),units="mm"),
#         axis.title.y = element_text(angle=-90))
# saveRDS(g2_all,paste0(outpath,"/tree_all.rds"))
g2_all <- readRDS("bcat_unaligned///dfs_15000_nolength/tree_filtered_ggplot.rds")
# xlab("Upstream tree")+
# geom_tippoint(aes(col=Selected),size=1.4)+
# scale_color_manual(values=tip_color)
g_merge_all <- g1_all %>% insert_top(g3_all,height = 0.2) %>% insert_top(g2_all,height = 0.4) %>% insert_bottom(g4_all,height = 0.2) 
ggsave(paste0(outpath,"/Fig5c/","upstream_tree.tip15000.nolength.pdf"),g_merge_all,width = 11.6,height = 9.2,units = "cm")
ggsave(paste0(outpath,"/Fig5c/","upstream_tree.tip15000.nolength.jpeg"),g_merge_all,width = 11.6,height = 9.2,units = "cm")
# ggsave("~/work/yachie/Konno_FRACTAL/Revise data/data/Yusuke//epPCR/BET002/BET002_upstream_tree.tip3000.jpeg",g_merge_all,width = 11.6,height = 9.2,units = "cm",dpi = 600)

###############################
#Sub plot
#Depth first search to get number of tips assigned to each node, stack
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


# d = fortify(as.phylo(tree_list$tree))
# saveRDS(d,paste0(outrds,"/fortify.nolength.rds"))
d <- readRDS(paste0(outrds,"/fortify.nolength.rds"))

d = subset(d, isTip)
ordered_tips <- with(d, label[order(y, decreasing=T)])
clade_to_uptree_node <- filter(tree_list$tree,!node %in% parent)$node; names(clade_to_uptree_node) <- filter(tree_list$tree,!node %in% parent)$label
ordered_tips <- clade_to_uptree_node[ordered_tips]

#BET002_aligned: 16793
#bcat_aligned: 43407
#bcat_aligned/15000: 19380
#bet002_aligned/15000: 7087

#bcat_unaligned/filt: 18056
#bet002_unaligned/filt: 2610
#bcat_unaligned/unfiltered: 22230
#bet002_unaligned/unfiltered: 7327

#bcat_unaligned/filt/nolength: 12130
#bet002_unaligned/filt: 2610
plt.check <- function(use_sub_clade,
                      offs_table_again=offs_table_again,
                      ordered_tips=ordered_tips,
                      count_long=count_long){
  selected_nodes <- na.omit(offs_table_again[[use_sub_clade]])$node
  selected_clades <- na.omit(offs_table_again[[use_sub_clade]])$label
  tip_category <- rep("No",length(ordered_tips))
  names(tip_category) <- ordered_tips
  tip_category[as.character(selected_nodes)] <- "Yes1"
  tip_col_df <- data.frame(node=names(tip_category) %>% as.integer,Selected=tip_category)
  tip_color <- c(adjustcolor("white",alpha.f = 0),adjustcolor("white",alpha.f = 0),"orangered")
  names(tip_color) <- c("No","Yes1","Yes2")
  
  
  ggplot(count_long %>% filter(Clade %in% selected_clades),aes(x=Clade,y=value,fill=name))+
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
}
available_nodes
use_sub_clade <- "12130"
plt.check(use_sub_clade,offs_table_again=offs_table_again,ordered_tips=ordered_tips,count_long=count_long)

# [1] "11689" "12025" "12074" "12130" "12370" "12624" "12790" "13038" "13088" "13153" "13258" "13513" "13781" "13970"
# [15] "14214" "14263" "14300" "14605" "14728" "14965" "15241" "15294" "15563" "15741" "15901" "16142" "16279" "16438"
# [29] "16652" "16862" "16967" "17171" "17342" "17618" "17952" "18056" "18383" "18514" "18693" "18776" "18943" "19060"
# [43] "19211" "19474" "19527" "19588" "19857" "20032" "20081" "20197" "20247" "20318" "20372" "20465" "20647" "20773"
# [57] "20824" "20867" "21132" "21241" "21298" "21435" "21630" "21679" "21821" "21870" "22016" "22080" "22319" "22524"
# [71] "22573" "22612" "22658" "22747" "22824" "22927" "23131"

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
# extracted_node_category[clade_to_extracted_node[clade_use]] <- "Yes2"

tip_dol_df_extracted <- data.frame(node=as.integer(names(extracted_node_category)),Selected=extracted_node_category,stringsAsFactors = F)
g2 <- ggtree(as.phylo(extracted_clade_tbl),size=0.1,branch.length="none") %<+% tip_dol_df_extracted+
  geom_tiplab(size=0)+
  scale_x_reverse()+
  coord_flip()+
  theme(plot.margin = unit(c(0,0,0,0),units="mm"),
        axis.title = element_blank())
# xlab("Upstream tree")+
# geom_tippoint(aes(col=Selected),size=1.4)+
# scale_color_manual(values = tip_color)

library(ggforce)
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
  # annotation_logticks(outside = T,sides = "l")
# scale_y_continuous(trans = trans_reverser('log10'),breaks = c(1,3000))

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
# scale_y_reverse(breaks=c(0,5),limits=c(6,0))
library(aplot)
g_merge <- g1 %>% insert_top(g3,height = 0.2) %>% insert_top(g2,height = 0.4) %>% insert_bottom(g4,height = 0.2) 
ggsave(paste0(outpath,"/Fig5d/upstream_tree.subtree.nolength.pdf"),g_merge,width = 2.52,height = 9.2,units = "cm")
###########


## Mark the clade
use_tip_df <- tree_list$tree %>% filter(! node %in% parent)
use_tip_df$use <- "0"
subclade.20105 <- offspring(tree_list$tree,as.integer(use_sub_clade)) %>% na.omit
use_tip_df$use[use_tip_df$node %in% subclade.20105$node] <- "1"

g2_all_marked <- ggtree(as.phylo(tree_list$tree),size=.2,branch.length = "none") %<+% use_tip_df +
  geom_tiplab(size=0)+
  scale_x_reverse()+
  coord_flip()+
  theme(plot.margin = unit(c(0,0,0,0),units="mm"),
        axis.title.y = element_text(angle=-90))+
  geom_tippoint(aes(col=use),size=0.01)+
  scale_color_manual(values=c(adjustcolor("white",alpha.f = 0),"red"))
ggsave(paste0(outpath,"/Fig5c/upstream_tree.subclade.nolength.Marked.pdf"),g2_all_marked,width = 11.6,height = 3,units = "cm")
saveRDS(g2_all_marked,paste0(outrds,"/tree_all_marked.nolength.rds"))

(count_table %>% filter(Clade %in% selected_clades))[,c(1,ncol(count_table))]
#Show clade: Clade5973 for bet002_aligned
#Show clade: Clade7373 for bcat_aligned/15000
#Show clade: Clade12195 for bcat_unaligned/15000
#Show clade: Clade2923 for bet002_aligned/15000
#Show clade: Clade2032 for bet002_aligned/15000

#Show clade: Clade7724 for bcat_unaligned/filt
#Show clade: Clade1073 for bcat_unaligned/filt

#Show clade: Clade3115 for bcat_unaligned/filt_nolength

clade_use <- "Clade3115"

show_node <- (extracted_clade_tbl %>% filter(label==clade_use))$node
tip_dol_df_extracted[tip_dol_df_extracted$node==show_node,2] <- "Yes2"
g2 <- ggtree(as.phylo(extracted_clade_tbl),size=0.2,branch.length="none") %<+% tip_dol_df_extracted+
  geom_tiplab(size=0)+
  scale_x_reverse()+
  coord_flip()+
  theme(plot.margin = unit(c(0,0,0,0),units="mm"),
        axis.title = element_blank())+
  geom_tippoint(aes(col=Selected),size=1.4)+
  scale_color_manual(values = c(adjustcolor("white",alpha.f=0),"blue"))
g1 %>% insert_top(g3,height = 0.2) %>% insert_top(g2,height = 0.4) %>% insert_bottom(g4,height = 0.2)



## Entropy distribution, N>=1000
count_table.df <- data.frame(Clade=count_table$Clade,total_count=count_table$total_count)
entropy_count_merge <- merge(count_table.df,clade_entropy,by.x=1,by.y=1)
entropy_count_merge <- entropy_count_merge %>% filter(total_count>=1000)


count_exp <- filter(count_exp,! node %in% parent)
set.seed(seed)
count_exp$col4 <- sample(count_exp$col4,nrow(count_exp))
count_exp <- count_exp[count_exp$Clade %in% entropy_count_merge$Clade,]
scr.count_table <- aggregate(count_exp$col4,
                         by=list(count_exp$Clade),table) %>% as.matrix %>% as.data.frame()
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
# entropy_col2 <- c("Real"="red","Scramble"=adjustcolor("white",alpha.f = 0))
g_hist <- ggplot(entropy_count_merge.scr.long,aes(x=value,fill=name))+
  geom_histogram(bins = 50,size=0.2,position = "identity")+
  # scale_y_log10()+
  # geom_vline(xintercept = fnc_entropy(c(1,rep(0,95))),size=0.2,linetype="dashed")+
  # geom_vline(xintercept = fnc_entropy(rep(1/96,96)),size=0.2,linetype="dashed")+
  # scale_y_continuous(trans=scales::pseudo_log_trans(sigma = 100,base = 10),breaks=c(0,10,100))+
  # coord_cartesian(clip="off") +
  # annotation_logticks(sides = "l",
  #                     short = unit(-0.05, "cm"),
  #                     mid = unit(-0.1, "cm"),
  #                     long = unit(-0.15, "cm"),
  #                     size=0.2)+
  theme_classic()+
  # scale_color_manual(values = entropy_col2)+
  scale_fill_manual(values = entropy_col)+
  theme(axis.title=element_blank(),
        axis.line=element_line(size=0.2),
        axis.text = element_text(size=7),
        axis.ticks = element_line(size=0.2),
        plot.background = element_rect(fill = "transparent",color = NA))+
  NoLegend()

ggsave(paste0(outpath,"/Fig5e/entropy_hist.nolength.pdf"),g_hist,width=2.7,height=2.1,units = "cm")

library(brunnermunzel)
t <- brunnermunzel.test(entropy_count_merge.scr$Real,entropy_count_merge.scr$Scramble)
saveRDS(t,paste0(outpath,"/Fig5e/entropy_hist.brunnermunzel.nolength.rds"))





## DelSub plot
clade_to_original_node <- as.integer(names(clade_vec)); names(clade_to_original_node) <- clade_vec
node_use <- clade_to_original_node[clade_use]
tree_use <- tree_phylo %>% extract.clade(node=node_use %>% as.integer)
tree_use_tbl <- tree_use %>% as_tibble
write.table(filter(tree_use_tbl, ! node %in% parent)$label,
            paste0(outrds,"/",clade_use,"_labels.nolength.txt"),
            quote=F,
            row.names = F,
            col.names = F)
# tip_label_use <- tree_use_tbl$label[str_detect(tree_use_tbl$label,"Plate")]


delsub   <- "~/work/yachie/Konno_FRACTAL/Revise round 1//data/Yusuke/epPCR_reanalysis/data_for_delsum/bcat_Clade3115_nolength_delsub_per_cell.tsv.gz"
ins      <- "~/work/yachie/Konno_FRACTAL/Revise round 1//data/Yusuke/epPCR_reanalysis/data_for_delsum/bcat_Clade3115_nolength_insertion.tsv.gz"
# well_cor <- "local/PCR_10022020/analysis_test_bigfile/BET002/upstream_tree/bet002_header_well_correspondence.tsv"


del_pat <- read.table(delsub,sep="\t",header = F)
colnames(del_pat) <- c("Cell",as.character(seq(1,ncol(del_pat)-1)))
del_pat <- filter(del_pat,Cell %in% filter(tree_use_tbl,!node %in% parent)$label)
# write.table(del_pat,
#             paste0("local/PCR_10022020/analysis_test_bigfile/upstream_tree/bcat_rmdup_delsub_per_cell.",clade_use,".tsv"),
#             col.names=F,
#             row.names=F,
#             quote=F)

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

library(RColorBrewer)
source("~/.Rprofile")
col_ATGCN[c("D","=")] <- c("firebrick","grey90")
# col_ATGCN <- col_ATGCN[-which(names(col_ATGCN) %in% "-")]
colnames(well_tbl)<- c("Cell","Well")
well_tbl$Cell <- factor(well_tbl$Cell,levels=well_tbl$Cell)
well_tbl <- filter(well_tbl,Cell %in% del_pat$Cell)

saveRDS(del_pat_long,paste0(outrds,"/",clade_use,"_delsub_pattern.nolength.rds"))
saveRDS(well_tbl,paste0(outrds,"/",clade_use,"_wells.nolength.rds"))
saveRDS(ins_pat_long,paste0(outrds,"/",clade_use,"_ins_pattern.nolength.rds"))

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
  # theme_tree2()+
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
  # geom_nodepoint(aes(col=is_root),size=1.4)+
  scale_color_manual(values = c(adjustcolor("white",alpha.f = 0),"orangered"))+
  theme(plot.margin = unit(c(0,0,0,0),units="mm"))
# ggsave("manuscript_alias_Yusuke/epPCR/hoge.jpeg",gtree,width=10,height = 3,dpi=600)
gmerge_pat <- g %>% insert_right(g_ins, width=.4) %>% insert_left(g_annot,width = 0.15) %>% insert_left(gtree,width = .5)

ggsave(paste0(outpath,"/Fig5f/selectedClade.delsubins.nolength.jpeg"),gmerge_pat,width = 8,height = 14.25,units = "cm",dpi=600)
ggsave(paste0(outpath,"/Fig5f/selectedClade.delsubins.nolength.pdf"),gmerge_pat,width = 8,height = 14.25,units = "cm")

