library(ggplot2)
library(ggtree)
library(tidyr)
library(igraph)
library(ggraph)
library(oak)
library(Seurat)
library(treeio)

fnc_recall_tesst <- function(tbl_tree,out_tree=data.frame(),n_node,rootnode){
  if(nrow(out_tree)>n_node){
    return(out_tree)
  }else{
    out_tree <- subset(tbl_tree,parent %in% rootnode)
    return(Recall(tbl_tree,out_tree=out_tree,n_node,rootnode=out_tree$node))
    }
}


fnc_get_availabel_node <- function(tbl_tree,ntips_total=10,ntips_each=3){
  node_use <- c()
  for(p in tbl_tree$parent %>% unique){
    leaves <- offspring(tbl_tree,p,tiponly = T)
    if(nrow(leaves)>=ntips_total&&nrow(child(tbl_tree,p))==2){
      c_tbl <- child(tbl_tree,p)
      ofs_1 <- offspring(tbl_tree,c_tbl$node[1],tiponly=T)
      ofs_2 <- offspring(tbl_tree,c_tbl$node[2],tiponly=T)
      if(nrow(ofs_1)>=ntips_each && nrow(ofs_2)>=ntips_each){node_use <- c(node_use,p)}
    }
  }
  return(node_use)
}

fnc_calc_BBI <- function(tbl_tree,ntips_total=10,ntips_each=3){
  node_use <- c()
  bbi <- c()
  for(p in (tbl_tree$parent %>% unique)){
    leaves <- offspring(tbl_tree,p,tiponly=T)
    if(nrow(leaves)>=ntips_total&&nrow(child(tbl_tree,p))==2){
      c_tbl <- child(tbl_tree,p)
      ofs_1 <- offspring(tbl_tree,c_tbl$node[1],tiponly=T)
      ofs_2 <- offspring(tbl_tree,c_tbl$node[2],tiponly=T)
      if(nrow(ofs_1)>=ntips_each && nrow(ofs_2)>=ntips_each){
        node_use <- c(node_use,p)
        bbi[as.character(p)] <- min(nrow(ofs_1),nrow(ofs_2))/max(nrow(ofs_1),nrow(ofs_2))
      }
    }
  }
  
  if(oak::is.empty(node_use)){
    return(NA)
  }else{
    return(bbi)
  }
}

setwd("~/work/yachie/YUSUKE_KIJIMA.LAB/Experiments/Projects/fractal_gestalt/")

for(i in 1:3){
  delsub   <- paste0("local/GESTALT_mutation_call_0624/zf",i,"_gapOpen_10_delsub_per_cell.tsv")
  ins      <- paste0("local/GESTALT_mutation_call_0624/zf",i,"_gapOpen_10_insertion_per_cell.tsv")
  treepath <- paste0("local/naoki_tree_multifur/zf",i,"_0702_newick.txt")
  
  outfile  <- paste0("~/work/yachie/Konno_FRACTAL/Revise data/data/Yusuke/GESTALT//v3_ZF",i,"_lineage.jpeg")
  outfile_leg  <- paste0("~/work/yachie/Konno_FRACTAL/Revise data/data/Yusuke/GESTALT/v3_ZF",i,"_lineage.legend.pdf")
  outfile_pdf  <- paste0("~/work/yachie/Konno_FRACTAL/Revise data/data/Yusuke/GESTALT/v3_ZF",i,"_lineage.pdf")
  outfile_hist  <- paste0("~/work/yachie/Konno_FRACTAL/Revise data/data/Yusuke/GESTALT/v3_ZF",i,"_lineage.BBIhist.pdf")
  
  outfile_hist_lowThresh  <- paste0("~/work/yachie/Konno_FRACTAL/Revise data/data/Yusuke/GESTALT/v3_ZF",i,"_lineage.lowThresh.BBIhist.pdf")
  outfile_lowThresh  <- paste0("~/work/yachie/Konno_FRACTAL/Revise data/data/Yusuke/GESTALT/v3_ZF",i,"_lineage.lowThresh.jpeg")
  outfile_leg_lowThresh  <- paste0("~/work/yachie/Konno_FRACTAL/Revise data/data/Yusuke/GESTALT/v3_ZF",i,"_lineage.lowThresh.legend.pdf")

  
  #start
  del_pat <- read.table(delsub,sep="\t",header = F)
  ins_pat <- read.table(ins,sep="\t",header = F)
  tree_nwk <- read.tree(treepath)
  
  colnames(del_pat) <- c("Cell",as.character(seq(1,ncol(del_pat)-1)))
  colnames(ins_pat) <- c("Cell",as.character(seq(1,ncol(ins_pat)-1)))
  del_pat_long <- pivot_longer(del_pat,
                               cols = -Cell,
                               names_to = "pos",
                               values_to = "pattern_raw")
  del_pat_long$pos <- as.integer(del_pat_long$pos)
  ins_pat_long <- pivot_longer(ins_pat,
                               cols = -Cell,
                               names_to = "pattern_index",
                               values_to = "pattern_raw")
  del_pat_long$Del_Sub_pattern <- factor(del_pat_long$pattern_raw,
                                         levels = c("A","T","G","C","N","D","="))
  ins_pat_long$Insertion_pattern <- factor(ins_pat_long$pattern_raw,
                                           levels = c(1,0))
  library(RColorBrewer)
  # my_col_pstel <- c(brewer.pal(4,"Set2"),"Red","grey90")
  col_ATGCN[c("D","=")] <- c("firebrick","grey90")
  # col_ATGCN <- col_ATGCN[-which(names(col_ATGCN) %in% "-")]
  
  library(aplot)
  tree_imported <- read.newick(treepath)
  tree_phylo <- as.phylo(tree_imported)
  tbl_tree <- as_tibble(tree_imported)
  tbl_tree$label <- NA
  
  # for(condition in c("normal","low")){
  for(condition in c("normal")){
    if(condition=="normal"){
      ntips_total=10
      ntips_each=3
      out_jpeg <-outfile
      out_pdf <- outfile_pdf
      out_legend <- outfile_leg
      out_hist <- outfile_hist
    }else{
      ntips_total=2
      ntips_each=1
      out_jpeg <-outfile_lowThresh
      out_legend <- outfile_leg_lowThresh
      out_hist <- outfile_hist_lowThresh
    }
    
    node_use <- fnc_get_availabel_node(tbl_tree,ntips_total = ntips_total,ntips_each = ntips_each)
    
    bbi <- fnc_calc_BBI(tbl_tree,ntips_total = ntips_total,ntips_each = ntips_each)
    bbi_df <- data.frame(node=tbl_tree$node,bbi=NA,available="no",stringsAsFactors = F)
    for(i in as.integer(names(bbi))){
      bbi_df$bbi[bbi_df$node==i] <- bbi[as.character(i)]
      bbi_df$available[bbi_df$node==i] <- "yes"
    }
    
    gtree <- ggtree(tree_phylo,size=.1,branch.length = "none") +
      theme(axis.line.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position="top",
            plot.margin = unit(c(0,0,0,0),units="mm"))

    if(all(!is.na(bbi))){
      gtree <- gtree %<+% bbi_df +
        # geom_tiplab(size=0,
        #             align = T,
        #             linesize = 0.15)+
        geom_nodepoint(aes(size=bbi),shape=16,col=adjustcolor("orange1",alpha.f=0.4),na.rm = T)+
        geom_nodepoint(aes(col = available),shape=16,size=0.5)+
        scale_color_manual(values=c(adjustcolor("white",alpha.f = 0),"orangered"))+
        scale_size(range = c(1,8),limits = c(0,1) ,name="BBI")
      
      g_hist <- na.omit(bbi_df) %>% ggplot(aes(x=bbi))+
        geom_histogram(bins = 10,breaks=seq(0,1,by=0.1),fill=adjustcolor("orange1",alpha.f=0.6))+
        # xlim(c(0,1))+ 
        theme_classic()+
        scale_x_continuous(breaks = seq(0,1,by=0.1))+
        scale_y_continuous(breaks = seq(0,10,by=1))
      # ggsave(out_hist,g_hist,width = 9,height = 5,units = "cm")
    }
    
    g <- ggplot(del_pat_long,aes(x=pos,y=Cell,fill=Del_Sub_pattern))+
      theme_classic()+
      geom_tile()+
      scale_fill_manual(values = col_ATGCN)+
      # theme_tree2()+
      theme() +
      theme(axis.title = element_blank(),
            axis.text.y  = element_blank(),
            axis.text.x = element_text(size=7),
            axis.ticks.y = element_blank(),
            axis.ticks.x = element_line(size=0.2),
            axis.line = element_blank(),
            legend.position = "right",
            plot.margin = unit(c(1,1,1,1),units="mm"))
      
    
    g2 <- ggplot(ins_pat_long,aes(x=pattern_index,y=Cell,fill=Insertion_pattern))+
      theme_classic()+
      geom_tile()+
      scale_fill_manual(values = c("blue","grey90"))+
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.line = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            plot.margin = unit(c(0.2,0.2,0.2,0.2),units="mm"))
    
    gtree_merge <- g %>% insert_left(gtree,width=.4) %>% insert_right(g2, width=.5)
    ggsave(out_legend,gtree_merge,width = 20,height = 15,units = "cm")
    
    gtree_merge <- (g+NoLegend()) %>% insert_left(gtree+NoLegend(),width=.4) %>% insert_right(g2+NoLegend(), width=.5)
    ggsave(out_pdf,gtree_merge,width = 8,height = 15,units = "cm")
    
    gtree_merge <- (g+NoLegend()) %>% insert_left(gtree+NoLegend(),width=.4) %>% insert_right(g2+NoLegend(), width=.5)
    ggsave(out_jpeg,gtree_merge,width = 8,height = 15,units = "cm",dpi=600)
    
    
    print(paste0("Cell Number: ",nrow(del_pat)))
    print(paste0("BBI: ",
                 round(mean(fnc_calc_BBI(tbl_tree,ntips_total = ntips_total,ntips_each = ntips_each)),digits = 2),
                 "±",
                 round(sd(fnc_calc_BBI(tbl_tree,ntips_total = ntips_total,ntips_each = ntips_each)),digits = 2)))
  }
}


#GESTALT
treepath <- paste0("local/naoki_tree_multifur/original_GESTALT_newick.txt")
outfile  <- paste0("~/work/yachie/fractal_gestalt/manuscript_alias_Yusuke/GESTALT/v3_ORIGINAL_GESTALT.tree.pdf")
outfile_hist  <- paste0("~/work/yachie/fractal_gestalt/manuscript_alias_Yusuke/GESTALT/v3_ORIGINAL_GESTALT_lineage.BBIhist.pdf")

outfile_lowThresh  <- paste0("~/work/yachie/fractal_gestalt/manuscript_alias_Yusuke/GESTALT/v3_ORIGINAL_GESTALT.tree.lowThresh.pdf")
outfile_hist_lowThresh  <- paste0("~/work/yachie/fractal_gestalt/manuscript_alias_Yusuke/GESTALT/v3_ORIGINAL_GESTALT_lineage.lowThresh.BBIhist.pdf")

tree_imported <- read.newick(treepath)
tree_phylo <- as.phylo(tree_imported)
tbl_tree <- as_tibble(tree_imported)
tbl_tree$label <- NA

for(condition in c("normal","low")){
  if(condition=="normal"){
    ntips_total=10
    ntips_each=3
    out_jpeg <-outfile
    out_hist <- outfile_hist
  }else{
    ntips_total=2
    ntips_each=1
    out_jpeg <-outfile_lowThresh
    out_hist <- outfile_hist_lowThresh
  }
  
  node_use <- fnc_get_availabel_node(tbl_tree,ntips_total = ntips_total,ntips_each = ntips_each)
  
  bbi <- fnc_calc_BBI(tbl_tree,ntips_total = ntips_total,ntips_each = ntips_each)
  bbi_df <- data.frame(node=tbl_tree$node,bbi=NA,available="no",stringsAsFactors = F)
  for(i in as.integer(names(bbi))){
    bbi_df$bbi[bbi_df$node==i] <- bbi[as.character(i)]
    bbi_df$available[bbi_df$node==i] <- "yes"
  }
  
  gtree <- ggtree(tree_phylo,size=.1,branch.length = "none") +
    theme(axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position="top")+
    NoLegend()
  
  if(all(!is.na(bbi))){
    gtree <- gtree %<+% bbi_df +
      # geom_tiplab(size=0,
      #             align = T,
      #             linesize = 0.15)+
      geom_nodepoint(aes(size=bbi),shape=16,col=adjustcolor("skyblue",alpha.f=0.6),na.rm = T)+
      # scale_color_manual(values=bubble_color)+
      scale_size(range = c(1,8),limits = c(0,1) ,name="BBI")
    
    g_hist <- na.omit(bbi_df) %>% ggplot(aes(x=bbi))+
      geom_histogram(bins = 10,breaks=seq(0,1,by=0.1),fill=adjustcolor("skyblue",alpha.f=0.6))+
      # xlim(c(0,1))+ 
      theme_classic()+
      scale_x_continuous(breaks = seq(0,1,by=0.1))+
      scale_y_continuous(breaks = seq(0,10,by=1))
    ggsave(out_hist,g_hist,width = 9,height = 5,units = "cm")
  }
  ggsave(out_jpeg,gtree,width = 9,height = 12,units = "cm")
  
  print(paste0("BBI: ",
               round(mean(fnc_calc_BBI(tbl_tree,ntips_total = ntips_total,ntips_each = ntips_each)),digits = 2),
               "±",
               round(sd(fnc_calc_BBI(tbl_tree,ntips_total = ntips_total,ntips_each = ntips_each)),digits = 2)))
}
