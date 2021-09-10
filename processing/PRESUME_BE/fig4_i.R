library(tidyverse)
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

nodes <- c("18186078","23485373","33059400")

for(node in nodes){
  tbl <- readRDS(paste0("~/work/yachie/fractal_gestalt/local/BEmodel/subtree_re//extracted_clade_BEmodel_",node,".rds"))
  # label_tip <- filter(tbl,! node %in% parent)$label
  # write.table(label_tip,paste0("~/work/yachie/fractal_gestalt/local/BEmodel/subtree/tip_labels_",node,"txt"),
  #             col.names = F,
  #             row.names=F,
  #             quote=F)
  tree_phylo <- as.phylo(tbl)
  
  g_list <- list()
  for(n in c(1,2,3,4,100)){
    delsub   <- paste0("~/work/yachie/fractal_gestalt/local/BEmodel/vis_data_re//node",
                       node,
                       "_chunk",
                       n,
                       "_delsub_per_cell.tsv.gz")
    
    del_pat <- read.table(delsub,sep="\t",header = F,stringsAsFactors = T,colClasses = c("character"))
    colnames(del_pat) <- c("Cell",as.character(seq(1,ncol(del_pat)-1)))
    del_pat <- filter(del_pat,Cell %in% filter(tbl,!node %in% parent)$label)
    
    del_pat_long <- pivot_longer(del_pat,
                                 cols = -Cell,
                                 names_to = "pos",
                                 values_to = "pattern_raw")
    del_pat_long$pos <- as.integer(del_pat_long$pos)
    # del_pat_long$Cell <- factor(del_pat_long$Cell,levels = well_tbl$Cell)
    del_pat_long$Del_Sub_pattern <- factor(del_pat_long$pattern_raw,
                                           levels = c("A","T","G","C","N","D","="))
    library(RColorBrewer)
    col_ATGCN[c("D","=")] <- c("firebrick","grey90")
    if(n!=100){
      g <- ggplot(del_pat_long,aes(x=Cell,y=pos,fill=Del_Sub_pattern))+
        theme_classic()+
        geom_tile()+
        scale_fill_manual(values = col_ATGCN)+
        # theme_tree2()+
        scale_y_reverse()+
        theme(axis.title = element_blank(),
              axis.text  = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              legend.position = "right",
              plot.margin = unit(c(0,0,0,0),units="mm"),
              legend.title = element_text(size=9),
              legend.text = element_text(size=9),
              panel.background = element_rect(fill = "transparent",color = NA),
              plot.background = element_rect(fill = "transparent",color = NA))+
        NoLegend()
      g_list[[as.character(n)]] <- g
    }else{
      g <- ggplot(del_pat_long,aes(x=Cell,y=pos,fill=Del_Sub_pattern))+
        theme_classic()+
        geom_tile()+
        scale_fill_manual(values = col_ATGCN)+
        # theme_tree2()+
        scale_y_reverse()+
        theme(axis.title = element_blank(),
              axis.text  = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              legend.position = "right",
              plot.margin = unit(c(5,0,0,0),units="mm"),
              legend.title = element_text(size=9),
              legend.text = element_text(size=9),
              panel.background = element_rect(fill = "transparent",color = NA),
              plot.background = element_rect(fill = "transparent",color = NA))+
        NoLegend()
      g_list[[as.character(n)]] <- g
    }
    
  }
  
  g_blank <- ggplot(data.frame(x=del_pat$Cell),aes(x=x,y="1"))+
    theme_void()+
    theme(plot.margin = unit(c(0,0,0,0),units="mm"))+
    geom_tile(fill="white")+
    coord_fixed()
                      
    
  gtree <- ggtree(tree_phylo,size=.2,branch.length = "none")+
    theme(plot.margin = unit(c(0,0,0,0),units="mm"),
          panel.background = element_rect(fill = "transparent",color = NA),
          plot.background = element_rect(fill = "transparent",color = NA))+
    scale_x_reverse()+
    coord_flip()
  
  gmerge_pat <- g_list[["1"]] %>% 
    insert_top(gtree,height=.5) %>%
    insert_bottom(g_list[["2"]]) %>% 
    insert_bottom(g_list[["3"]]) %>%
    insert_bottom(g_list[["4"]]) %>%
    # insert_bottom(g_blank) %>%
    insert_bottom(g_list[["100"]]) 
  
  ggsave(paste0("~/work/yachie/Konno_FRACTAL/Revise data/data/Yusuke/BEmodel/BEmodel_new_high_ratio_",node,".jpeg"),
         gmerge_pat,
         width = Ntip(tree_phylo)/10,
         height = 19.2,
         units = "cm",
         bg = "transparent",
         dpi=600)
  ggsave(paste0("~/work/yachie/Konno_FRACTAL/Revise data/data/Yusuke/BEmodel/BEmodel_new_high_ratio_",node,".pdf"),
         gmerge_pat,
         width = Ntip(tree_phylo)/10,
         height = 19.2,
         units = "cm",
         bg = "transparent",
         dpi=600)
  
}

#ROOT
array_nuc <- paste0(strrep("N",12),paste0(rep(paste0(strrep("N",24),"PPP",collaplse=""),9),collapse=""),collapse="")
g_root_list <- list()
col_new <- c("grey40","orangered")
names(col_new) <- c("N","P")
for(n in c(1,2,3,4,100)){
  root <- data.frame(root="root",nuc=array_nuc,stringsAsFactors = F)
  root_df <- substring(root$nuc,seq(1,nchar(root$nuc),1),seq(1,nchar(root$nuc),1)) %>% data.frame
  root_df$pos <- 1:nrow(root_df)
  colnames(root_df) <- c("nuc","pos")
  if(n==100){
    g <- ggplot(root_df,aes(x="X",y=pos,fill=nuc))+
      theme_classic()+
      geom_tile()+
      scale_fill_manual(values = col_new)+
      # theme_tree2()+
      scale_y_reverse()+
      theme(axis.title = element_blank(),
            axis.text  = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            legend.position = "right",
            plot.margin = unit(c(5,0,0,0),units="mm"),
            legend.title = element_text(size=9),
            legend.text = element_text(size=9))+
      NoLegend()
    g_root_list[[as.character(n)]] <- g
  }else{
    g <- ggplot(root_df,aes(x="X",y=pos,fill=nuc))+
      theme_classic()+
      geom_tile()+
      scale_fill_manual(values = col_new)+
      # theme_tree2()+
      scale_y_reverse()+
      theme(axis.title = element_blank(),
            axis.text  = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            legend.position = "right",
            plot.margin = unit(c(0,0,0,0),units="mm"),
            legend.title = element_text(size=9),
            legend.text = element_text(size=9))+
      NoLegend()
    g_root_list[[as.character(n)]] <- g
  }
  
}

gmerge_root <- g_root_list[["1"]] %>% 
  insert_bottom(g_root_list[["2"]]) %>% 
  insert_bottom(g_root_list[["3"]]) %>%
  insert_bottom(g_root_list[["4"]]) %>%
  insert_bottom(g_root_list[["100"]]) 

ggsave(paste0("~/work/yachie/Konno_FRACTAL/Revise data/data/Yusuke/BEmodel/BEmodel_high_ratio_root.jpeg"),
       gmerge_root,
       width = 0.73,
       height = 17.35,
       units = "cm",
       dpi=600)



