library(tidyverse)
library(viridis)
library(ggpubr)
library(Seurat)

#############################
#########CATNNB##############
#############################

tbl <- read.table("/path/to//bcat_summary_singletophit.tsv",header=F)
well.col <- readRDS("/path/to/well.COLORCODE.96.bcat.thick.rds")

cannot_detected <- setdiff(unique(tbl$V3), unique(tbl$V1))
tbl <- filter(tbl,! V3 %in% cannot_detected)

tip_order <- unique(tbl$V3) %>% sort %>% as.character 
tbl$V1 <- factor(tbl$V1,levels = tip_order)
tbl$V3 <- factor(tbl$V3,levels = rev(tip_order))

g_annot <- ggplot(tbl,aes(x=1,y=V3,fill=V3))+
  geom_tile()+
  theme_classic()+
  NoLegend()+
  scale_fill_manual(values=well.col)+
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.margin = unit(c(1,0.5,1,1),
                           units="mm"))

g1 <- ggplot(tbl,aes(x=V1,y=V3,fill=V2))+
  geom_tile()+
  theme_classic()+
  gradient_fill(c("black","red","yellow","white"))+
  theme(axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=7),
        plot.margin = unit(c(1,0.5,1,0.5),units="mm"))

g_merge <- ggarrange(g_annot,g1,ncol=2,widths = c(1,10),align = "h")

ggsave("/path/to/output_dir/heatmmap.jpeg",
       g_merge,
       width = 7.5,
       height = 5,
       units = "cm",
       dpi=600)




#############################
#########BET002##############
#############################

tbl <- read.table("/path/to/bet002_summary_singletophit.tsv",header=F)
well.col <- readRDS("/path/to/well.COLORCODE.96.bet002.thick.rds")

tip_order <- unique(tbl$V3) %>% sort %>% as.character
tbl$V1 <- factor(tbl$V1,levels = tip_order)
tbl$V3 <- factor(tbl$V3,levels = rev(tip_order))

g_annot <- ggplot(tbl,aes(x=1,y=V1,fill=V1))+
  geom_tile()+
  theme_classic()+
  NoLegend()+
  scale_fill_manual(values=well.col)+
  theme(axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.margin = unit(c(1,0.5,1,1),
                           units="mm"))

g1 <- ggplot(tbl,aes(x=V1,y=V3,fill=V2))+
  geom_tile()+
  theme_classic()+
  gradient_fill(c("black","red","yellow","white"))+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        legend.title = element_blank(),
        plot.margin = unit(c(1,0.5,1,0.5),units="mm"))

g_merge <- ggarrange(g_annot,g1,ncol=2,widths = c(1,10),align = "h")
ggsave("/path/to/output_dir/heatmap.jpeg",
       g_merge,
       width = 7.5,
       height = 5,
       units = "cm",
       dpi=600)

