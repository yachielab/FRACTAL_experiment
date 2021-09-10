library(tidyverse)
library(viridis)
library(ggpubr)
library(Seurat)

tbl <- read.table("~/work/yachie/fractal_gestalt/local/PCR_10022020/re_analysis/fractal_tree/tip_to_parent_lev/bcat_summary_singletophit.tsv",header=F)
well.col <- readRDS("~/work/yachie/fractal_gestalt/local/PCR_10022020/re_analysis/fractal_tree/well.COLORCODE.96.bcat.thick.rds")

# well.col.add <- c("black","grey20","grey30","grey50","grey60")
# names(well.col.add) <- setdiff(unique(tbl$V3),names(well.col))
# well.col <- c(well.col,well.col.add)

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
  # scale_fill_gradient2(low = "black",mid = "yellow",high = "white",midpoint = 0.4)+
  # scale_fill_viridis(option="B")+
  # scale_fill_brewer(palette = "YlGn")+
  gradient_fill(c("black","red","yellow","white"))+
  theme(# axis.ticks = element_line(size=0.2),
    # axis.ticks.length = unit(0.5,units = "mm"),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    # legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size=7),
    plot.margin = unit(c(1,0.5,1,0.5),units="mm"))

g_merge <- ggarrange(g_annot,g1,ncol=2,widths = c(1,10),align = "h")
# ggsave("~/work/yachie/Konno_FRACTAL/Revise data/Figures/Fig5 - PCR CTNNB1/Fig5b/bcat_tip_to_parent.nolength.legend.pdf",
#        g1,
#        width = 7.5,
#        height = 5,
#        units = "cm")
ggsave("~/work/yachie/Konno_FRACTAL/Revise round 1//Figures - raw data//Fig5 - PCR CTNNB1/Fig5b/bcat_tip_to_parent.singletophit.pdf",
       g_merge,
       width = 7.5,
       height = 5,
       units = "cm")
ggsave("~/work/yachie/Konno_FRACTAL/Revise round 1//Figures - raw data//Fig5 - PCR CTNNB1/Fig5b/bcat_tip_to_parent.singletophit.jpeg",
       g_merge,
       width = 7.5,
       height = 5,
       units = "cm",
       dpi=600)




#############################
#############################
#########BET002##############
#############################

tbl <- read.table("~/work/yachie/fractal_gestalt/local/PCR_10022020/re_analysis/fractal_tree/tip_to_parent_lev/bet002_summary_singletophit.tsv",header=F)
well.col <- readRDS("~/work/yachie/fractal_gestalt/local/PCR_10022020/re_analysis/fractal_tree/well.COLORCODE.96.unfiltered.thick.rds")

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
  # scale_fill_gradient2(low = "black",mid = "yellow",high = "white",midpoint = 0.4)+
  # scale_fill_viridis(option="B")+
  # scale_fill_brewer(palette = "YlGn")+
  gradient_fill(c("black","red","yellow","white"))+
  theme(axis.text = element_blank(),
        # axis.ticks = element_line(size=0.2),
        # axis.ticks.length = unit(0.5,units = "mm"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        # legend.position = "top",
        legend.title = element_blank(),
        # legend.text = element_text(size=7),
        plot.margin = unit(c(1,0.5,1,0.5),units="mm"))

g_merge <- ggarrange(g_annot,g1,ncol=2,widths = c(1,10),align = "h")
# ggsave("~/work/yachie/Konno_FRACTAL/Revise round 1//Figures - raw data//FigS17 - PCR/BET002_filtered/correspond_to_fig5b/bet002_tip_to_parent.nolength.legend.pdf",
#        g1,
#        width = 7.5,
#        height = 5,
#        units = "cm")
ggsave("~/work/yachie/Konno_FRACTAL/Revise round 1//Figures - raw data//FigS17 - PCR/BET002_filtered/correspond_to_fig5b/bet002_tip_to_parent.singletophit.pdf",
       g_merge,
       width = 7.5,
       height = 5,
       units = "cm")
ggsave("~/work/yachie/Konno_FRACTAL/Revise round 1//Figures - raw data//FigS17 - PCR/BET002_filtered/correspond_to_fig5b/bet002_tip_to_parent.singletophit.jpeg",
       g_merge,
       width = 7.5,
       height = 5,
       units = "cm",
       dpi=600)

