library(tidyverse)

setwd("~/work/yachie/fractal_gestalt/local/PCR_10022020/revise/")

bcat <- read.table("bcat.judge.tsv",sep="\t",header = FALSE)
bet002 <- read.table("bet002.judge.tsv",sep="\t",header = F)
g <- ggplot(bcat,aes(x=V2,fill=V3))+
  geom_histogram(bins = 40)+
  theme_classic()+
  scale_y_log10()

g <- ggplot(bet002,aes(x=V2,fill=V3))+
  geom_histogram(bins = 40)+
  theme_classic()+
  scale_y_log10()


well_design <- data.frame(label=paste0(rep(LETTERS[1:8],each=12),rep(1:12,8)),
                          x=rep(1:12,8),
                          y=rep(8:1,each=12))
well.unexpected.count <- read.table("bcat.drop.unexpected.wellcount",header=F,stringsAsFactors = F)
well_design <- merge(well_design,well.unexpected.count,by.x=1,by.y=1)
ggplot(well_design,aes(x,y,col=V2))+
  geom_point(size=5)+
  scale_color_viridis_c(option="C")+
  theme_light()+
  scale_x_continuous(breaks = 1:12)+
  scale_y_continuous(breaks = 1:8,labels = rev(LETTERS[1:8]))

count_dat <- read.table("bcat.drop.unexpected",header=F,stringsAsFactors = F)
count_dat <- count_dat[order(count_dat$V2,decreasing = T),]
count_dat$rank <- 1:nrow(count_dat)
ggplot(count_dat,aes(x=rank,y=V2))+
  geom_bar(stat = "identity")+
  theme_classic()

dist_mat <- read.table("bcat_parent_pair_dist.tsv",header=F,stringsAsFactors = F)
dist_mat$V1 <- factor(dist_mat$V1,levels = as.character(count_dat$V1))
dist_mat$V2 <- factor(dist_mat$V2,levels = as.character(count_dat$V1))

ggplot(dist_mat,aes(V1,V2,fill=V3))+
  geom_raster()+
  scale_fill_viridis_c()+
  theme_classic()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())



#contami bcat
library(Seurat)
bcat.contami <- read.table("bcat.contatmi.summary",sep="\t",header = F,stringsAsFactors = F)
tmp <- subset(bcat.contami,V3=="Drop (unexpected parental sequence)")
bcat.contami$V3 <- factor(bcat.contami$V3,levels = c("PASS","Drop (unexpected parental sequence)","Drop (redundant best-hit parental sequence)"),
                          labels = c("PASS","Drop\n(unexpected parental\nsequence)","Drop\n(redundant best-hit\nparental sequences)"))

bcat.contami$V1 <- factor(bcat.contami$V1,levels = c(as.character(tmp[order(tmp$V2,decreasing = T),]$V1),"NULL"))
col_code <- readRDS("../re_analysis/fractal_tree/well.COLORCODE.96.bcat.thick.rds")
col_code["NULL"] <- "black"
g <- ggplot(bcat.contami,aes(V3,V2,fill=V1))+
  theme_classic()+
  geom_bar(stat = "identity")+
  NoLegend()+
  scale_fill_manual(values = col_code)+
  theme(axis.title = element_blank(),
        axis.text = element_text(size=7),
        axis.line = element_line(size=0.2),
        axis.ticks = element_line(size=0.2))

ggsave("~/work/yachie/Konno_FRACTAL/Revise round 2 - NBT-RA50592C/Figures - raw data/FigS18 - PCR contamination/CATNNB_contamination.sorted.pdf",g,width = 8,height = 10,units = "cm")



bet002.contami <- read.table("bet002.contatmi.summary",sep="\t",header = F,stringsAsFactors = F)
tmp <- subset(bet002.contami,V3=="Drop (unexpected parental sequence)")
bet002.contami$V3 <- factor(bet002.contami$V3,levels = c("PASS","Drop (unexpected parental sequence)","Drop (redundant best-hit parental sequence)"),
                          labels = c("PASS","Drop\n(unexpected parental\nsequence)","Drop\n(redundant best-hit\nparental sequences)"))
bet002.contami$V1 <- factor(bet002.contami$V1,levels = c(as.character(tmp[order(tmp$V2,decreasing = T),]$V1),"NULL"))
col_code <- readRDS("../re_analysis/fractal_tree/well.COLORCODE.96.unfiltered.thick.rds")
col_code["NULL"] <- "black"
g <- ggplot(bet002.contami,aes(V3,V2,fill=V1))+
  theme_classic()+
  geom_bar(stat = "identity")+
  NoLegend()+
  scale_fill_manual(values = col_code)+
  theme(axis.title = element_blank(),
        axis.text = element_text(size=7),
        axis.line = element_line(size=0.2),
        axis.ticks = element_line(size=0.2))

ggsave("~/work/yachie/Konno_FRACTAL/Revise round 2 - NBT-RA50592C/Figures - raw data/FigS18 - PCR contamination/BET002_contamination.sorted.pdf",g,width = 8,height = 10,units = "cm")
