source("~/work/yachie/fractal_gestalt/FRACTAL_mutation_call/script/R/errorPCR_toolkit.R")

setwd("/Users/yusukekijima/work/yachie/fractal_gestalt/local/PCR_10022020/re_analysis/fractal_tree/")

# in_df_split <- readRDS("~/work/yachie/Konno_FRACTAL/Revise round 1//data/Yusuke/epPCR_reanalysis/RDS/bcat_filtered//bcat_NOLENGTH_in_df_vln_fullcombination.rds")
# in_df_split <- readRDS("~/work/yachie/Konno_FRACTAL/Revise round 1//data/Yusuke/epPCR_reanalysis/RDS/bet002_filtered/bet002_NOLENGTH_in_df_vln_fullcombination.rds")
# in_df_split <- as_tibble(in_df_split)

# in_df <- read.table("../distance_calc_withoutJackPot/nolength/bcat_ntips_nolength_dist.well.tsv.gz",stringsAsFactors = F)
in_df <- read.table("../distance_calc_withoutJackPot/nolength/bet002_ntips_nolength_dist.well.tsv.gz",stringsAsFactors = F)

parent_path <- "../../analysis_test_bigfile/BET002/bet002_parent_summary.tsv" #need new version of entropy?
# parent_path <- "parent_thresh/original/catb_parent_summary.tsv"
seed <- 1
n <- 1000

compensate <- function(in_df_split){
  possible_cmbn <- expand.grid(unique(in_df_split$V1),c("intra","inter")) %>% as_tibble()
  possible_cmbn <- paste0(possible_cmbn$Var1,"_",possible_cmbn$Var2)
  need_to_fill <- setdiff(possible_cmbn,in_df_split$typemerge)
  print(need_to_fill)
  
  for(i in need_to_fill){
    print(i)
    well=strsplit(i,"_")[[1]][1]
    type=strsplit(i,"_")[[1]][2]
    
    df_tmp <- filter(in_df_split,V2==well&type==type)
    df_tmp$V2 <- df_tmp$V1
    df_tmp$V1 <- well
    df_tmp$typemerge <- i
    in_df_split <- rbind(in_df_split,df_tmp)
  }
  return(in_df_split)
}

# in_df <- in_df[,c(4,5)] #for filtered
in_df <- in_df[,c(2,4)] ; colnames(in_df) <- c("V5","V4") #for unfiltered

in_df_split <- str_split_fixed(in_df$V5,pattern = ":",2)
in_df_split <- as_tibble(in_df_split)
in_df_split$val <- in_df$V4
in_df_split$type <- ifelse(in_df_split$V1==in_df_split$V2,"intra","inter")
in_df_split$typemerge <- paste0(in_df_split$V1,"_",in_df_split$type)

in_df_split <- compensate(in_df_split)

#Distances are measured with subsampled tree
parent_df <- read.table(parent_path)
# available_tips <- read.table("parent_thresh/bet002_parent_knee.60.available.header")
# available_tips <- read.table("parent_thresh/bcat_parent_knee.50.available.header")
parent_df$seq <- paste0(parent_df$V3,"_",sub("Pattern","",parent_df$V2))
# parent_df <- parent_df %>% filter(seq %in% available_tips$V1)

use_well_list <- in_df_split$V1 %>% as.character %>% unique
parent_df <- filter(parent_df,V3 %in% use_well_list)

l <- setdiff(unique(in_df_split$V1),parent_df$V3)
if(length(l)>0){
  in_df_split <- filter(in_df_split,(! V1 %in% l) & (! V2 %in% l))
}

#sampling
cnt <- in_df_split$typemerge %>% table
combi.small <- cnt[cnt<n] %>% names
in_df_split.small <- filter(in_df_split,typemerge %in% combi.small)
in_df_split <- filter(in_df_split,! typemerge %in% combi.small)

set.seed(seed)
in_df_split <- in_df_split %>% group_by(typemerge) %>% sample_n(n)
in_df_split <- rbind(as.data.frame(in_df_split),as.data.frame(in_df_split.small)) %>% as_tibble


#normalize
sum_aggr <- aggregate(parent_df$V1,list(parent_df$V3),sum)
sum_v <- sum_aggr$x; names(sum_v) <- sum_aggr$Group.1
parent_df$V1 <- parent_df$V1 / sum_v[as.character(parent_df$V3)]

#Entropy
parent_entropy <- aggregate(parent_df[,1],list(parent_df$V3),fnc_entropy)
well_parent_entropy_sorted <- as.character(parent_entropy[order(parent_entropy$x),]$Group.1)

#Re-label
count_pat <- 0
prev <- "a"
for(i in 1:nrow(parent_df)){
  if(parent_df[i,3]!=prev){
    count_pat <- 1
    prev <- parent_df[i,3]
  }else{
    count_pat <- count_pat + 1
  }
  
  parent_df[i,2] <- paste0("Pattern",count_pat)
}

#Ordering
in_df_split$V1 <- factor(in_df_split$V1 ,levels = well_parent_entropy_sorted)
parent_df$V3 <- factor(parent_df$V3,levels = well_parent_entropy_sorted)
parent_df$V2 <- factor(parent_df$V2,
                       levels = rev(paste0("Pattern",1:length(unique(parent_df$V2)))),
                       labels = rev(paste0("Top",1:length(unique(parent_df$V2)))))

#Plotting
# g1 <- ggplot(in_df_split,aes(x=V1,y=val))+
#   geom_violin(data=filter(in_df_split,type=="intra"),
#               aes(x=V1,y=val),
#               fill="steelblue4",
#               col="steelblue4",
#               alpha=0.7,
#               size=0,
#               scale="width",
#               position=position_dodge(.9))+
#   geom_violin(data=filter(in_df_split,type=="inter"),
#               aes(x=V1,y=val),
#               fill=adjustcolor("white",alpha.f = 0),
#               col="black",
#               size=0.1,
#               scale="width",
#               position=position_dodge(.9))+
#   theme_classic()+
#   NoLegend()+
#   scale_y_continuous(breaks = c(0,0.5,1),labels = c("0.0","0.5","1.0"))+
#   theme(axis.ticks.x = element_blank(),
#         axis.text.x = element_blank(),
#         legend.position = "top",
#         axis.title=element_blank(),
#         axis.line.x=element_blank(),
#         axis.line.y=element_line(size=0.2),
#         axis.ticks.y=element_line(size=0.2),
#         axis.text.y=element_text(size=7))+
#   stat_summary(data=filter(in_df_split,type=="intra"),
#                fun=median,
#                geom="crossbar",
#                size=0.3,
#                col="blue4",
#                position=position_dodge(.9))+
#   stat_summary(data=filter(in_df_split,type=="inter"),
#                fun=median,
#                geom="crossbar",
#                size=0.3,
#                col="gray50",
#                position=position_dodge(.9))

my_col <- brewer.pal(length(unique(parent_df$V2)),"Accent")
names(my_col) <- unique(parent_df$V2) %>% as.character
g2 <- ggplot(parent_df,aes(x=V3,y=V1,fill=V2))+
  geom_bar(stat="identity",position = "stack",width=0.8)+
  theme_classic()+
  NoLegend()+
  scale_fill_manual(values = my_col)+
  scale_y_continuous(breaks=c(0,0.5,1))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title=element_blank(),
        axis.line.y=element_line(size=0.2),
        axis.line.x=element_blank(),
        axis.ticks.y=element_line(size=0.2),
        axis.text.y=element_text(size=7))+
  ylab("Parental\nsequence proportion")

# saveRDS(in_df_split,"~/work/yachie/Konno_FRACTAL/Revise round 1//data/Yusuke/epPCR_reanalysis/RDS/bet002_filtered//bet002_NOLENGTH_in_df_vln_fullcombination.well.rds")
saveRDS(in_df_split,"~/work/yachie/Konno_FRACTAL/Revise round 1//data/Yusuke/epPCR_reanalysis/RDS/bet002_filtered//bet002_NOLENGTH_in_df_vln_fullcombination.well.rds")

in_df_split <- readRDS("~/work/yachie/Konno_FRACTAL/Revise round 1 - NBT-RA50592B/Data/Yusuke/epPCR_reanalysis/RDS/bcat_filtered/bcat_NOLENGTH_in_df_vln_fullcombination.rds")
# in_df_split <- readRDS("~/work/yachie/Konno_FRACTAL/Revise data/data/Yusuke/epPCR_reanalysis/RDS/bcat_filtered/bet002_NOLENGTH_in_df_vln_fullcombination.well.rds")


#Statistical test, non-parametric
library(brunnermunzel)
utest_vec <- c()
stats <- c()
updn <- c()
for(well in (in_df_split$V1 %>% unique)){
  dist_intra <- filter(in_df_split,V1==well&type=="intra")
  dist_inter <- filter(in_df_split,V1==well&type=="inter")
  v1 <- dist_intra$val
  v2 <- dist_inter$val
  u <- brunnermunzel.test(v1,v2)
  stats <- c(stats,unname(u$estimate))
  updn <- c(updn,ifelse(u$estimate>0.5,"dn","up"))
  u <- u$p.value
  utest_vec[well] <- u
}
utest_df <- data.frame(well=names(utest_vec),p.val=utest_vec)
utest_df$p.va.adj <- p.adjust(utest_df$p.val,method="bonferroni")
utest_df$well <- factor(utest_df$well,levels = well_parent_entropy_sorted)
utest_df$stats <- stats
utest_df$updn <- updn


utest_vec
# saveRDS(utest_df,"~/work/yachie/Konno_FRACTAL/Revise data/data/Yusuke/epPCR_reanalysis/RDS/bcat_unaligned_filtered_vln.rds")
# utest_df <- readRDS("~/work/yachie/Konno_FRACTAL/Revise data/data/Yusuke/RDS/fig5_violin/brunnermunzel_bcat.rds")
# g1_tested <- g1+
#   geom_point(data = utest_df %>% filter(p.va.adj<.05),
#              aes(x=well,y=1.1),
#              shape=16,
#              size=0.7,
#              col=adjustcolor("orange",alpha.f=0.8))
# 
# library(ggpubr)
# g_merge <- ggarrange(g2,g1_tested,nrow=2,align = "v",heights = c(1,2.5))
# ggsave("~/work/yachie/Konno_FRACTAL/Revise data/Figures/FigS17 - PCR/BET002_filtered/correspond_to_fig5h//bet002_unaligned_nolength_filtfullcombi.pdf",
#        g_merge,
#        width=25,
#        height=6.5
#        ,units = "cm")




#split violin

GeomSplitViolin <- ggproto(
  "GeomSplitViolin",
  GeomViolin,
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <- transform(data,
                      xminv = x - violinwidth * (x - xmin),
                      xmaxv = x + violinwidth * (xmax - x)
    )
    grp <- data[1, "group"]
    newdata <- plyr::arrange(
      transform(data, x = if (grp %% 2 == 1) xminv else xmaxv),
      if (grp %% 2 == 1) y else -y
    )
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
      quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
      aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
      aesthetics$alpha <- rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <- GeomPath$draw_panel(both, ...)
      ggplot2:::ggname(
        "geom_split_violin",
        grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob)
      )
    } else {
      ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
    }
  }
)

geom_split_violin <- function(mapping = NULL,
                              data = NULL,
                              stat = "ydensity",
                              position = "identity", ...,
                              draw_quantiles = NULL,
                              trim = TRUE,
                              scale = "area",
                              na.rm = FALSE,
                              show.legend = NA,
                              inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomSplitViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      draw_quantiles = draw_quantiles,
      na.rm = na.rm, ...
    )
  )
}

in_df_split$type <- factor(in_df_split$type,levels = c("intra","inter"))
in_df_split$logval <- log10(in_df_split$val)
colsets <- c("intra"=adjustcolor("blue2",alpha.f = 0.6),
             "inter"=adjustcolor("grey60",alpha.f = 0.6))
colsets_med <- c("intra"=adjustcolor("navy",alpha.f = 1),
             "inter"=adjustcolor("grey20",alpha.f = 1))
g1 <- ggplot(in_df_split,aes(V1,logval,fill=type))+
  theme_classic()+
  stat_summary(data=filter(in_df_split),
               aes(color=type),
               fun=median,
               geom="crossbar",
               size=0.2,
               # col="blue4",
               position="dodge")+
  # stat_summary(data=filter(in_df_split,type=="inter"),
  #              fun=median,
  #              geom="crossbar",
  #              size=0.2,
  #              # col="gray50",
  #              position=position_dodge(.9))+
  geom_split_violin(scale="width",size=0,col=adjustcolor("white",alpha.f = 0))+
  scale_fill_manual(values = colsets)+
  scale_color_manual(values = colsets_med)+
  NoLegend()+
  scale_y_continuous(breaks = c(-4,-3,-2,-1,0),labels = c("1e-4","1e-3","1e-2","1e-1","1"))+
  # scale_y_log10(breaks = c(1e-4,2e-4,1e-3,1e-2,1e-1,1),labels = c("1e-4","","1e-3","1e-2","1e-1","1"))+
  # scale_y_log10()+
  # annotation_logticks(T)+
  coord_cartesian(clip="off") +
  annotation_logticks(sides = "l",
                      short = unit(-0.05, "cm"),
                      mid = unit(-0.1, "cm"),
                      long = unit(-0.15, "cm"),
                      size=0.2)+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title=element_blank(),
        axis.line.x=element_blank(),
        axis.line.y=element_line(size=0.2),
        axis.ticks.y=element_line(size=0.2),
        axis.text.y=element_text(size=7))

# utest_df$type="0"
g1_tested <- g1+
  geom_point(data = utest_df %>% filter(p.va.adj<.05),
             aes(x=well,y=.3,fill=NA,shape=updn),
             size=0.7)+
  scale_shape_manual(values = c("dn"=6,"up"=2))

library(ggpubr)
g_merge <- ggarrange(g2,g1_tested,nrow=2,align = "v",heights = c(1,2.5))
ggsave("~/work/yachie/Konno_FRACTAL/Revise round 1//Figures/FigS17 - PCR/BET002_filtered/correspond_to_fig5h/bet002_unaligned_nolength_filtfullcombi.log.split.well.pdf",
       g_merge,
       width=25,
       height=6.5
       ,units = "cm")
# ggsave("~/work/yachie/Konno_FRACTAL/Revise round 1//Figures/Fig5 - PCR CTNNB1/Fig5h/bcat_unaligned_nolength_filtfullcombi.log.split.well.pdf",
#        g_merge,
#        width=25,
#        height=6.5
#        ,units = "cm")
