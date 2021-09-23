library(Seurat)
contami <- read.table("contamination.summary",sep="\t",header = F,stringsAsFactors = F)
tmp <- subset(contami,V3=="Drop (unexpected parental sequence)")
contami$V3 <- factor(contami$V3,levels = c("PASS","Drop (unexpected parental sequence)","Drop (redundant best-hit parental sequence)"),
                          labels = c("PASS","Drop\n(unexpected parental\nsequence)","Drop\n(redundant best-hit\nparental sequences)"))

contami$V1 <- factor(contami$V1,levels = c(as.character(tmp[order(tmp$V2,decreasing = T),]$V1),"NULL"))
col_code <- readRDS("COLORCODE.rds") #CATNNB or BET002
col_code["NULL"] <- "black"
g <- ggplot(contami,aes(V3,V2,fill=V1))+
  theme_classic()+
  geom_bar(stat = "identity")+
  NoLegend()+
  scale_fill_manual(values = col_code)+
  theme(axis.title = element_blank(),
        axis.text = element_text(size=7),
        axis.line = element_line(size=0.2),
        axis.ticks = element_line(size=0.2))
plot(g)
