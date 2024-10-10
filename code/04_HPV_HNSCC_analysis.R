setwd("~/HNSCC_PDO_Project")
library(reshape2)
library(ggplot2)
library(ggsci)
library(dplyr)

my_theme=theme_classic()+theme(plot.background=element_blank(), panel.background=element_blank(),panel.grid.minor=element_blank(),
                               panel.grid.major=element_blank(),axis.title.y = element_text(size=20,color="black"),
                               axis.title.x = element_text(size=20,color="black"),
                               axis.text.x=element_text(size=18,color="black"),axis.text.y=element_text(size=18,color="black"),
                               legend.title = element_text(hjust=0.5, face="bold",size=15),
                               legend.text = element_text(size=15,color="black"))

# meta program
modules <- readRDS("output/Fig5/nmf_metaprograms_sigtop40_human.RDS")  
modules <- lapply(modules, names) 
RHP8=modules$RHP8

#1) TCGA data
TCGA_corr=read.csv("data/tcga.areg.correlation.csv",row.names = 1,header = T) # TCGA: AREG correlation with other genes
TCGA_corr_positive=TCGA_corr[TCGA_corr$cor>=0,]
TCGA_corr_positive$corr_order=rank(TCGA_corr_positive$cor)

TCGA_corr_positive$label=ifelse(TCGA_corr_positive$column%in%RHP8,TCGA_corr_positive$column,"")

p1 <- ggplot() +
  geom_point(
    data = TCGA_corr_positive[!TCGA_corr_positive$column %in% RHP8,], 
    aes(x = corr_order, y = -log10(p.adj)),
    color = "grey",
    shape = 10,
    size = 0.1,
    stroke = 1.5) +
  geom_point(
    data = TCGA_corr_positive[TCGA_corr_positive$column %in% RHP8,],  
    aes(x = corr_order, y = -log10(p.adj)),
    color = "red",
    shape = 10,
    size = 1,
    stroke = 1.5) +
  geom_vline(xintercept=6148,lty=4,col="#666666",lwd=0.8) +
  geom_hline(yintercept = -log10(0.1),lty=4,col="#666666",lwd=0.8) +
  my_theme +
  xlab("Correlation Coefficient (Rank order)") +
  ylab("-log10(q-value)")

p1=p1 + ggrepel::geom_text_repel(
  data = TCGA_corr_positive,
  mapping = aes(x = corr_order, y = -log10(p.adj), label = label),
  size = 2,
  box.padding = 0.6,
  point.padding = 1,
  min.segment.length = 0.5,
  segment.color = "black",
  show.legend = FALSE,
  max.overlaps = 1000)

ggsave("output/Fig8/Fig8B.pdf",plot = p1,width = 6.5, height = 5.5)
write.csv(TCGA_corr_positive,"output/Fig8/TCGA_corr_positive_order.csv")


#2) HPA data
HPA_corr=read.csv("data/HPA.areg.correlation.csv",row.names = 1,header = T) # CCLE/Human protein atlas: AREG correlation with other genes
HPA_corr_positive=HPA_corr[HPA_corr$cor>=0,]
HPA_corr_positive$corr_order=rank(HPA_corr_positive$cor)

HPA_corr_positive$label=ifelse(HPA_corr_positive$column%in%RHP8,HPA_corr_positive$column,"")

p2 <- ggplot() +
  geom_point(
    data = HPA_corr_positive[!HPA_corr_positive$column %in% RHP8,],  
    aes(x = corr_order, y = -log10(p.adj)),
    color = "grey",
    shape = 20,
    size = 0.5,
    stroke = 1.5) +
  geom_point(
    data = HPA_corr_positive[HPA_corr_positive$column %in% RHP8,],  
    aes(x = corr_order, y = -log10(p.adj)),
    color = "red",
    shape = 20,
    size = 2,
    stroke = 1.5) +
  geom_vline(xintercept=7955,lty=4,col="#666666",lwd=0.8) +
  geom_hline(yintercept = -log10(0.1),lty=4,col="#666666",lwd=0.8) +
  my_theme +
  xlab("Correlation Coefficient (Rank order)") +
  ylab("-log10(q-value)")


p2=p2 + ggrepel::geom_text_repel(
  data = HPA_corr_positive,
  mapping = aes(x = corr_order, y = -log10(p.adj), label = label),
  size = 2,
  box.padding = 0.6,
  point.padding = 1,
  min.segment.length = 0.5,
  segment.color = "black",
  show.legend = FALSE,
  max.overlaps = 1000)

ggsave("output/Fig8/Fig8C.pdf",plot = p2,width = 6.5, height = 5.5)
write.csv(HPA_corr_positive,"output/Fig8/HPA_corr_positive_order.csv")


