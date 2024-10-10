setwd("~/HNSCC_PDO_Project")

##### Module 1. Calculate the 8 RHP scores for each cell #####
# Method and function'seurat_functions_public.R' come from article: DOI: 10.1038/s41588-022-01141-9
source('code/seurat_functions_public.R') 

# my data
# meta program
modules <- readRDS("output/Fig5/nmf_metaprograms_sigtop40_human.RDS") 
modules <- lapply(modules, names) 

# load HNSCC PDO merged seurat objects
srat_harmony <- readRDS("data/srat_harmony_dims50_res0.5.RDS")

# colors.module
colors.module <- list("RHP1"="#f59999","RHP2"="skyblue","RHP3"="slateblue2","RHP4"="#fcbe6e","RHP5"="#b4d789","RHP6"="#1f78b4","RHP7"="#f47e20","RHP8"="#189d77")

srt <- srat_harmony 
modules_rand <- MakeRand(srt, db = modules, assay = "RNA",nrand = 3, nbin = 25) 

ini <- matrix(0,nrow = ncol(srt), ncol = length(modules))
rownames(ini) <- colnames(srt)
colnames(ini) <- names(modules)
srt@meta.data[,names(modules)]<- as.data.frame(ini)

for (m in names(modules)){
  print(m)
  tryCatch(expr = {
    srt = GeneToEnrichment(srt, db = modules[m], method = 'rand', db_rand = modules_rand[m]) 
  }, error = function(e){c()})
}
saveRDS(srt, file = 'output/Fig6/srt.scored_RHP.RDS')

scores <- srt@meta.data[,names(modules)]
scores$sample <- paste0("S",str_sub(row.names(scores),1,3))
scores$state_TCGA <- srt$state_TCGA
scores$cluster <- srt$seurat_clusters
saveRDS(scores,file = "output/Fig6/scores.RHP.RDS")

plot.list <- lapply(names(modules), function(m){
  h = FeaturePlot(srt, features = m, pt.size = 0.1, cols = c('grey',colors.module[m])) +
    NoAxes() + NoLegend() +
    theme(plot.title = element_text(size = 10))
})

pdf("output/Fig6/Fig6A.pdf",width = 10, height = 5)
print(CombinePlots(plot.list, ncol = 4)) 
dev.off()


##### Module 2. Calculation of RHP cell ratios between samples #####
library(reshape2)
library(ggplot2)
library(ggsci)
library(dplyr)

srt <- readRDS("output/Fig6/srt.scored_RHP.RDS")
scores <- readRDS("output/Fig6/scores.RHP.RDS")

my_theme <- theme_classic()+theme(plot.background=element_blank(), panel.background=element_blank(),panel.grid.minor=element_blank(),
                                  panel.grid.major=element_blank(),axis.title.y = element_text(size=20,color="black"),
                                  axis.title.x = element_text(size=20,color="black"),
                                  axis.text.x=element_text(size=18,color="black"),axis.text.y=element_text(size=18,color="black"),
                                  legend.title = element_text(hjust=0.5, face="bold",size=15),
                                  legend.text = element_text(size=15,color="black"))

cols=c("#f59999","skyblue","slateblue2","#fcbe6e","#b4d789","#1f78b4","#f47e20","#189d77")

data <- scores[,-c(10,11)]
data1 <- reshape2::melt(data,id.vars = 'sample',variable.name = "RHP",value.name = "scores")

# DotPlot
cellratio <- data1 %>%
  dplyr::group_by(sample,RHP) %>%
  dplyr::summarise('Percent expressed'=mean(scores>0.5,na.rm = TRUE),`Average expression`=mean(scores))
cellratio$sample=factor(cellratio$sample,levels = c("S37T","S42T","S46T","S49T","S51T","S56T","S56L","S58T"))

p <- ggplot(cellratio,aes(x=RHP,y=sample))+
  geom_point(aes(size=`Percent expressed`,
                 color=`Average expression`))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size=12,color="black",angle = 30,vjust = 0.85,hjust = 0.75),
        axis.text.y = element_text(size=12,color="black"),
        legend.title = element_text(hjust=0.5, size=14),
        legend.text = element_text(size=10,color="black"))+
  scale_color_gradient(low="lightgrey",high="blue")+
  labs(x=NULL,y=NULL)
ggsave("output/Fig6/Fig6B.pdf",plot = p,width = 7,height = 3.5)

# Barplot
p2 <- ggplot(cellratio,aes(x=sample,y=`Percent expressed`,fill=RHP))+
  geom_bar(stat = 'identity', 
           position = position_dodge(width=0.9), 
           width = 0.9,     
           color='black')+
  labs(x=NULL,y="cellratio")+ 
  theme_bw(base_size = 18)+  
  theme(axis.text = element_text(colour = 'black'))+
  scale_fill_manual(values = rep(cols,3))
ggsave("output/Fig6/FigS7B.pdf",plot = p2,width = 8,height = 3.5)



##### Module 3. RHP and drug resistance #####
library(dplyr)
library(ggplot2)
library(reshape2)
library(rstatix)
library(ggsci)

srt <- readRDS("output/Fig6/srt.scored_RHP.RDS")
scores <- readRDS("output/Fig6/scores.RHP.RDS")

# Resistant vs Sensitive
###### 1) Differences in RHPscores between resistant and sensitive groups ######
###### mean scores ######
data <- scores[,-c(10,11)]
data$group <- srt$Group

data2 <- reshape::melt(data,id.vars = c('sample','group'),variable.name = 'RHP',value.name = 'scores')
data2 <- dplyr::rename(data2,'RHP'='variable','scores'='value')
data2 <- data2%>%dplyr::group_by(sample,group,RHP)%>%
  dplyr::summarise(avg=mean(scores))

set.seed(123)
data2 %>% sample_n_by(group, RHP, size = 1)

# Paired T-test
stat.test <- data2 %>%
  group_by(RHP) %>%
  pairwise_t_test(
    avg ~ group, paired = F, 
    p.adjust.method = "bonferroni",
    alternative="less"
  ) 
stat.test

bxp <- ggboxplot(
  data2, x = "RHP", y = "avg",
  color = "group"
)

stat.test <- stat.test %>% add_xy_position(x = "RHP")
bxp + stat_pvalue_manual(
  stat.test, label = "p.adj.signif", 
  step.increase = 0.08
)
bxp <- bxp + stat_pvalue_manual(
  stat.test, label = "p.adj.signif", 
  step.increase = 0.08, hide.ns = TRUE  #, tip.length = 0
)
bxp <- bxp+xlab("")+ylab("RHP scores")+theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75))

pdf("output/Fig6/Fig6C.pdf",height = 3,width = 5)
bxp
dev.off()


###### cellratio ######
data <- scores[,-c(10,11)]
data$group <- srt$Group
data1 <- reshape2::melt(data,id.vars = c('sample','group'),variable.name = "RHP",value.name = "scores")

# cellratio
cellratio <- data1 %>%
  dplyr::group_by(group,RHP,sample) %>%
  dplyr::summarise('Percent expressed'=mean(scores>0.5,na.rm = TRUE),`Average expression`=mean(scores)) 

# Paired T-test
colnames(cellratio) <- c("group","RHP","sample","Percent_expressed","Average_expression")
stat.test <- cellratio %>%
  group_by(RHP) %>%
  pairwise_t_test(
    Percent_expressed ~ group, paired = F, 
    p.adjust.method = "bonferroni",
    alternative="less"
  ) 
stat.test

p <- ggboxplot(
  cellratio, x = "RHP", y = "Percent_expressed",
  color = "group"
)

stat.test <- stat.test %>% add_xy_position(x = "RHP")
p <- p + stat_pvalue_manual(
  stat.test, label = "p.adj.signif", 
  step.increase = 0.08, hide.ns = TRUE  
)
p <- p+xlab("")+ylab("RHP cellratio")+theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75))

pdf("output/Fig6/Fig6D.pdf",height = 3,width = 5)
p
dev.off()


###### 2) Correlation between RHP and IC50 in drug-resistant/sensitive groups ######
library(reshape2)
library(ggplot2)
library(ggsci)
library(dplyr)
library(ggpubr)

data <- scores[,-c(10,11)]
data$group <- srt$Group # resistant or sensitive
data1 <- reshape2::melt(data,id.vars = c('sample','group'),variable.name = "RHP",value.name = "scores")

# cellratio
cellratio <- data1 %>%
  dplyr::group_by(group,RHP,sample) %>%
  dplyr::summarise('Percent expressed'=mean(scores>0.5,na.rm = TRUE),`Average expression`=mean(scores))

###### mean scores ######
cellratio$IC50 <- ifelse(cellratio$sample%in%c("S37T"),8.992,ifelse(cellratio$sample%in%c("S42T"),19.92,ifelse(cellratio$sample%in%c("S46T"),7.111,ifelse(cellratio$sample%in%c("S51T"),1.166,ifelse(cellratio$sample%in%c("S56T"),20.92,ifelse(cellratio$sample%in%c("S58T"),15.2,ifelse(cellratio$sample%in%c("S49T"),28.89,2.778)))))))

my_theme <- theme_classic()+theme(plot.background=element_blank(), panel.background=element_blank(),panel.grid.minor=element_blank(),
                                  panel.grid.major=element_blank())

test1 <- cellratio[which(cellratio$RHP == "RHP1"),]
test2 <- cellratio[which(cellratio$RHP == "RHP2"),]
test3 <- cellratio[which(cellratio$RHP == "RHP3"),]
test4 <- cellratio[which(cellratio$RHP == "RHP4"),]
test5 <- cellratio[which(cellratio$RHP == "RHP5"),]
test6 <- cellratio[which(cellratio$RHP == "RHP6"),]
test7 <- cellratio[which(cellratio$RHP == "RHP7"),]
test8 <- cellratio[which(cellratio$RHP == "RHP8"),]

p1 <- ggplot()+ geom_point(data=test1,aes(x=`Average expression`,y=IC50,color=group),size=1,shape=15)+geom_smooth(data=test1,aes(x=`Average expression`,y=IC50),method=lm)+stat_cor(data=test1,aes(x=`Average expression`,y=IC50), method = "pearson",label.x = 0.1,label.y = 30,size=5)+xlab("RHP1 scores")+my_theme
p2 <- ggplot()+ geom_point(data=test2,aes(x=`Average expression`,y=IC50,color=group),size=1,shape=15)+geom_smooth(data=test2,aes(x=`Average expression`,y=IC50),method=lm)+stat_cor(data=test2,aes(x=`Average expression`,y=IC50), method = "pearson",label.x = 0.1,label.y = 30,size=5)+xlab("RHP2 scores")+my_theme
p3 <- ggplot()+ geom_point(data=test3,aes(x=`Average expression`,y=IC50,color=group),size=1,shape=15)+geom_smooth(data=test3,aes(x=`Average expression`,y=IC50),method=lm)+stat_cor(data=test3,aes(x=`Average expression`,y=IC50), method = "pearson",label.x = 0.1,label.y = 30,size=5)+xlab("RHP3 scores")+my_theme
p4 <- ggplot()+ geom_point(data=test4,aes(x=`Average expression`,y=IC50,color=group),size=1,shape=15)+geom_smooth(data=test4,aes(x=`Average expression`,y=IC50),method=lm)+stat_cor(data=test4,aes(x=`Average expression`,y=IC50), method = "pearson",label.x = 0.1,label.y = 30,size=5)+xlab("RHP4 scores")+my_theme
p5 <- ggplot()+ geom_point(data=test5,aes(x=`Average expression`,y=IC50,color=group),size=1,shape=15)+geom_smooth(data=test5,aes(x=`Average expression`,y=IC50),method=lm)+stat_cor(data=test5,aes(x=`Average expression`,y=IC50), method = "pearson",label.x = 0.1,label.y = 30,size=5)+xlab("RHP5 scores")+my_theme
p6 <- ggplot()+ geom_point(data=test6,aes(x=`Average expression`,y=IC50,color=group),size=1,shape=15)+geom_smooth(data=test6,aes(x=`Average expression`,y=IC50),method=lm)+stat_cor(data=test6,aes(x=`Average expression`,y=IC50), method = "pearson",label.x = 0.1,label.y = 30,size=5)+xlab("RHP6 scores")+my_theme
p7 <- ggplot()+ geom_point(data=test7,aes(x=`Average expression`,y=IC50,color=group),size=1,shape=15)+geom_smooth(data=test7,aes(x=`Average expression`,y=IC50),method=lm)+stat_cor(data=test7,aes(x=`Average expression`,y=IC50), method = "pearson",label.x = 0.1,label.y = 30,size=5)+xlab("RHP7 scores")+my_theme
p8 <- ggplot()+ geom_point(data=test8,aes(x=`Average expression`,y=IC50,color=group),size=1,shape=15)+geom_smooth(data=test8,aes(x=`Average expression`,y=IC50),method=lm)+stat_cor(data=test8,aes(x=`Average expression`,y=IC50), method = "pearson",label.x = 0.1,label.y = 30,size=5)+xlab("RHP8 scores")+my_theme

pdf("output/Fig6/Fig6E.pdf",width = 16,height = 7)
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, ncol = 4)
dev.off()


###### cellratio ######
test1 <- cellratio[which(cellratio$RHP == "RHP1"),]
test2 <- cellratio[which(cellratio$RHP == "RHP2"),]
test3 <- cellratio[which(cellratio$RHP == "RHP3"),]
test4 <- cellratio[which(cellratio$RHP == "RHP4"),]
test5 <- cellratio[which(cellratio$RHP == "RHP5"),]
test6 <- cellratio[which(cellratio$RHP == "RHP6"),]
test7 <- cellratio[which(cellratio$RHP == "RHP7"),]
test8 <- cellratio[which(cellratio$RHP == "RHP8"),]

p1 <- ggplot()+ geom_point(data=test1,aes(x=`Percent expressed`,y=IC50,color=group),size=1,shape=15)+geom_smooth(data=test1,aes(x=`Percent expressed`,y=IC50),method=lm)+stat_cor(data=test1,aes(x=`Percent expressed`,y=IC50), method = "pearson",label.x = 0.1,label.y = 30,size=5)+xlab("RHP1 cellratio")+my_theme
p2 <- ggplot()+ geom_point(data=test2,aes(x=`Percent expressed`,y=IC50,color=group),size=1,shape=15)+geom_smooth(data=test2,aes(x=`Percent expressed`,y=IC50),method=lm)+stat_cor(data=test2,aes(x=`Percent expressed`,y=IC50), method = "pearson",label.x = 0.1,label.y = 30,size=5)+xlab("RHP2 cellratio")+my_theme
p3 <- ggplot()+ geom_point(data=test3,aes(x=`Percent expressed`,y=IC50,color=group),size=1,shape=15)+geom_smooth(data=test3,aes(x=`Percent expressed`,y=IC50),method=lm)+stat_cor(data=test3,aes(x=`Percent expressed`,y=IC50), method = "pearson",label.x = 0.1,label.y = 30,size=5)+xlab("RHP3 cellratio")+my_theme
p4 <- ggplot()+ geom_point(data=test4,aes(x=`Percent expressed`,y=IC50,color=group),size=1,shape=15)+geom_smooth(data=test4,aes(x=`Percent expressed`,y=IC50),method=lm)+stat_cor(data=test4,aes(x=`Percent expressed`,y=IC50), method = "pearson",label.x = 0.1,label.y = 30,size=5)+xlab("RHP4 cellratio")+my_theme
p5 <- ggplot()+ geom_point(data=test5,aes(x=`Percent expressed`,y=IC50,color=group),size=1,shape=15)+geom_smooth(data=test5,aes(x=`Percent expressed`,y=IC50),method=lm)+stat_cor(data=test5,aes(x=`Percent expressed`,y=IC50), method = "pearson",label.x = 0.1,label.y = 30,size=5)+xlab("RHP5 cellratio")+my_theme
p6 <- ggplot()+ geom_point(data=test6,aes(x=`Percent expressed`,y=IC50,color=group),size=1,shape=15)+geom_smooth(data=test6,aes(x=`Percent expressed`,y=IC50),method=lm)+stat_cor(data=test6,aes(x=`Percent expressed`,y=IC50), method = "pearson",label.x = 0.1,label.y = 30,size=5)+xlab("RHP6 cellratio")+my_theme
p7 <- ggplot()+ geom_point(data=test7,aes(x=`Percent expressed`,y=IC50,color=group),size=1,shape=15)+geom_smooth(data=test7,aes(x=`Percent expressed`,y=IC50),method=lm)+stat_cor(data=test7,aes(x=`Percent expressed`,y=IC50), method = "pearson",label.x = 0.1,label.y = 30,size=5)+xlab("RHP7 cellratio")+my_theme
p8 <- ggplot()+ geom_point(data=test8,aes(x=`Percent expressed`,y=IC50,color=group),size=1,shape=15)+geom_smooth(data=test8,aes(x=`Percent expressed`,y=IC50),method=lm)+stat_cor(data=test8,aes(x=`Percent expressed`,y=IC50), method = "pearson",label.x = 0.1,label.y = 30,size=5)+xlab("RHP8 cellratio")+my_theme

pdf("output/Fig6/Fig6F.pdf",width = 16,height = 7)
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, ncol = 4)
dev.off()


###### 3) Comparison of RHP8 gene in drug-resistant/sensitive groups ######
library(ggplot2)
library(ggunchained)
modules <- readRDS("output/Fig5/nmf_metaprograms_sigtop40_human.RDS") 
modules <- lapply(modules, names) 

targetGeneExp <- srt@assays$RNA@data[modules$RHP8,]
targetGeneExp <- data.frame(targetGeneExp, check.names = F)
targetGeneExp <- data.frame(t(targetGeneExp), check.names = F)
targetGeneExp <- cbind(srt@meta.data[,c("seq_folder", "Group")], targetGeneExp)
targetGeneExp <- targetGeneExp[!targetGeneExp$Group%in%"others",]
targetGeneExp$Group <- factor(targetGeneExp$Group, levels=c("Sensitive", "Resistant"))
plotdata <- targetGeneExp[,c(2:ncol(targetGeneExp))]
plotdata <- reshape2::melt(plotdata)
colnames(plotdata) <- c("group", "gene", "expression")
plotdata$cells <- "cells"
p <- ggplot(plotdata, aes(x = cells, y = expression, fill = group, color=group)) +geom_split_violin()+facet_wrap(~gene,nrow=2)+theme_classic()+xlab("")
p <- p+ylab("log (TPM+1)")+scale_fill_manual(values=c("blue", "red"))+scale_color_manual(values=c("blue", "red"))
p <- p+theme(axis.text = element_text(size=12, color="black"), axis.title = element_text(size=14, color="black"))

pdf("output/FigS7/FigS7C.pdf", width=16, height=5)
print(p)
dev.off()



######  4) RHP8 genes expression vs IC50 ######
library(ggplot2)
library(dplyr)
library(ggpubr)

my_theme=theme_classic()+theme(plot.background=element_blank(), panel.background=element_blank(),panel.grid.minor=element_blank(),
                               panel.grid.major=element_blank(),axis.title.y = element_text(size=20,color="black"),
                               axis.title.x = element_text(size=20,color="black"),
                               axis.text.x=element_text(size=18,color="black"),axis.text.y=element_text(size=18,color="black"),
                               legend.title = element_text(hjust=0.5, face="bold",size=18),
                               legend.text = element_text(size=15,color="black"))

modules <- readRDS("output/Fig5/nmf_metaprograms_sigtop40_human.RDS") 
modules <- lapply(modules, names)

# Correlation between average expression of RHP8 genes and IC50
# read the matrix of single cell data after log-normalization
HNSCC_expr <- readRDS("data/HNSCC_expr.RDS")

common_genes  <- Reduce(intersect, lapply(HNSCC_expr, rownames))
HNSCC_expr <-  lapply(HNSCC_expr, function(x) x[common_genes,])
gene_expr <- lapply(HNSCC_expr, rowMeans)
data <- as.data.frame(gene_expr)

RHP8genes <- modules$RHP8
target_expr <- data[RHP8genes,]

Plot_list <- list()
for (x in 1:32) {
  gene=RHP8genes[x]
  expr=target_expr[gene,]
  expr=reshape2::melt(expr,variable.name = "sample",value.name = "expression")
  expr$IC50=ifelse(expr$sample%in%c("S37TO"),8.992,ifelse(expr$sample%in%c("S42TO"),19.92,ifelse(expr$sample%in%c("S46TO"),7.111,ifelse(expr$sample%in%c("S51TO"),1.166,ifelse(expr$sample%in%c("S56TO"),20.92,ifelse(expr$sample%in%c("S58TO"),15.2,ifelse(expr$sample%in%c("S49TO"),28.89,2.778)))))))
  cor(expr$expression,expr$IC50,method="pearson")
  p=ggplot(expr,aes(x=expression,y=IC50))+ geom_point(size=1,shape=15)+geom_smooth(method=lm)+stat_cor(data=expr, method = "pearson",label.x = min(expr$expression),label.y = 35,size=5)+xlab(RHP8genes[x])+my_theme
  Plot_list[[gene]]=p
}

pdf("output/Fig7/Fig7B.pdf",width = 32,height = 16)
p_mean <- do.call(grid.arrange, c(Plot_list, ncol = 8))
dev.off()

expr_ratio <- list()
for (i in names(HNSCC_expr)){
  a=rowMeans(HNSCC_expr[[i]]>0)
  expr_ratio[[i]]=a
}
expr_ratio_df=as.data.frame(expr_ratio)

RHP8genes <- modules$RHP8
target_expr_ratio <- expr_ratio_df[RHP8genes,]

Plot_list1 <- list()
for (x in 1:32) {
  gene=RHP8genes[x]
  expr=target_expr_ratio[gene,]
  expr=reshape2::melt(expr,variable.name = "sample",value.name = "expression")
  expr$IC50=ifelse(expr$sample%in%c("S37TO"),8.992,ifelse(expr$sample%in%c("S42TO"),19.92,ifelse(expr$sample%in%c("S46TO"),7.111,ifelse(expr$sample%in%c("S51TO"),1.166,ifelse(expr$sample%in%c("S56TO"),20.92,ifelse(expr$sample%in%c("S58TO"),15.2,ifelse(expr$sample%in%c("S49TO"),28.89,2.778)))))))
  cor(expr$expression,expr$IC50,method="pearson")
  p=ggplot(expr,aes(x=expression,y=IC50))+ geom_point(size=1,shape=15)+geom_smooth(method=lm)+stat_cor(data=expr, method = "pearson",label.x = min(expr$expression),label.y = 35,size=5)+xlab(RHP8genes[x])+my_theme
  Plot_list1[[gene]]=p
}
pdf("output/Fig7/Fig7B.pdf",width = 32,height = 16)
p_ratio <- do.call(grid.arrange, c(Plot_list1, ncol = 8))
dev.off()


###### 5) RHP8 scatter plot labeled AREG ######
my_theme=theme_classic()+theme(plot.background=element_blank(), panel.background=element_blank(),panel.grid.minor=element_blank(),
                               panel.grid.major=element_blank(),axis.title.y = element_text(size=20,color="black"),
                               axis.title.x = element_text(size=20,color="black"),
                               axis.text.x=element_text(size=18,color="black"),axis.text.y=element_text(size=18,color="black"),
                               legend.title = element_text(hjust=0.5, face="bold",size=18),
                               legend.text = element_text(size=15,color="black"))

#cellratio
cellratio=readxl::read_xlsx("output/Fig7/RHP8corr_log2FC.xlsx",sheet = "cellratio")  # 'avg_log2FC(Resistant_vs_Sensitive)' is calculated using the Seurat::FindMarkers function
cellratio=as.data.frame(cellratio)
rownames(cellratio) = cellratio$'...1'
cellratio = cellratio[,-1]

#cellratio$label=ifelse(cellratio$`avg_log2FC(Resistant_vs_Sensitive)` >= 1&cellratio$`R value`>0.7,rownames(cellratio),"")
cellratio$label=ifelse(rownames(cellratio)=="AREG",rownames(cellratio),"")
cellratio$label=ifelse(rownames(cellratio)%in%c("AREG","SFN","LAMA3","IL1B","KRT6A"),rownames(cellratio),"")

p1=ggplot(cellratio, aes(x=-log10(cellratio$`P value`), y=cellratio$`R value`)) + 
  geom_point(
    color="#223D6C",
    fill='#2873B3',
    shape=21,
    size=3.5,
    stroke =1.5)+
  my_theme+ 
  xlab("-Log10(Correlation P value (cellratio vs IC50))")+
  ylab("Correlation R value (cellratio vs IC50)")

p1=p1+ggrepel::geom_label_repel(
  mapping = aes(label=label),cellratio,
  size = 5, #
  box.padding = 0.6, 
  point.padding = 0.5, 
  min.segment.length = 0.3, 
  segment.color = "black", 
  show.legend = F)
ggsave("output/Fig7/Fig7A.pdf",plot = p1,width = 6.5, height = 5.5)


###### 6) AREG vs. HNSCC metastasis (GSE9349) ######
library(data.table)

exp <- read.table("data/GSE9349_series_matrix.txt", header = TRUE, row.names = 1, sep = "\t")
exp$probe_ID=rownames(exp)
probe=fread('data/GPL201-30390.txt',data.table = F)[,c(1,11)]
colnames(probe)=c("probe_ID","Gene")

data=merge(exp,probe,by="probe_ID",all = T)
data1=na.omit(data)
data1=data1[data1$Gene!="",]
data1=data1[,c(1,24,2:23)]

library(reshape)
data1=data1[,-1]
test1=reshape::melt(data1,id.vars = 'Gene',variable.name = "sample",value.name = "expression")
colnames(test1)=c('Gene',"sample","expression")

test1$group=ifelse(test1$sample%in%c("GSM237996", "GSM237997", "GSM237998", "GSM237999", "GSM238000", "GSM238001", "GSM238002", "GSM238003", "GSM238004", "GSM238005", "GSM238006"),"non-metastasized","metastasized")
test1$position=ifelse(test1$sample%in%c("GSM237996", "GSM237997", "GSM237998", "GSM238007", "GSM238008", "GSM238009", "GSM238010"), "mond", ifelse(test1$sample%in%c("GSM237999", "GSM238000", "GSM238001", "GSM238002", "GSM238011","GSM238012","GSM238013"),"orofarynx","larynx"))
table(test1$sample,test1$group,test1$position)

TargetGene_expr=test1[test1$Gene=="AREG",]
TargetGene_expr$group=factor(TargetGene_expr$group,levels=c("non-metastasized","metastasized"))
TargetGene_expr$position=factor(TargetGene_expr$position,levels = c("mond","orofarynx","larynx"))
#T-test
library(rstatix)
library(ggpubr)
stat.test <- TargetGene_expr %>%
  pairwise_t_test(
    #t_test(
    expression ~ group, paired = F, 
    p.adjust.method = "bonferroni",
    alternative="greater"
  ) 
stat.test

bxp <- ggboxplot(
  TargetGene_expr, x = "group", y = "expression",
  color = "group", palette = "jco"
)
stat.test <- stat.test %>% add_xy_position(x = "group")
bxp=bxp + stat_pvalue_manual(
  stat.test, label = "p.adj.signif", 
  step.increase = 0.08
)
bxp=bxp + stat_pvalue_manual(
  stat.test, label = "p.adj.signif", 
  step.increase = 0.08, hide.ns = TRUE  
)
bxp=bxp+xlab("")+ylab("log2(AREG normalized value)")+ggtitle("HNSCC tumors(GSE9349)") +
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none")+ 
  theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75))
bxp

#T-test
stat.test <- TargetGene_expr %>%
  group_by(position) %>%
  pairwise_t_test(
    #t_test(
    expression ~ group, paired = F, 
    p.adjust.method = "bonferroni",
    alternative="greater"
  ) 
stat.test

bxp <- ggboxplot(
  TargetGene_expr, x = "position", y = "expression",
  color = "group", palette = "jco"
)
stat.test <- stat.test %>% add_xy_position(x = "position")
bxp=bxp + stat_pvalue_manual(
  stat.test, label = "p.adj.signif", 
  step.increase = 0.08
)
bxp=bxp + stat_pvalue_manual(
  stat.test, label = "p.adj.signif", 
  step.increase = 0.08, hide.ns = TRUE  
)
bxp=bxp+xlab("")+ylab("log2(AREG RMA-normalized value)")+ggtitle("HNSCC tumors(GSE9349)") +
  theme(plot.title = element_text(hjust = 0.5))+ 
  theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75))

pdf("output/Fig7/Fig7G.pdf",height = 4.5,width = 6)
bxp
dev.off()

