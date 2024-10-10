setwd("~/HNSCC_PDO_Project")

my_theme <- theme_classic()+theme(plot.background=element_blank(), panel.background=element_blank(),panel.grid.minor=element_blank(),
                                  panel.grid.major=element_blank(),axis.title.y = element_text(size=20,color="black"),
                                  axis.title.x = element_text(size=20,color="black"),
                                  axis.text.x=element_text(size=18,color="black"),axis.text.y=element_text(size=18,color="black"),
                                  legend.title = element_text(hjust=0.5, face="bold",size=18),
                                  legend.text = element_text(size=15,color="black"))



##### Module 1. Marker genes expression #####-----------------------------------------------------------------------------------------------------------
library(Seurat)
library(SeuratObject)
# load merged seurat object
srat_harmony <- readRDS("data/srat_harmony_dims50_res0.5.RDS")

# (1) Sample distribution
cols <- c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#b4d789","#F08080","#1E90FF")
srat_harmony$seq_folder <- factor(srat_harmony$seq_folder,levels = c("S37TO","S42TO","S46TO","S49TO","S51TO","S56LNO","S56TO","S58TO"))
srat_harmony <- SetIdent(srat_harmony,value = "seq_folder")
p1.1 <- DimPlot(srat_harmony, reduction = "umap", cols=cols) + my_theme
ggsave(plot = p1.1, filename = "output/Fig4/Fig4C.pdf", width = 6.5, height = 4.5)

# (2) Marker genes expression
marker <- c("EPCAM","KRT14","KRT5","SPRR1B","MKI67","LAMA3")
p1.2 <- FeaturePlot(object = srat_harmony, features = marker, min.cutoff = "q9", order = T,ncol = 3) & theme(legend.position = "right")
ggsave(plot = p1.2, filename = "output/Fig4/Fig4D.pdf", width = 10.5, height = 6)

marker_celltype <- c('COL1A2', 'COL3A1', 'VWF', 'CDH5', 'CD19', 'CD79A', 'CD3D', 'CD3E', 'LYZ', 'CLEC4C')

p1.3 <- FeaturePlot(object = srat_harmony, features = marker_celltype, min.cutoff = "q9", order = T,ncol = 4)
ggsave(plot = p1.3, filename = "output/Fig7/FigS7A.pdf", width = 12, height = 7)



##### Module 2. TCGA subtype analysis - bulk RNA-seq #####---------------------------------------------------------------------------------------------------------------------
###### 1) Define molecular subtyping of 13 HNSCC PDOs based on bulk RNA-seq ######
library(GSVA)
library(ComplexHeatmap)
library(FactoMineR)
library(factoextra)
# load data
# bulk RNA-seq expression matrix
exprSet<-read.csv(file = 'data/bulk_expr.csv',header = T,row.names = 1) 
exprSet<-exprSet[,(1:13)]
exprSet <- exprSet[rowMeans(exprSet)>1,]
colnames(exprSet) <- c("S04TO","S06TO","S36TO","S37TO","S41TO","S42TO","S46TO","S49TO","S51TO","S55TO","S56TO","S56LNO","S58TO")

# TCGA subtype gene list
TCGA_HPV <- readRDS("data/TCGA_HPV.RDS")
TCGA_HPV <- TCGA_HPV[-1]

subtype_cols <- c("Basal"="#ff9389","Atypical"="#00d55c","Classical"="#83b7ff")
anot_cols <- list(subtype=subtype_cols)

# Scoring the expression matrix with TCGA subtype gene list using gsva
expr_marix <- as.matrix(exprSet) 
expr_marix <- expr_marix[,c("S04TO","S06TO","S58TO","S36TO","S42TO","S51TO","S55TO","S56LNO","S37TO","S41TO","S46TO","S49TO","S56TO")]
TCGA_gsva <- gsva(expr_marix, TCGA_HPV, kcdf="Gaussian",method = "gsva",parallel.sz=10) 
# meta
meta <- list()
meta$subtype <- c(rep("Classical",3),rep("Atypical",5),rep("Basal",5))
meta <- as.data.frame(meta) 
rownames(meta) <- colnames(expr_marix)
meta$subtype <- factor(meta$subtype,levels = c("Basal","Atypical","Classical"))

p2.1 <- pheatmap(TCGA_gsva,annotation_col = meta, annotation_colors = anot_cols ,cluster_rows = F,cluster_cols = F,name = "scores")
pdf("output/Fig4/Fig4A.pdf",width = 10,height = 3)
p2.1
dev.off()


###### 2) Principal component analysis of subtypes ######
# metadata: Information about whether the sample is resistant to drugs
samples <- read.table("data/HNSCC_samples.all.txt",header = T)

cols=c("#65b9fe","#00d759","#ff938d")

pd <- samples[,c("sample","condition")]
pd$condition <- factor(pd$condition,levels = c("Resistant","Sensitive"))
pd$subtype <- ifelse(pd$sample%in%c("S37TO","S41TO","S46TO","S49TO","S56TO"),"Basal",ifelse(pd$sample%in%c("S36TO","S42TO","S51TO","S55TO","S56LNO"),"Atypical","Classical"))
group_list <- factor(pd$subtype,levels = c("Basal","Atypical","Classical"))
names(group_list) <- pd$sample
group_list

dat <- as.data.frame(t(exprSet))
dat.pca <- FactoMineR::PCA(dat, graph = FALSE)
p3 <- factoextra::fviz_pca_ind(dat.pca, geom.ind = "point",
                           col.ind = group_list, 
                           addEllipses = TRUE, ellipse.type = "convex",
                           legend.title = "Groups",
                           pointsize = 3,
                           add.ind.names = FALSE,
                           palette = cols
)
ggsave(plot = p3, filename = "output/Fig4/Fig4B.pdf",width = 5.5,height = 4)



##### Module 3. TCGA subtype analysis - single cell RNA-seq #####--------------------------------------------------------------------------------------
###### 1) Subtype annotation of single cells ######
library(reshape2)
library(ggplot2)
library(ggsci)
library(dplyr)
library(msigdbr)
library(Seurat)
library(GSVA)
library(pheatmap)
library(patchwork)
library(openxlsx)
source('code/seurat_functions_public.R')

# TCGA subtype
TCGA_HPV <- readRDS("data/TCGA_HPV.RDS")

# load HNSCC PDO merged seurat objects
srat_harmony <- readRDS("data/srat_harmony_dims50_res0.5.RDS")
# colors.module
colors.module <- list("Basal"="#65b9fe","Atypical"="#00d759","Classical"="#ff938d")

srt <- srat_harmony
modules <- TCGA_HPV[-1]

modules_rand <- MakeRand(srt, db = modules, assay = "RNA",nrand = 3, nbin = 25) 

ini = matrix(0,nrow = ncol(srt), ncol = length(modules))
rownames(ini) = colnames(srt)
colnames(ini) = names(modules)
srt@meta.data[,names(modules)] = as.data.frame(ini)

for (m in names(modules)){
  print(m)
  tryCatch(expr = {
    srt = GeneToEnrichment(srt, db = modules[m], method = 'rand', db_rand = modules_rand[m]) 
  }, error = function(e){c()})
}
scores = srt@meta.data[,names(modules)]
frequency = colMeans(scores > 0.5, na.rm = TRUE) 
save(srt, file = 'output/Fig4/srt_scored_subtype.RData')

scores$sample <- paste0("S",str_sub(row.names(scores),1,3))
scores$state_TCGA <- srt$state_TCGA
scores$cluster <- srt$seurat_clusters
saveRDS(scores,file = "output/Fig4/scores_subtype.RDS")

plot.list = lapply(names(modules), function(m){
  h = FeaturePlot(srt, features = m, pt.size = 0.1, cols = c('grey',colors.module[m])) +
    NoAxes() + NoLegend() +
    theme(plot.title = element_text(size = 10))
})

pdf("output/Fig4/Fig4F.pdf",width = 6, height = 2)
print(CombinePlots(plot.list, ncol = 3)) 
dev.off()

h = DimPlot(srt, group.by = 'state_TCGA', cols = colors.module)
print(h)


###### 2) TCGA subtype score for each sample ######
###### Method 1: scRNA-seq scoring (using 100 gene random backgroud) ######
scores <- readRDS("output/Fig4/scores_subtype.RDS")

data <- scores[,-c(5,6)]
scores_sc <- reshape2::melt(data,id.vars = 'sample',variable.name = "subtype",value.name = "scores")

# heatmap
scores_sc <- scores_sc%>%dplyr::group_by(sample,subtype)%>%
  dplyr::summarise(avg=mean(scores))
scores_sc <- dcast(scores_sc, sample~subtype)
rownames(scores_sc)=scores_sc$sample
scores_sc <- scores_sc[,-1]
scores_sc <- scores_sc[c("S58T","S42T","S51T","S56L","S37T","S46T","S49T","S56T"),]
p3.1 <- ComplexHeatmap::pheatmap(t(scores_sc),cluster_rows = F,cluster_cols = F,name = "scores",main = "scRNA-seq scores")


###### Method 2: scRNA-seq gsva ######
# get single cell expression matrix
# srat_harmony <- readRDS("data/srat_harmony_dims50_res0.5.RDS")
Idents(srat_harmony) <- "seq_folder"
expr <- AverageExpression(srat_harmony,assays = "RNA",slot = "data")[[1]]
expr <- expr[rowSums(expr)>0.1,] 
expr_sc <- as.matrix(expr)
colnames(expr_sc) <- substr(colnames(expr_sc),1,4)
expr_sc <- expr_sc[,c(1,2,3,4,5,8,7,6)]

# TCGA subtype
TCGA_HPV <- readRDS("data/TCGA_HPV.RDS")
TCGA_HPV=TCGA_HPV[-1]

# GSVA analysis
scRNA_gsva <- gsva(expr_sc, TCGA_HPV, kcdf="Gaussian",method = "gsva",parallel.sz=10)
scRNA_gsva <- scRNA_gsva[,c("S58T","S42T","S51T","S56L","S37T","S46T","S49T","S56T")]

p3.2 <- ComplexHeatmap::pheatmap(scRNA_gsva,cluster_rows = F,cluster_cols = F,name = "scores",main = "scRNA-seq (GSVA)")


###### Method 3: bulk RNA-seq gsva ######
# get bulk RNA-seq expression matrix
exprSet<-read.csv(file = 'data/bulk_expr.csv',header = T,row.names = 1) 
exprSet<-exprSet[,(1:13)]
colnames(exprSet) <- c("S04TO","S06TO","S36TO","S37TO","S41TO","S42TO","S46TO","S49TO","S51TO","S55TO","S56TO","S56LNO","S58TO")
exprSet <- exprSet[,c("S37TO","S42TO","S46TO","S49TO","S51TO","S56LNO","S56TO","S58TO")]
colnames(exprSet) <- substr(colnames(exprSet),1,4)
exprSet <- exprSet[rowMeans(exprSet)>1,] 
expr_bulk <- as.matrix(exprSet) 

# TCGA subtype
TCGA_HPV <- readRDS("data/TCGA_HPV.RDS")
TCGA_HPV=TCGA_HPV[-1]

# GSVA analysis
bulk_gsva <- gsva(expr_bulk, TCGA_HPV, kcdf="Gaussian",method = "gsva",parallel.sz=10) 
bulk_gsva <- bulk_gsva[,c("S58T","S42T","S51T","S56L","S37T","S46T","S49T","S56T")]

p3.3 <- ComplexHeatmap::pheatmap(bulk_gsva,cluster_rows = F,cluster_cols = F,name = "scores",main = "bulk RNA-seq (GSVA)")

pdf("output/Fig4/Fig4E.pdf", height = 3, width = 6)
plot(p3.2)
plot(p3.3)
plot(p3.1)
dev.off()

gsva <- cbind(t(scRNA_gsva),t(bulk_gsva),scores_sc)
colnames(gsva) <- c("Basal_scRNA","Atypical_scRNA","Classical_scRNA","Basal_bulk","Atypical_bulk","Classical_bulk","Basal_scoring","Atypical_scoring","Classical_scoring")
gsva$sample <- rownames(gsva)
write.csv(gsva,file = "output/Fig4/bulk_vs_scRNA_gsva.csv")


###### 3) Comparing the three methods ######
###### scRNA_gsva vs bulk_gsva ######
p1.1 <- ggplot(gsva,aes(x=Basal_bulk,y=Basal_scRNA))+ geom_point(size=2,shape=15)+geom_smooth(method=lm)+stat_cor(method = "pearson",label.x = min(gsva$Basal_bulk),label.y = 1,size=5)+xlab("bulk RNA-seq")+ylab("scRNA-seq")+ggtitle("Basal")+my_theme+theme(plot.title = element_text(size = 28, face = "bold",hjust =0.4))+ guides(colour=guide_legend(title=NULL))
p1.2 <- ggplot(gsva,aes(x=Atypical_bulk,y=Atypical_scRNA))+ geom_point(size=2,shape=15)+geom_smooth(method=lm)+stat_cor(method = "pearson",label.x = min(gsva$Atypical_bulk),label.y = max(gsva$Atypical_bulk),size=5)+xlab("bulk RNA-seq")+ylab("scRNA-seq")+ggtitle("Atypical")+my_theme+theme(plot.title = element_text(size = 28, face = "bold",hjust =0.4))+ guides(colour=guide_legend(title=NULL))
p1.3 <- ggplot(gsva,aes(x=Classical_bulk,y=Classical_scRNA))+ geom_point(size=2,shape=15)+geom_smooth(method=lm)+stat_cor(method = "pearson",label.x = min(gsva$Classical_bulk),label.y = 0.6,size=5)+xlab("bulk RNA-seq")+ylab("scRNA-seq")+ggtitle("Classical")+my_theme+theme(plot.title = element_text(size = 28, face = "bold",hjust =0.4))+ guides(colour=guide_legend(title=NULL))
pdf("output/Fig4/Fig4G.pdf",width = 16,height = 5)
p1.1+p1.2+p1.3
dev.off()

###### scRNA_scoring vs bulk_gsva ######
p2.1 <- ggplot(gsva,aes(x=Basal_bulk,y=Basal_scoring))+ geom_point(size=2,shape=15)+geom_smooth(method=lm)+stat_cor(method = "pearson",label.x = min(gsva$Basal_bulk),label.y = 1.1,size=5)+xlab("bulk RNA-seq")+ylab("scRNA-seq scores")+ggtitle("Basal")+my_theme+theme(plot.title = element_text(size = 28, face = "bold",hjust =0.4))+ guides(colour=guide_legend(title=NULL))
p2.2 <- ggplot(gsva,aes(x=Atypical_bulk,y=Atypical_scoring))+ geom_point(size=2,shape=15)+geom_smooth(method=lm)+stat_cor(method = "pearson",label.x = min(gsva$Atypical_bulk),label.y = 0.4,size=5)+xlab("bulk RNA-seq")+ylab("scRNA-seq scores")+ggtitle("Atypical")+my_theme+theme(plot.title = element_text(size = 28, face = "bold",hjust =0.4))+ guides(colour=guide_legend(title=NULL))
p2.3 <- ggplot(gsva,aes(x=Classical_bulk,y=Classical_scoring))+ geom_point(size=2,shape=15)+geom_smooth(method=lm)+stat_cor(method = "pearson",label.x = min(gsva$Classical_bulk),label.y = 0.6,size=5)+xlab("bulk RNA-seq")+ylab("scRNA-seq scores")+ggtitle("Classical")+my_theme+theme(plot.title = element_text(size = 28, face = "bold",hjust =0.4))+ guides(colour=guide_legend(title=NULL))
p2=p2.1+p2.2+p2.3

###### scRNA_scoring vs scRNA_gsva ######
p3.1 <- ggplot(gsva,aes(x=Basal_scoring,y=Basal_scRNA))+ geom_point(size=2,shape=15)+geom_smooth(method=lm)+stat_cor(method = "pearson",label.x = min(gsva$Basal_scoring),label.y = 1,size=5)+xlab("scRNA-seq scores")+ylab("scRNA-seq")+ggtitle("Basal")+my_theme+theme(plot.title = element_text(size = 28, face = "bold",hjust =0.4))+ guides(colour=guide_legend(title=NULL))
p3.2 <- ggplot(gsva,aes(x=Atypical_scoring,y=Atypical_scRNA))+ geom_point(size=2,shape=15)+geom_smooth(method=lm)+stat_cor(method = "pearson",label.x = min(gsva$Atypical_scoring),label.y = 0.6,size=5)+xlab("scRNA-seq scores")+ylab("scRNA-seq")+ggtitle("Atypical")+my_theme+theme(plot.title = element_text(size = 28, face = "bold",hjust =0.4))+ guides(colour=guide_legend(title=NULL))
p3.3 <- ggplot(gsva,aes(x=Classical_scoring,y=Classical_scRNA))+ geom_point(size=2,shape=15)+geom_smooth(method=lm)+stat_cor(method = "pearson",label.x = min(gsva$Classical_scoring),label.y = 0.7,size=5)+xlab("scRNA-seq scores")+ylab("scRNA-seq")+ggtitle("Classical")+my_theme+theme(plot.title = element_text(size = 28, face = "bold",hjust =0.4))+ guides(colour=guide_legend(title=NULL))
pdf("output/Fig4/Fig4H.pdf",width = 16,height = 5)
p3.1+p3.2+p3.3
dev.off()


