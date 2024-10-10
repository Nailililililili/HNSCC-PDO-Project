setwd("~/HNSCC_PDO_Project")

##### Module 1. Identifying heterogeneous processes using nmf #####
source("code/nmf_programs.R")

# read the matrix of single cell data after log-normalization
HNSCC_expr <- readRDS("data/HNSCC_expr.RDS")

# perform NMF with ranks ranging from 6 to 9 
w_basis_human <- list() # nmf gene scores
h_coef_human <- list() # nmf cell scores
for(i in names(HNSCC_expr)) {
  w <- NULL
  h <- NULL
  for(j in 6:9) {
    print(paste0(i, " ", j))
    n <- nmf_programs(HNSCC_expr[[i]],is.log= T, rank=j)
    colnames(n$w_basis) <- paste0(i, "_", j, ".", 1:j)
    colnames(n$h_coef) <- paste0(i, "_", j, ".", 1:j)
    w <- cbind(w, n$w_basis)
    h <- cbind(h, n$h_coef)
  }
  w_basis_human[[i]] <- w
  h_coef_human[[i]] <- h
}
# save output
saveRDS(w_basis_human, "output/Fig5/nmf_w_basis_human.RDS")
saveRDS(h_coef_human, "output/Fig5/nmf_h_coef_human.RDS")


##### Module 2. Identifying recurrent heterogeneous programs (RHPs) in PDOs #####
library(reshape2)
library(ggplot2)
library(scales)
library(egg)
source("code/custom_magma.R") 
source("code/robust_nmf_programs.R") 
source("code/nmf_cell_class.R") 

# read the matrix of single cell data after log-normalization
HNSCC_expr <- readRDS("data/HNSCC_expr.RDS")

# read heterogeneity programs identified using NMF
nmf_programs_genes_human <- readRDS("output/Fig5/nmf_w_basis_human.RDS") # nmf gene scores
nmf_programs_cells_human <- readRDS("output/Fig5/nmf_h_coef_human.RDS") # nmf cell scores

# Identifying recurrent heterogeneous programs (RHPs) in PDOs - NMF
# get gene programs (top 50 genes by NMF score)
nmf_programs_sig_human <- lapply(nmf_programs_genes_human, function(x) apply(x, 2, function(y) names(sort(y, decreasing = T))[1:50])) 

# for each POD, select robust NMF programs (i.e. obseved using different ranks in the same PDO), remove redundancy due to multiple ranks, and apply a filter based on the similarity to programs from other PDOs. 
nmf_filter_human <- robust_nmf_programs(nmf_programs_sig_human, intra_min = 35, intra_max = 10, inter_filter=T, inter_min = 10)
nmf_programs_sig_human <- lapply(nmf_programs_sig_human, function(x) x[, is.element(colnames(x), nmf_filter_human),drop=F])
nmf_programs_sig_human <- do.call(cbind, nmf_programs_sig_human) 

# calculate similarity between programs 
nmf_intersect_human <- apply(nmf_programs_sig_human , 2, function(x) apply(nmf_programs_sig_human , 2, function(y) length(intersect(x,y)))) 

# hierarchical clustering of the similarity matrix 
nmf_intersect_hc_human <- hclust(as.dist(50-nmf_intersect_human), method="average") 
nmf_intersect_hc_human <- reorder(as.dendrogram(nmf_intersect_hc_human), colMeans(nmf_intersect_human))
nmf_intersect_human <- nmf_intersect_human[order.dendrogram(nmf_intersect_hc_human), order.dendrogram(nmf_intersect_hc_human)]

# plot similarity matrix heatmap   
nmf_intersect_meltI_human <- reshape2::melt(nmf_intersect_human) 

p1 <- ggplot(data = nmf_intersect_meltI_human, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + 
  geom_tile() + 
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +                                
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), 
         panel.background = element_blank(),  axis.line = element_blank(), 
         axis.text = element_text(size = 8), axis.title = element_text(size = 12), 
         legend.title = element_text(size=11), legend.text = element_text(size = 10), 
         legend.text.align = 0.5, legend.justification = "bottom") + 
  theme(axis.text.x = element_text(angle=90))+
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1))+
  scale_x_discrete(name="\nPrograms", breaks=unique(nmf_intersect_meltI_human$Var1)[seq(10, length(unique(nmf_intersect_meltI_human$Var1)), by=10)], labels= seq(10, length(unique(nmf_intersect_meltI_human$Var1)), by = 10)) + 
  scale_y_discrete(name="\nPrograms", breaks=unique(nmf_intersect_meltI_human$Var2)[seq(10, length(unique(nmf_intersect_meltI_human$Var2)), by=10)], labels= seq(10, length(unique(nmf_intersect_meltI_human$Var2)), by = 10))

pdf("output/Fig5/Fig5A.pdf", height = 6, width = 7.5, onefile=FALSE) 
p1
dev.off()

#save programs genes
saveRDS(nmf_programs_sig_human, "output/Fig5/nmf_programs_sig_human.RDS") 


# manually define metaprograms (genes observed in at least 40% of programs composing the respective metacluster)
nmf_meta1_human <- sort(table(nmf_programs_sig_human[,colnames(nmf_intersect_human)[2:8]])/length(2:8), decreasing=T)
nmf_meta1_human_programs <- colnames(nmf_intersect_human)[2:8]

nmf_meta2_human <- sort(table(nmf_programs_sig_human[,colnames(nmf_intersect_human)[9:16]])/length(9:16), decreasing=T) 
nmf_meta2_human_programs <- colnames(nmf_intersect_human)[9:16]

nmf_meta3_human <- sort(table(nmf_programs_sig_human[,colnames(nmf_intersect_human)[17:20]])/length(17:20), decreasing=T) 
nmf_meta3_human_programs <- colnames(nmf_intersect_human)[17:20]

nmf_meta4_human <- sort(table(nmf_programs_sig_human[,colnames(nmf_intersect_human)[23:27]])/length(23:27), decreasing=T) 
nmf_meta4_human_programs <- colnames(nmf_intersect_human)[23:27]

nmf_meta5_human <- sort(table(nmf_programs_sig_human[,colnames(nmf_intersect_human)[29:32]])/length(29:32), decreasing=T) 
nmf_meta5_human_programs <- colnames(nmf_intersect_human)[29:32]

nmf_meta6_human <- sort(table(nmf_programs_sig_human[,colnames(nmf_intersect_human)[36:40]])/length(36:40), decreasing=T) 
nmf_meta6_human_programs <- colnames(nmf_intersect_human)[36:40]

nmf_meta7_human <- sort(table(nmf_programs_sig_human[,colnames(nmf_intersect_human)[41:44]])/length(41:44), decreasing=T) 
nmf_meta7_human_programs <- colnames(nmf_intersect_human)[41:44]

nmf_meta8_human <- sort(table(nmf_programs_sig_human[,colnames(nmf_intersect_human)[48:50]])/length(48:50), decreasing=T) 
nmf_meta8_human_programs <- colnames(nmf_intersect_human)[48:50]

nmf_meta_all_human <- list(RHP1=nmf_meta1_human, RHP2=nmf_meta2_human, RHP3=nmf_meta3_human, RHP4=nmf_meta4_human, RHP5=nmf_meta5_human, RHP6=nmf_meta6_human, RHP7=nmf_meta7_human, RHP8=nmf_meta8_human)       
nmf_meta_all_human_top40 <- lapply(nmf_meta_all_human, function(x) x[x>=0.40]) 
RHP=lapply(nmf_meta_all_human_top40,names)

saveRDS(nmf_meta_all_human,"output/Fig5/nmf_metaprograms_sig_human.RDS")
saveRDS(nmf_meta_all_human_top40,"output/Fig5/nmf_metaprograms_sigtop40_human.RDS") 
saveRDS(list(RHP1=nmf_meta1_human_programs, RHP2=nmf_meta2_human_programs, RHP3=nmf_meta3_human_programs, RHP4=nmf_meta4_human_programs, RHP5=nmf_meta5_human_programs, RHP6=nmf_meta6_human_programs, RHP7=nmf_meta7_human_programs, RHP8=nmf_meta8_human_programs),"output/Fig5/nmf_metaprograms_programs_human.RDS")                                                                         


# plot NMF gene scores
nmf_genes <- do.call("c",lapply(nmf_meta_all_human_top40, function(x) rev(names(x))))
nmf_genes <- sapply(unlist(unname(lapply(nmf_programs_genes_human, function(x) apply(x, 2, list))), recursive=F)[colnames(nmf_intersect_human)], function(x) unlist(x)[nmf_genes])
rownames(nmf_genes) <- NULL
nmf_genes[is.na(nmf_genes)] <- 0

nmf_genes_melt <- reshape2::melt(nmf_genes)

p2 <- ggplot(data = nmf_genes_melt, aes(x=Var2, y=Var1, fill=value*100)) + 
  geom_tile() + 
  scale_fill_gradient2(limits=c(0,8), low= "#b5cce3", mid=  "white", high = "darkred", midpoint = 4,   oob=squish, name=expression(paste("NMF score (10"^"2", ")", sep = ""))) +
  scale_color_gradient2(limits=c(0,8), low= "#b5cce3", mid=  "white", high = "darkred", midpoint = 4,   oob=squish, name=expression(paste("NMF score (10"^"2", ")", sep = ""))) + 
  theme( axis.ticks = element_blank(), panel.background = element_blank(), panel.border=element_rect(colour = "black", size = 0.4, fill=F),   axis.line = element_blank(), axis.title = element_text(size = 12), axis.text = element_text(size = 10), legend.title = element_text(size=10), legend.text = element_text(size = 9), legend.text.align = 0.5, legend.direction = "horizontal",  plot.margin = unit(c(0,3,0.5,0.5), "cm")) +
  scale_x_discrete(name="\nPrograms", breaks=unique(nmf_intersect_meltI_human$Var1)[seq(10, length(unique(nmf_intersect_meltI_human$Var1)), by=10)], labels= seq(10, length(unique(nmf_intersect_meltI_human$Var1)), by = 10)) +
  labs(x="\nPrograms", y="Genes") +
  theme(axis.text.x = element_text(angle=90))+
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5, barwidth = 4.4, barheight = 0.6)) 

p1=p1+theme(axis.text.x = element_blank())
# combine plots
pdf("output/Fig5/Fig5C.pdf", height = 9.5, width = 8, onefile=FALSE)
p <- egg::ggarrange( p1,p2, nrow = 2, heights = c(9,7))
dev.off()                               



##### Module 3. Comparing RHPs found in SCC tumor samples #####
library(magrittr)
library(tibble)
library(tidyr)

# read the matrix of single cell data after log-normalization
HNSCC_expr <- readRDS("data/HNSCC_expr.RDS")

# calculate average gene expression 
common_genes  <- Reduce(intersect, lapply(HNSCC_expr, rownames))
HNSCC_expr <-  lapply(HNSCC_expr, function(x) x[common_genes,]) 

ave_expr_human <- rowMeans(do.call(cbind, HNSCC_expr))
gene_uni_human <- names(sort(ave_expr_human, decreasing = T)[1:7000]) 

# read  metaprograms from SCC tumors
HNSCC <- readRDS("data/RHP_reference.RDS")
ESCC <- readRDS("data/epithelial_program.RDS")
CSCC <- readRDS("data/CSCC_MPlist.rds")

# read  metaprograms from PODs
RHP <- readRDS("output/Fig5/nmf_metaprograms_sigtop40_human.RDS")
RHP <- lapply(RHP, names)
RHP <- RHP[-c(3,4)]

# Compare metaprograms signatures from cell lines and tumors
# (jaccard index and hypergeometric test)                             
# hypergeometric test and jaccard index
RHP_filt <- lapply(RHP, function(x){ x[is.element(x, gene_uni_human)]}) 
HNSCC_filt <- lapply(HNSCC, function(x){ x[is.element(x, gene_uni_human)]})
ESCC_filt <- lapply(ESCC, function(x){ x[is.element(x, gene_uni_human)]})
CSCC_filt <- lapply(CSCC, function(x){ x[is.element(x, gene_uni_human)]})


###### RHP vs HNSCC ######
meta_jaccard <- sapply(RHP_filt, function(x) sapply(HNSCC_filt, function(y) length(intersect(x,y))/length(union(x,y))))             
meta_phyper  <- sapply(RHP_filt, function(x) sapply(HNSCC_filt, function(y) phyper(q=length(intersect(x,y)), m= length(x), 7000-length(x), k=length(y), lower.tail = F)))   
meta_phyper  <- apply(meta_phyper, 2, function(x) p.adjust(x, n=length(meta_phyper), "fdr"))                                     

# remove tumor metaprograms with low similarity to PDO metaprograms
meta_jaccard <- meta_jaccard[apply(meta_jaccard, 1, function(x) length(which(x > 0.05)))!=0,]                      
meta_phyper <- meta_phyper[rownames(meta_jaccard), colnames(meta_jaccard)]

df <- melt(meta_jaccard) %>% 
  mutate(pvalue = melt(meta_phyper)[, 3],
         p_signif = symnum(pvalue, corr = FALSE, na = FALSE,  
                           cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                           symbols = c("***", "**", "*", "", " "))) %>% 
  set_colnames(c("Ref", "RHP", "r", "p", "p_signif"))

# plot heatmap                                     
meta_jaccard_melt <- reshape2::melt(meta_jaccard)                      
red_palette <- brewer.pal(9, "Reds")   
df$value=meta_jaccard_melt$value

pdf("output/Fig5/Fig5E.pdf", width = 6, height = 4) 
ggplot(df, aes(x=RHP, y=Ref, fill=value, color=value)) +
  geom_tile() + 
  geom_text(label = as.character(df$p_signif), color = "black") +
  labs(x="", y="") +                                 
  scale_fill_gradient2(limits=c(0.01, 0.16), midpoint = 0.085, low= c( "white",red_palette[1:3]), mid= red_palette[4:6], high = red_palette[7:9] ,   oob=squish, name="Jaccard\nindex")+
  scale_color_gradient2(limits=c(0.01, 0.16), midpoint = 0.085, low= c("white",red_palette[1:3]), mid= red_palette[4:6], high = red_palette[7:9] ,   oob=squish, name="Jaccard\nindex") +
  theme(axis.text.x=element_text(angle=45, hjust=1), panel.background = element_blank(), panel.border = element_rect(fill=F), axis.text=element_text(size=11), axis.title=element_text(size=12))
dev.off()    


###### RHP vs ESCC ######
ESCC_filt$Epi2.0 <- NULL

meta_jaccard <- sapply(RHP_filt, function(x) sapply(ESCC_filt, function(y) length(intersect(x,y))/length(union(x,y))))             
meta_phyper  <- sapply(RHP_filt, function(x) sapply(ESCC_filt, function(y) phyper(q=length(intersect(x,y)), m= length(x), 7000-length(x), k=length(y), lower.tail = F)))   
meta_phyper  <- apply(meta_phyper, 2, function(x) p.adjust(x, n=length(meta_phyper), "fdr"))                                     

# remove tumor metaprograms with low similarity to PDO metaprograms
meta_jaccard <- meta_jaccard[apply(meta_jaccard, 1, function(x) length(which(x > 0.05)))!=0,]                      
meta_phyper <- meta_phyper[rownames(meta_jaccard), colnames(meta_jaccard)]

df <- melt(meta_jaccard) %>% 
  mutate(pvalue = melt(meta_phyper)[, 3],
         p_signif = symnum(pvalue, corr = FALSE, na = FALSE,  
                           cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                           symbols = c("***", "**", "*", "", " "))) %>% 
  set_colnames(c("Ref", "RHP", "r", "p", "p_signif"))

# plot heatmap                                     
meta_jaccard_melt <- reshape2::melt(meta_jaccard)                      
red_palette <- brewer.pal(9, "Reds")   
df$value=meta_jaccard_melt$value

pdf("output/Fig5/Fig5F.pdf", width = 6, height = 4) 
ggplot(df, aes(x=RHP, y=Ref, fill=value, color=value)) +
  geom_tile() + 
  geom_text(label = as.character(df$p_signif), color = "black") +
  labs(x="", y="") +                                 
  scale_fill_gradient2(limits=c(0.01, 0.16), midpoint = 0.085, low= c( "white",red_palette[1:3]), mid= red_palette[4:6], high = red_palette[7:9] ,   oob=squish, name="Jaccard\nindex")+
  scale_color_gradient2(limits=c(0.01, 0.16), midpoint = 0.085, low= c("white",red_palette[1:3]), mid= red_palette[4:6], high = red_palette[7:9] ,   oob=squish, name="Jaccard\nindex") +
  theme(axis.text.x=element_text(angle=45, hjust=1), panel.background = element_blank(), panel.border = element_rect(fill=F), axis.text=element_text(size=11), axis.title=element_text(size=12))
dev.off()    


###### RHP vs CSCC ######
meta_jaccard <- sapply(RHP_filt, function(x) sapply(CSCC_filt, function(y) length(intersect(x,y))/length(union(x,y))))             
meta_phyper  <- sapply(RHP_filt, function(x) sapply(CSCC_filt, function(y) phyper(q=length(intersect(x,y)), m= length(x), 7000-length(x), k=length(y), lower.tail = F)))   
meta_phyper  <- apply(meta_phyper, 2, function(x) p.adjust(x, n=length(meta_phyper), "fdr"))                                     

# remove tumor metaprograms with low similarity to PDO metaprograms
meta_jaccard <- meta_jaccard[apply(meta_jaccard, 1, function(x) length(which(x > 0.05)))!=0,]                      
meta_phyper <- meta_phyper[rownames(meta_jaccard), colnames(meta_jaccard)]

df <- melt(meta_jaccard) %>% 
  mutate(pvalue = melt(meta_phyper)[, 3],
         p_signif = symnum(pvalue, corr = FALSE, na = FALSE,  
                           cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                           symbols = c("***", "**", "*", "", " "))) %>% 
  set_colnames(c("Ref", "RHP", "r", "p", "p_signif"))

# plot heatmap                                     
meta_jaccard_melt <- reshape2::melt(meta_jaccard)                      
red_palette <- brewer.pal(9, "Reds")   
df$value=meta_jaccard_melt$value

pdf("output/Fig5/Fig5G.pdf", width = 6, height = 4) 
ggplot(df, aes(x=RHP, y=Ref, fill=value, color=value)) +
  geom_tile() + 
  geom_text(label = as.character(df$p_signif), color = "black") +
  labs(x="", y="") +                                 
  scale_fill_gradient2(limits=c(0.01, 0.28), midpoint = 0.145, low= c( "white",red_palette[1:3]), mid= red_palette[4:6], high = red_palette[7:9] ,   oob=squish, name="Jaccard\nindex")+
  scale_color_gradient2(limits=c(0.01, 0.28), midpoint = 0.145, low= c("white",red_palette[1:3]), mid= red_palette[4:6], high = red_palette[7:9] ,   oob=squish, name="Jaccard\nindex") +
  theme(axis.text.x=element_text(angle=45, hjust=1), panel.background = element_blank(), panel.border = element_rect(fill=F), axis.text=element_text(size=11), axis.title=element_text(size=12))
dev.off()    


###### Compare individual programs from PDOs to tumor metaprograms #####
library(reshape2)
library(ggplot2)
library(scales)
source("code/custom_magma.R") 
source("code/control_geneset.R") 

# read the matrix of single cell data after log-normalization
HNSCC_expr <- readRDS("data/HNSCC_expr.RDS")

# calculate average gene expression 
common_genes  <- Reduce(intersect, lapply(HNSCC_expr, rownames))
HNSCC_expr <-  lapply(HNSCC_expr, function(x) x[common_genes,]) 

ave_expr_human <- rowMeans(do.call(cbind, HNSCC_expr))
gene_uni_human <- names(sort(ave_expr_human, decreasing = T)[1:7000])

# read  metaprograms from tumors
MP <- readRDS("data/MPgenelist.RDS")

# read  metaprograms from PDOs
RHP=readRDS("output/Fig5/nmf_metaprograms_sigtop40_human.RDS")
RHP=lapply(RHP, names)

# read nmf programs from PDOs
nmf_programs_sig_human <- readRDS("output/Fig5/nmf_programs_sig_human.RDS")
nmf_meta_programs_human <- readRDS("output/Fig5/nmf_metaprograms_programs_human.RDS")


# Compare individual programs (not metaprograms) from PDOs to tumor metaprograms 
# (mean jaccard index and mean correlation of program scores)               
# calculate the jaccard index between nmf programs from cell lines and metaprograms from tumors
indprog_jaccard_RHPvsMP <- sapply(MP, function(x) apply(nmf_programs_sig_human[,unlist(nmf_meta_programs_human)], 2, function(y) length(intersect(x,y))/length(union(x,y))))

# calculate program scores
nmf_meta_sig_RHP_filtII <- apply(nmf_programs_sig_human[,unlist(nmf_meta_programs_human)], 2, function(x) x[is.element(x, common_genes)])
nmf_meta_sig_RHP_filtII_ctr <- lapply(nmf_meta_sig_RHP_filtII, function(x) control_geneset(ave_tpm = ave_expr_human, program = x, bins = 42))

scores_RHP <- list()
for(i in names(nmf_meta_sig_RHP_filtII)) {
  a <- gsub(".{4}$", "", i)
  b <- HNSCC_expr[[a]]
  scores_RHP[[i]] <- colMeans(b[ nmf_meta_sig_RHP_filtII[[i]],]) - colMeans(b[nmf_meta_sig_RHP_filtII_ctr[[i]],])
}

# calculate program scores - MP metaprograms
meta_sig_MP_all_filtII <- lapply(MP, function(x) x[is.element(x, common_genes)])
meta_sig_MP_all_filtII_ctr <- lapply(meta_sig_MP_all_filtII, function(x) control_geneset(ave_tpm = ave_expr_human, program = x, bins = 42))

scores_MP <- list()
for(i in names(nmf_meta_sig_RHP_filtII)) {
  a <- gsub(".{4}$", "", i)
  b <- HNSCC_expr[[a]]
  scores_MP[[i]] <- sapply(meta_sig_MP_all_filtII, function(x) colMeans(b[x,])) - sapply(meta_sig_MP_all_filtII_ctr, function(x) colMeans(b[x,]))
}

# calculate correlation between program scores - RHP vs MP 
indprog_corr_RHPvsMP <- sapply(unlist(nmf_meta_programs_human), function(x) cor(scores_RHP[[x]], scores_MP[[x]]))
rownames(indprog_corr_RHPvsMP) <- names(MP)                       
colnames(indprog_corr_RHPvsMP) <- unlist(nmf_meta_programs_human)

# calculating random correlation and jaccard index 
jaccard_perm <- list()
corr_perm <- list()
for(i in 1:100) {
  print(i)
  a <- lapply(nmf_meta_sig_RHP_filtII, function(x) control_geneset(ave_tpm = ave_expr_human, program = x, bins = 42, size=1, seed = i))
  a_crt <- lapply(a, function(x) control_geneset(ave_tpm = ave_expr_human, program = x, bins = 42))     
  b <- list()
  for(j in names(a)) {
    c <- gsub(".{4}$", "", j)
    d <- HNSCC_expr[[c]]
    b[[j]] <- colMeans(d[a[[j]],]) - colMeans(d[a_crt[[j]],])
  }
  jaccard_perm[[i]] <- sapply(MP, function(x) sapply(a, function(y) length(intersect(x,y))/length(union(x,y))))
  corr_perm[[i]] <- sapply(names(a), function(x) cor(b[[x]], scores_MP[[x]]))                                   
}

# aggregate results by metaprogram 
indprog_corr_RHPvsMP <- data.frame(aggregate(t(indprog_corr_RHPvsMP), list(melt(nmf_meta_programs_human)$L1), mean), row.names=1) 
indprog_jaccard_RHPvsMP <- data.frame(aggregate(indprog_jaccard_RHPvsMP, list(melt(nmf_meta_programs_human)$L1), mean), row.names=1) 

# plot
jaccard_RHPvsMP_plot <- data.frame(melt(as.matrix(indprog_corr_RHPvsMP)),  melt(as.matrix(indprog_jaccard_RHPvsMP)))                 

library(ggplot2)
library(ggrepel)
pdf("output/Fig5/Fig5B.pdf", width = 12, height = 8, onefile = F)                                                                                                
ggplot(jaccard_RHPvsMP_plot, aes(x=value, y=value.1)) +
  geom_point(size=4, shape=21, fill="gray90") +
  facet_wrap(facets = vars(Var1), nrow = 2) +
  geom_label_repel(  
    data = subset(jaccard_RHPvsMP_plot, value >= 0.65 & abs(value.1) >= 0.032),
    aes(label = Var2),
    size = 3.25,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")) +
  geom_hline(yintercept =  quantile(unlist(jaccard_perm), 0.999), linetype="dashed") +
  geom_vline(xintercept =  quantile(unlist(corr_perm), 0.999), linetype="dashed") +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=F), strip.text = element_text(size=13), strip.background = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=13)) +
  labs(y="RHPvsMP Mean similarity\n(Jaccard Index)", x="Mean correlation") +
  scale_x_continuous(breaks=seq(-0.5, 0.5, 0.5))
dev.off()



##### Module 4. Evaluating expression heterogeneity in HNSCC by analyzing human tumors and PDOs as a single dataset #####
library(reshape2)
library(ggplot2)
library(scales)

# read scRNA-seq data from PDOs and tumors
HNSCC_expr <- readRDS("data/HNSCC_expr.RDS")
expr_tumor <- readRDS("data/CCLE_heterogeneity_Rfiles/tumors_scRNAseq_logTPM.rds") # hnscc primary tumors

# process data                        
expr_tumor <- lapply(expr_tumor, function(x) {(x-rowMeans(x))/sd(x)})
HNSCC_expr <- lapply(HNSCC_expr, function(x) {(x-rowMeans(x))/sd(x)})

rowmeans_tumor <- lapply(expr_tumor, rowMeans)
rowmeans_human <- lapply(HNSCC_expr, rowMeans)

# select relevant PDOs and tumors 
nmf_metaprograms_programs_nc_tumor <- readRDS("data/CCLE_heterogeneity_Rfiles/nmf_metaprograms_programs_nc_tumor.rds") # programs in each tumor metaprogram                             
nmf_metaprograms_programs_human <- readRDS("output/Fig5/nmf_metaprograms_programs_human.RDS")

episen_emt_tumor <- unique(intersect(gsub(".{4}$", "",nmf_metaprograms_programs_nc_tumor$episen), gsub(".{4}$", "",nmf_metaprograms_programs_nc_tumor$emt))) # HNCSS tumors harboring both EpiSen and EMTII    
episen_emt_human=unique(intersect(gsub(".{4}$", "",nmf_metaprograms_programs_human$RHP7), gsub(".{4}$", "",nmf_metaprograms_programs_human$RHP8)))

# read metaprogram signatures from tumors                                          
meta_sig_tumor  <- unlist(list(read.table("data/CCLE_heterogeneity_Rfiles/metaprograms_tumors_literature.txt", sep = "\t", header = T,stringsAsFactors = F)), recursive=F)
meta_sig_tumor <- meta_sig_tumor[c("HNSCC.PEMT","HNSCC.Epidif.1" )]  
names(meta_sig_tumor) <- paste0(names(meta_sig_tumor), "_vivo")
meta_sig_tumor <- lapply(meta_sig_tumor, function(x) x[1:50]) 

# read metaprogram signatures from PDOs
meta_sig_human <- readRDS("output/Fig5/nmf_metaprograms_sigtop40_human.RDS")
meta_sig_human <- lapply(meta_sig_human, names)


# **************************************************************************
# Combine datasets and run PCA - HNSCC
# get expression data from selected PDOs and primary tumor
expr_emt_episen_vivo <- expr_tumor[episen_emt_tumor]  
expr_emt_episen_human <- HNSCC_expr[episen_emt_human]                                         

# select common genes in the cell line and tumor datasets
expr_emt_episen_human <- sapply(expr_emt_episen_human, function(x) x[Reduce(intersect, c(lapply(expr_emt_episen_vivo, rownames), lapply(expr_emt_episen_human, rownames))),]) 
expr_emt_episen_vivo <- sapply(expr_emt_episen_vivo, function(x) x[Reduce(intersect, c(lapply(expr_emt_episen_vivo, rownames), lapply(expr_emt_episen_human, rownames))),]) 

# unlist datasets 
expr_emt_episen_vivo <-  do.call(cbind, expr_emt_episen_vivo)                  
expr_emt_episen_human <-  do.call(cbind, expr_emt_episen_human)                  

# select top expressed genes among cell lines and tumors   
rowmeans_tumor_hnscc <- rowMeans(sapply(rowmeans_tumor[episen_emt_tumor], function(x) x[rownames(expr_emt_episen_vivo)]))
rowmeans_human_hnscc <- rowMeans(sapply(rowmeans_human[episen_emt_human], function(x) x[rownames(expr_emt_episen_human)]))

expr_emt_episen_vivo_filt <- expr_emt_episen_vivo[order(rowMeans(cbind(rowmeans_human_hnscc, rowmeans_tumor_hnscc)), decreasing=T)[1:4500]  ,]
expr_emt_episen_human_filt <- expr_emt_episen_human[order(rowMeans(cbind(rowmeans_tumor_hnscc, rowmeans_human_hnscc)), decreasing=T)[1:4500]  ,]

# combine the two datasets                                        
expr_emt_episen_final <- cbind(expr_emt_episen_vivo_filt,expr_emt_episen_human_filt)    

# run pca
expr_emt_episen_pca <- prcomp(t(expr_emt_episen_final))
saveRDS(expr_emt_episen_final,file = "output/Fig5/expr_emt_episen_final1.rds")
saveRDS(expr_emt_episen_pca,file = "output/Fig5/expr_emt_episen_pca1.rds")


# **************************************************************************
# Plot PCA results
# calculate program scores 
###### 1) human organoids RHP gene list ######
# store pca coordinates in  dataframe                        
expr_emt_episen_plot <- data.frame(expr_emt_episen_pca$x[,1:5])
# add  metadata     
expr_emt_episen_plot$type <- rep(c("Primary tumor", "PDO"), c(ncol(expr_emt_episen_vivo), ncol(expr_emt_episen_human)))                                          
expr_emt_episen_plot$type <- factor(expr_emt_episen_plot$type,levels = c("PDO", "Primary tumor"))

expr_emt_episen_plot$hEMT <- colMeans(expr_emt_episen_final[is.element(rownames(expr_emt_episen_final), meta_sig_human$RHP8),])
expr_emt_episen_plot$Episen <- colMeans(expr_emt_episen_final[is.element(rownames(expr_emt_episen_final), meta_sig_human$RHP7),])    

expr_emt_episen_plot$scores <- expr_emt_episen_plot$hEMT-expr_emt_episen_plot$Episen
hist(expr_emt_episen_plot$scores)
expr_emt_episen_plot$scores <- ifelse(expr_emt_episen_plot$scores>=2,2,ifelse(expr_emt_episen_plot$scores<=-2,-2,expr_emt_episen_plot$scores))

# plot pca coordinates   
# Relative Score
p1.1 <- ggplot(expr_emt_episen_plot, aes(x=PC5, y=PC4, color=scores)) +
  geom_point( size=0.5) +
  scale_color_gradient2(limits=c(-2, 2), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "gray95", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F") ,   oob=squish, name="Relative score\n(hEMT-Episen)") +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=F), legend.text.align = 0.5, axis.text=element_text(size=12), axis.title = element_text(size=13), legend.title.align = 0.5, legend.position = "top", legend.title = element_text(size=12, vjust = 1), legend.text = element_text(size=11), plot.margin=unit(c(0.1,2,0.1,0.1), "cm"))

p1.2 <- ggplot(expr_emt_episen_plot[expr_emt_episen_plot$type=="Primary tumor",] , aes(x=PC5, y=PC4, color=scores)) +
  geom_point( size=0.5) +
  scale_color_gradient2(limits=c(-2, 2), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "gray95", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F") ,   oob=squish, name="Relative score\n(hEMT-Episen)") +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=F), legend.text.align = 0.5, axis.text=element_text(size=12), axis.title = element_text(size=13), legend.title.align = 0.5, legend.position = "top", legend.title = element_text(size=12, vjust = 1), legend.text = element_text(size=11), plot.margin=unit(c(0.1,2,0.1,0.1), "cm"))

p1.3 <- ggplot(expr_emt_episen_plot[expr_emt_episen_plot$type=="PDO",] , aes(x=PC5, y=PC4, color=scores)) +
  geom_point( size=0.5) +
  scale_color_gradient2(limits=c(-2, 2), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "gray95", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F") ,   oob=squish, name="Relative score\n(hEMT-Episen)") +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=F), legend.text.align = 0.5, axis.text=element_text(size=12), axis.title = element_text(size=13), legend.title.align = 0.5, legend.position = "top", legend.title = element_text(size=12, vjust = 1), legend.text = element_text(size=11), plot.margin=unit(c(0.1,2,0.1,0.1), "cm"))

p1 <- p1.3+p1.2+p1.1
  

p2.1 <- ggplot(expr_emt_episen_plot[sample(1:nrow(expr_emt_episen_plot)),], aes(x=PC5, y=PC4, color=type)) +
  geom_point( size=0.5) +
  scale_color_manual(values=c("goldenrod", "mediumpurple"), name="") +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=F), legend.text.align = 0.5, axis.text=element_text(size=12), axis.title = element_text(size=13), legend.title.align = 0.5, legend.position = "top", legend.title = element_text(size=13, vjust = 0.5), legend.text = element_text(size=11), legend.key = element_blank(), plot.margin=unit(c(0.1,2,0.1,0.1), "cm")) +
  guides(color=guide_legend(override.aes = list(size=4)))  

p2.2 <- ggplot(expr_emt_episen_plot[expr_emt_episen_plot$type=="Primary tumor",], aes(x=PC5, y=PC4, color=type)) +
  geom_point(size=0.5) +
  scale_color_manual(values=c("mediumpurple"), name="") +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=F), legend.text.align = 0.5, axis.text=element_text(size=12), axis.title = element_text(size=13), legend.title.align = 0.5, legend.position = "top", legend.title = element_text(size=13, vjust = 0.5), legend.text = element_text(size=11), legend.key = element_blank(), plot.margin=unit(c(0.1,2,0.1,0.1), "cm")) +
  guides(color=guide_legend(override.aes = list(size=4)))  

p2.3 <- ggplot(expr_emt_episen_plot[expr_emt_episen_plot$type=="PDO",], aes(x=PC5, y=PC4, color=type)) +
  geom_point( size=0.5) +
  scale_color_manual(values=c("goldenrod", "mediumpurple","lightblue"), name="") +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=F), legend.text.align = 0.5, axis.text=element_text(size=12), axis.title = element_text(size=13), legend.title.align = 0.5, legend.position = "top", legend.title = element_text(size=13, vjust = 0.5), legend.text = element_text(size=11), legend.key = element_blank(), plot.margin=unit(c(0.1,2,0.1,0.1), "cm")) +
  guides(color=guide_legend(override.aes = list(size=4)))  

p2 <- p2.3+p2.2+p2.1

pdf("output/Fig5/Fig5H.pdf", width = 12, height = 4) 
print(p1)
print(p2)
dev.off()


###### 2) Puram in vivo HNSCC RHP gene list ######
# store pca coordinates in  dataframe                        
expr_emt_episen_plot <- data.frame(expr_emt_episen_pca$x[,1:5])
# add  metadata     
expr_emt_episen_plot$type <- rep(c("Primary tumor", "PDO"), c(ncol(expr_emt_episen_vivo), ncol(expr_emt_episen_human)))                                          
expr_emt_episen_plot$type <- factor(expr_emt_episen_plot$type,levels = c("PDO", "Primary tumor"))

expr_emt_episen_plot$PEMT <- colMeans(expr_emt_episen_final[is.element(rownames(expr_emt_episen_final), meta_sig_tumor$HNSCC.PEMT_vivo),])
expr_emt_episen_plot$Epidif=colMeans(expr_emt_episen_final[is.element(rownames(expr_emt_episen_final), meta_sig_tumor$HNSCC.Epidif.1_vivo),])

expr_emt_episen_plot$scores=expr_emt_episen_plot$PEMT-expr_emt_episen_plot$Epidif
hist(expr_emt_episen_plot$scores)
expr_emt_episen_plot$scores=ifelse(expr_emt_episen_plot$scores>=2,2,ifelse(expr_emt_episen_plot$scores<=-2,-2,expr_emt_episen_plot$scores))

# plot pca coordinates     
p3.1 <- ggplot(expr_emt_episen_plot, aes(x=PC5, y=PC4, color=scores)) +
  geom_point( size=0.5) +
  scale_color_gradient2(limits=c(-2, 2), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "gray95", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F") ,   oob=squish, name="Relative score\n(PEMT-Epidif)") +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=F), legend.text.align = 0.5, axis.text=element_text(size=12), axis.title = element_text(size=13), legend.title.align = 0.5, legend.position = "top", legend.title = element_text(size=12, vjust = 1), legend.text = element_text(size=11), plot.margin=unit(c(0.1,2,0.1,0.1), "cm"))

p3.2 <- ggplot(expr_emt_episen_plot[expr_emt_episen_plot$type=="Primary tumor",] , aes(x=PC5, y=PC4, color=scores)) +
  geom_point( size=0.5) +
  scale_color_gradient2(limits=c(-2, 2), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "gray95", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F") ,   oob=squish, name="Relative score\n(PEMT-Epidif)") +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=F), legend.text.align = 0.5, axis.text=element_text(size=12), axis.title = element_text(size=13), legend.title.align = 0.5, legend.position = "top", legend.title = element_text(size=12, vjust = 1), legend.text = element_text(size=11), plot.margin=unit(c(0.1,2,0.1,0.1), "cm"))

p3.3 <- ggplot(expr_emt_episen_plot[expr_emt_episen_plot$type=="PDO",] , aes(x=PC5, y=PC4, color=scores)) +
  geom_point( size=0.5) +
  scale_color_gradient2(limits=c(-2, 2), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "gray95", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F") ,   oob=squish, name="Relative score\n(PEMT-Epidif)") +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=F), legend.text.align = 0.5, axis.text=element_text(size=12), axis.title = element_text(size=13), legend.title.align = 0.5, legend.position = "top", legend.title = element_text(size=12, vjust = 1), legend.text = element_text(size=11), plot.margin=unit(c(0.1,2,0.1,0.1), "cm"))

p3 <- p3.3+p3.2+p3.1

pdf("output/Fig5/Fig5I.pdf", width = 12, height = 4) 
print(p3)
dev.off()



###### 3) shared RHP gene list ######
# store pca coordinates in  dataframe                        
expr_emt_episen_plot <- data.frame(expr_emt_episen_pca$x[,1:5])
# add  metadata     
expr_emt_episen_plot$type <- rep(c("Primary tumor", "PDO"), c(ncol(expr_emt_episen_vivo), ncol(expr_emt_episen_human)))                                          
expr_emt_episen_plot$type <- factor(expr_emt_episen_plot$type,levels = c("PDO", "Primary tumor"))

expr_emt_episen_plot$PEMT <- colMeans(expr_emt_episen_final[is.element(rownames(expr_emt_episen_final), intersect(meta_sig_tumor$HNSCC.PEMT_vivo, meta_sig_human$RHP8)),])
expr_emt_episen_plot$Epidif=colMeans(expr_emt_episen_final[is.element(rownames(expr_emt_episen_final), intersect(meta_sig_tumor$HNSCC.Epidif.1_vivo, meta_sig_human$RHP7)),])

expr_emt_episen_plot$scores=expr_emt_episen_plot$PEMT-expr_emt_episen_plot$Epidif
hist(expr_emt_episen_plot$scores)
expr_emt_episen_plot$scores=ifelse(expr_emt_episen_plot$scores>=3,3,ifelse(expr_emt_episen_plot$scores<=-3,-3,expr_emt_episen_plot$scores))

# plot pca coordinates     
p4.1 <- ggplot(expr_emt_episen_plot, aes(x=PC5, y=PC4, color=scores)) +
  geom_point( size=0.5) +
  scale_color_gradient2(limits=c(-3, 3), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "gray95", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F") ,   oob=squish, name="Relative score\n(h/pEMT-EpiSen/dif.1)") +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=F), legend.text.align = 0.5, axis.text=element_text(size=12), axis.title = element_text(size=13), legend.title.align = 0.5, legend.position = "top", legend.title = element_text(size=12, vjust = 1), legend.text = element_text(size=11), plot.margin=unit(c(0.1,2,0.1,0.1), "cm"))

p4.2 <- ggplot(expr_emt_episen_plot[expr_emt_episen_plot$type=="Primary tumor",] , aes(x=PC5, y=PC4, color=scores)) +
  geom_point( size=0.5) +
  scale_color_gradient2(limits=c(-3, 3), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "gray95", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F") ,   oob=squish, name="Relative score\n(h/pEMT-EpiSen/dif.1)") +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=F), legend.text.align = 0.5, axis.text=element_text(size=12), axis.title = element_text(size=13), legend.title.align = 0.5, legend.position = "top", legend.title = element_text(size=12, vjust = 1), legend.text = element_text(size=11), plot.margin=unit(c(0.1,2,0.1,0.1), "cm"))

p4.3 <- ggplot(expr_emt_episen_plot[expr_emt_episen_plot$type=="PDO",] , aes(x=PC5, y=PC4, color=scores)) +
  geom_point( size=0.5) +
  scale_color_gradient2(limits=c(-3, 3), midpoint = 0, low= c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE"   ), mid= c("#D1E5F0", "gray95", "#FDDBC7"), high = c( "#F4A582", "#D6604D" ,"#B2182B" ,"#67001F") ,   oob=squish, name="Relative score\n(h/pEMT-EpiSen/dif.1)") +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=F), legend.text.align = 0.5, axis.text=element_text(size=12), axis.title = element_text(size=13), legend.title.align = 0.5, legend.position = "top", legend.title = element_text(size=12, vjust = 1), legend.text = element_text(size=11), plot.margin=unit(c(0.1,2,0.1,0.1), "cm"))

p4 <- p4.3+p4.2+p4.1

pdf("output/Fig5/Fig5J.pdf", width = 12, height = 4) 
print(p4)
dev.off()

