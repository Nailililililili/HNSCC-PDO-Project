# ---------------------------
# Set up env and load library
# ---------------------------
setwd('/home/park/storage2/organoid_code/Figure1B')

library(maftools)
library(patchwork)
library(ggplot2)

# --------------
# Load MAF files
# --------------
medium.maf = data.table::fread("/home/park/storage2/organoid_code/Figure1B/input/medium.maf", header = T)

# ------------------
# Load clinical info
# ------------------
medium.clin = read.csv("/home/park/storage2/organoid_code/Figure1B/input/medium.csv")

# --------------------------
# Load TCGA pathway gene set
# --------------------------
pathways <- read.table("/home/park/storage2/organoid_code/Figure1B/input/TCGA_pathway.txt", header=TRUE, sep="\t")

# -----------------
# Making MAF object
# -----------------
medium = read.maf(maf = medium.maf, clinicalData = medium.clin, verbose = FALSE )

# -----------------
# Setting gene list
# -----------------
cosmic.genes <- read.csv("/home/park/storage4/project/organoid/04.maftool/MAFTOOLS/input/Census.csv")

cna.genes <- c('TP63', 'SOX2', 'PIK3CA', 'FAT1', 'NSD1', 'EGFR', 'ERBB2', 'FGFR1', 'CDKN2A', 'CDKN2B', 'NOTCH1', 
               'CCND1', 'FADD', 'FGF3', 'FGF4', 'CTTN', 'YAP1', 'SMAD4', 'BIRC2')

genes <- as.data.frame(c(pathways$Genes, cosmic.genes$Gene.Symbol, cna.genes))

colnames(genes)[1] <- "Gene.Symbol"
genes <- genes[!duplicated(genes$Gene.Symbol), ]

genes.list <- genes[genes %in% c(medium.maf$Hugo_Symbol, cna.genes)]

# ---------
# Oncoplots 
# ---------
pdf('oncoplot.S41.medium.pdf', height = 7, width = 2.5)
oncoplot(maf = medium, showTumorSampleBarcodes = TRUE, removeNonMutated = FALSE,
         sampleOrder = c('S41B', 'S41C'), 
         top = 60, gene_mar = 10, genes = c(genes.list, cna.genes),
         keepGeneOrder = TRUE)
dev.off()

pdf('oncoplot.S51.medium.pdf', height = 7, width = 2.7)
oncoplot(maf = medium, showTumorSampleBarcodes = TRUE, removeNonMutated = FALSE,
         sampleOrder = c('S51A', 'S51B', 'S51C'), 
         top = 60, gene_mar = 10, genes = c(genes.list, cna.genes),
         keepGeneOrder = TRUE)
dev.off()

# ------------------------------------
# Manually combine with cna results!!!
# ------------------------------------
