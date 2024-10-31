# ---------------------------
# set up env and load library
# ---------------------------
setwd('/home/park/storage2/organoid_code')

library(maftools)
library(patchwork)
library(ggplot2)

# -------------
# load maf file
# -------------
hnc.1.maf = data.table::fread("/home/park/storage2/organoid_code/Figure3/input/tumorVsEarly.maf", header = T)
hnc.2.maf = data.table::fread("/home/park/storage2/organoid_code/Figure3/input/earlyVsLate.maf", header = T)

# ------------------
# load clinical info
# ------------------
hnc.1.clin = read.csv("/home/park/storage2/organoid_code/Figure3/input/tumorVsEarly.txt")
hnc.2.clin = read.csv("/home/park/storage2/organoid_code/Figure3/input/earlyVsLate.txt")

# --------------------------
# load tcga pathway gene set
# --------------------------
pathways <- read.table("/home/park/storage2/organoid_code/Figure3/input/TCGA_pathway.txt", header=TRUE, sep="\t")

# -----------------
# making maf object
# -----------------
hnc.1 = read.maf(maf = hnc.1.maf, clinicalData = hnc.1.clin, verbose = FALSE )
hnc.2 = read.maf(maf = hnc.2.maf, clinicalData = hnc.2.clin, verbose = FALSE )

# -----------------
# setting gene list
# -----------------
tcga.genes <- c('CDKN2A', 'FAT1', 'TP53', 'CASP8', 'AJUBA', 'PIK3CA', 'NOTCH1', 'KMT2D', 'NSD1', 'HLA-A',
                'TGFBR2', 'HRAS', 'FBXW7', 'RB1' ,'PIK3R1', 'TRAF3', 'NFE2L2', 'CUL3', 'PTEN')

cna.genes <- c('TP63', 'SOX2', 'PIK3CA', 'FAT1', 'NSD1', 'EGFR', 'FGFR1', 'CDKN2A', 'CDKN2B', 'NOTCH1', 
               'CCND1', 'FADD', 'FGF3', 'FGF4', 'CTTN', 'YAP1', 'SMAD4')

genes <- as.data.frame(c(pathways$Genes, tcga.genes, cna.genes))

colnames(genes)[1] <- "Gene.Symbol"
genes <- genes[!duplicated(genes$Gene.Symbol), ]

genes.list <- genes[genes %in% c(hnc.1.maf$Hugo_Symbol, hnc.2.maf$Hugo_Symbol, cna.genes)]

# ----------------------------
# Tumor vs Early PDO Oncoplots 
# ----------------------------
pdf("oncoplot.tumorVsEarly.pdf", height = 15, width = 10)
oncoplot(maf = hnc.1, showTumorSampleBarcodes = TRUE, removeNonMutated = FALSE,
         sampleOrder = c('S04T', 'S04E', 'S06T', 'S06E', 'S36T', 'S36E', 'S37T', 'S37E',
                         'S41T', 'S41E', 'S42T', 'S42E', 'S46T', 'S46E', 'S49T', 'S49E', 
                         'S51T', 'S51E', 'S55T', 'S55E', 'S56T', 'S56E', 'S56LNT', 'S56LNE', 
                         'S58T', 'S58E'), 
         top = 60, gene_mar = 10, keepGeneOrder = TRUE, genes = c(genes.list))
dev.off()

# ---------------------------
# Early vs Late PDO Oncoplots 
# ---------------------------
pdf('oncoplot.earlyVsLate.pdf', height = 15, width = 10)
oncoplot(maf = hnc.2, showTumorSampleBarcodes = TRUE, removeNonMutated = FALSE,
         sampleOrder = c('S04E', 'S04L', 'S06E', 'S06L', 'S36E', 'S36L', 'S37E', 'S37L',
                         'S41E', 'S41L', 'S42E', 'S42L', 'S46E', 'S46L', 'S49E', 'S49L', 
                         'S51E', 'S51L', 'S55E', 'S55L', 'S56E', 'S56L', 'S56LNE', 'S56LNL', 
                         'S58E', 'S58L'), 
         top = 60, gene_mar = 10, keepGeneOrder = TRUE, genes = c(genes.list))
dev.off()

# -------------------------------------
# Manually combine with CNA results !!!
# -------------------------------------