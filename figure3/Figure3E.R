# -----------------------------
# set up env and load libraries
# -----------------------------
setwd("/home/park/storage2/organoid_code/Figure3")

library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(circlize)

# ---------
# load data
# ---------
countData <- read.csv('./input/all_thresholded.by_genes.txt', sep = '\t') # sequenza result files
countData <- countData[, -c(2,3)]
countData <- tibble::column_to_rownames(countData, var = 'Gene.Symbol')

# ---------------------------
# preparing data for plotting
# ---------------------------

# tumor and early
countData <- countData[, c("S04T", "S04E", "S06T", "S06E", "S36T", "S36E", 
                           "S37T", "S37E", "S41T", "S41E", "S42T", "S42E", 
                           "S46T", "S46E", "S49T", "S49E", "S51T", "S51E", 
                           "S55T", "S55E", "S56T", "S56E", "S56LNT", "S56LNE", 
                           "S58T", "S58E")]

# early and late
countData <- countData[, c("S04E", "S04L", "S06E", "S06L", "S36E", "S36L",
                           "S37E", "S37L", "S41E", "S41L", "S42E", "S42L",
                           "S46E", "S46L", "S49E", "S49L", "S51E", "S51L",
                           "S55E", "S55L", "S56E", "S56L", "S56LNE", "S56LNL",
                           "S58E", "S58L")]

colnames(countData) <- gsub('S', '', colnames(countData))

# -------------------------------------
# Plot sample-sample clustering heatmap
# -------------------------------------
corr <- cor(countData)

# --------------------------------------------------
# pearson correlation heatmap with top annotation
# --------------------------------------------------
col_fun = circlize::colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))

pdf("TumorVsEarly.pdf", width = 40, height = 30)
p <- Heatmap(corr, show_row_names = T, show_column_names = T, cluster_rows = T, cluster_columns = T, 
             column_dend_height = unit(50, "mm"), row_dend_width = unit(50, "mm"), col = col_fun,
             name = "corr", column_title_gp = gpar(fontsize = 50), row_title_gp = gpar(fontsize = 20),
             column_names_gp = gpar(fontsize = 85), column_names_rot = 90, column_names_max_height = unit(100, 'mm'), 
             row_names_gp = gpar(fontsize = 85),
             heatmap_legend_param = list(grid_width = unit(0.5, "cm"), grid_height = unit(10, "cm"),
                                         title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 60)))
p = draw(p)
decorate_heatmap_body(heatmap = "corr", code = {
  grid.rect(x = 0.038, y = 0.963, width = 0.076, height = 0.076, gp = gpar(col = "black", fill = NA, lty = 2, lwd = 10))
})
decorate_heatmap_body(heatmap = "corr", code = {
  grid.rect(x = 0.115, y = 0.885, width = 0.076, height = 0.076, gp = gpar(col = "black", fill = NA, lty = 2, lwd = 10))
})
decorate_heatmap_body(heatmap = "corr", code = {
  grid.rect(x = 0.192, y = 0.807, width = 0.076, height = 0.076, gp = gpar(col = "black", fill = NA, lty = 2, lwd = 10))
})
decorate_heatmap_body(heatmap = "corr", code = {
  grid.rect(x = 0.269, y = 0.729, width = 0.076, height = 0.076, gp = gpar(col = "black", fill = NA, lty = 2, lwd = 10))
})
decorate_heatmap_body(heatmap = "corr", code = {
  grid.rect(x = 0.384, y = 0.617, width = 0.076, height = 0.076, gp = gpar(col = "black", fill = NA, lty = 2, lwd = 10))
})
decorate_heatmap_body(heatmap = "corr", code = {
  grid.rect(x = 0.460, y = 0.539, width = 0.076, height = 0.076, gp = gpar(col = "black", fill = NA, lty = 2, lwd = 10))
})
decorate_heatmap_body(heatmap = "corr", code = {
  grid.rect(x = 0.537, y = 0.461, width = 0.076, height = 0.076, gp = gpar(col = "black", fill = NA, lty = 2, lwd = 10))
})
decorate_heatmap_body(heatmap = "corr", code = {
  grid.rect(x = 0.616, y = 0.384, width = 0.076, height = 0.076, gp = gpar(col = "black", fill = NA, lty = 2, lwd = 10))
})
decorate_heatmap_body(heatmap = "corr", code = {
  grid.rect(x = 0.729, y = 0.270, width = 0.076, height = 0.076, gp = gpar(col = "black", fill = NA, lty = 2, lwd = 10))
})
decorate_heatmap_body(heatmap = "corr", code = {
  grid.rect(x = 0.807, y = 0.192, width = 0.076, height = 0.076, gp = gpar(col = "black", fill = NA, lty = 2, lwd = 10))
})
decorate_heatmap_body(heatmap = "corr", code = {
  grid.rect(x = 0.884, y = 0.114, width = 0.076, height = 0.076, gp = gpar(col = "black", fill = NA, lty = 2, lwd = 10))
})
decorate_heatmap_body(heatmap = "corr", code = {
  grid.rect(x = 0.961, y = 0.037, width = 0.076, height = 0.076, gp = gpar(col = "black", fill = NA, lty = 2, lwd = 10))
})
dev.off() 

pdf("EarlyVsLate.pdf", width = 40, height = 30)
p <- Heatmap(corr, show_row_names = T, show_column_names = T, cluster_rows = T, cluster_columns = T, 
        column_dend_height = unit(50, "mm"), row_dend_width = unit(50, "mm"), col = col_fun,
        name = "corr", column_title_gp = gpar(fontsize = 50), row_title_gp = gpar(fontsize = 20),
        column_names_gp = gpar(fontsize = 85), column_names_rot = 90, column_names_max_height = unit(100, 'mm'), 
        row_names_gp = gpar(fontsize = 85),
        heatmap_legend_param = list(grid_width = unit(0.5, "cm"), grid_height = unit(10, "cm"),
                                  title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 60)))
p = draw(p)
decorate_heatmap_body(heatmap = "corr", code = {
  grid.rect(x = 0.038, y = 0.963, width = 0.076, height = 0.076, gp = gpar(col = "black", fill = NA, lty = 2, lwd = 10))
})
decorate_heatmap_body(heatmap = "corr", code = {
  grid.rect(x = 0.115, y = 0.885, width = 0.076, height = 0.076, gp = gpar(col = "black", fill = NA, lty = 2, lwd = 10))
})
decorate_heatmap_body(heatmap = "corr", code = {
  grid.rect(x = 0.192, y = 0.807, width = 0.076, height = 0.076, gp = gpar(col = "black", fill = NA, lty = 2, lwd = 10))
})
decorate_heatmap_body(heatmap = "corr", code = {
  grid.rect(x = 0.269, y = 0.729, width = 0.076, height = 0.076, gp = gpar(col = "black", fill = NA, lty = 2, lwd = 10))
})
decorate_heatmap_body(heatmap = "corr", code = {
  grid.rect(x = 0.346, y = 0.653, width = 0.076, height = 0.076, gp = gpar(col = "black", fill = NA, lty = 2, lwd = 10))
})
decorate_heatmap_body(heatmap = "corr", code = {
  grid.rect(x = 0.423, y = 0.575, width = 0.076, height = 0.076, gp = gpar(col = "black", fill = NA, lty = 2, lwd = 10))
})
decorate_heatmap_body(heatmap = "corr", code = {
  grid.rect(x = 0.500, y = 0.497, width = 0.076, height = 0.076, gp = gpar(col = "black", fill = NA, lty = 2, lwd = 10))
})
decorate_heatmap_body(heatmap = "corr", code = {
  grid.rect(x = 0.577, y = 0.423, width = 0.076, height = 0.076, gp = gpar(col = "black", fill = NA, lty = 2, lwd = 10))
})
decorate_heatmap_body(heatmap = "corr", code = {
  grid.rect(x = 0.654, y = 0.345, width = 0.076, height = 0.076, gp = gpar(col = "black", fill = NA, lty = 2, lwd = 10))
})
decorate_heatmap_body(heatmap = "corr", code = {
  grid.rect(x = 0.729, y = 0.270, width = 0.076, height = 0.076, gp = gpar(col = "black", fill = NA, lty = 2, lwd = 10))
})
decorate_heatmap_body(heatmap = "corr", code = {
  grid.rect(x = 0.807, y = 0.192, width = 0.076, height = 0.076, gp = gpar(col = "black", fill = NA, lty = 2, lwd = 10))
})
decorate_heatmap_body(heatmap = "corr", code = {
  grid.rect(x = 0.884, y = 0.114, width = 0.076, height = 0.076, gp = gpar(col = "black", fill = NA, lty = 2, lwd = 10))
})
decorate_heatmap_body(heatmap = "corr", code = {
  grid.rect(x = 0.961, y = 0.037, width = 0.076, height = 0.076, gp = gpar(col = "black", fill = NA, lty = 2, lwd = 10))
})
dev.off() 
