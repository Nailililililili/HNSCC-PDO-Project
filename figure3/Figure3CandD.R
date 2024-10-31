# ------------------
# set up environment
# ------------------
setwd('/home/park/storage2/organoid_code/Figure3')

library(tidyverse)
library(hrbrthemes)
library(tm)
library(proustr)
library(VennDiagram)
library(patchwork)

# ----------------------------
# Common versus tumor plotting
# ----------------------------
dat <- list('early' = c(27, 11, 26, 56, 31, 239, 53, 17, 77, 62, 7, 15, 24),
            'common' = c(28, 63, 47, 57, 69, 2, 74, 84, 55, 23, 48, 43, 91),
            'tumor' = c(10, 6, 5, 20, 8, 7, 22, 43, 5, 26, 6, 5, 7))

dat <- as.data.frame(dat)
dat$ID <- c('04', '06', '36', '37', '41', '42', '46', '49', '51', '55', '56', '56LN', '58')
dat$concordance_total <- dat$common / (dat$tumor + dat$common + dat$early)
dat$concordance_tumor <- dat$common / (dat$tumor + dat$common)

pdf('commonVstumor.pdf')
ggplot(dat, aes(x = reorder(ID, -concordance_tumor), y = concordance_tumor)) + 
  geom_bar(stat = 'identity', width = 0.5) +
  geom_col(mapping = aes(fill = ID), position = 'dodge', width = 0.5) +
  xlab("PDO") +
  ylab("SNVs concordance \n(common / tumor)") +
  theme(axis.title = element_text(size = 15, face = 'bold', color = 'black'),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.text.x = element_text(size = 15, face = 'bold', color = 'black', angle = 90),
        axis.text.y = element_text(size = 15, face = 'bold', color = 'black'),
        axis.line = element_line(linewidth = 0.5),
        legend.position = "none",
        panel.background = element_blank())
dev.off()

# ----------------------------
# Common versus early plotting
# ----------------------------
dat <- list('late' = c(8, 15, 7, 14, 6, 154, 14, 7, 5, 10, 11, 28, 7),
            'common' = c(62, 68, 67, 69, 89, 179, 121, 96, 106, 77, 47, 54, 114),
            'early' = c(3, 6, 6, 44, 11, 62, 6, 5, 26, 8, 8, 4, 1))

dat <- as.data.frame(dat)
dat$ID <- c('04', '06', '36', '37', '41', '42', '46', '49', '51', '55', '56', '56LN', '58')
dat$concordance_total <- dat$common / (dat$early + dat$common + dat$late)
dat$concordance_early <- dat$common / (dat$early + dat$common)

pdf('commonVsearly.pdf')
ggplot(dat, aes(x = reorder(ID, -concordance_early), y = concordance_early)) + 
  geom_bar(stat = 'identity', width = 0.5) +
  geom_col(mapping = aes(fill = ID), position = 'dodge', width = 0.5) +
  xlab("PDO") +
  ylab("SNVs concordance \n(common / early)") +
  theme(axis.title = element_text(size = 15, face = 'bold', color = 'black'),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.text.x = element_text(size = 15, face = 'bold', color = 'black', angle = 90),
        axis.text.y = element_text(size = 15, face = 'bold', color = 'black'),
        axis.line = element_line(linewidth = 0.5),
        legend.position = "none",
        panel.background = element_blank())
dev.off()
