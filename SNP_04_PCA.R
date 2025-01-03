# =========================== #
#
# Maerl Whole Genome Re-sequencing Project 2024
#
# SNP Principal Component Analysis
#
# Species:
# Phymatolithon calcareum
# Lithothamnion corallioides
#
# =========================== #

# In RStudio set working directory to the path where this R script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
set.seed(123)

# Load packages
library(vcfR)
library(adegenet)
library(stringr)
library(dplyr)
library(tidyr)
library(tibble)
library(mapmixture)
library(poppr)
library(hierfstat)
library(mmod)
library(pegas)
library(dartR)
library(RColorBrewer)
library(scales)
library(patchwork)


# ----------------- #
# Phymatolithon calcareum: PCA ####
# ----------------- #

# Read in Pcalcareum SNPs
pcal_genind <- vcfR2genind(read.vcfR("outputs/pcalcareum_SNPs_Britain.vcf.gz"))
pcal_genind

# Merge the following samples into the same group (population)
# • St Mawes coarse maerl from 2011 and 2022
pcal_genind$pop <- as.factor(sub("\\_.*", "", indNames(pcal_genind)))
pcal_genind$pop <- as.factor(str_replace(pcal_genind$pop, "Maw22C|Maw11C", "MawC"))
# pcal_genind$pop <- as.factor(str_replace(pcal_genind$pop, "Maw15|Maw22", "Maw"))
pcal_genind$pop <- as.factor(str_replace(pcal_genind$pop, "Gri|AusII", "Aus"))
summary(pcal_genind$pop)

# Colour palette
library(randomcoloR)
palette_pcal <- distinctColorPalette(nPop(pcal_genind))

# Principal component analysis
pca_pcal <- dudi.pca(pcal_genind, scannf = FALSE, nf = 10)

# Visualise PC1 vs. PC2
(pca_pcal1_plt <- scatter_plot(
  dataframe = pca_pcal$li, group_ids = pcal_genind$pop, labels = pcal_genind$pop,
  type = "labels", axes = c(1,2), xlab = "PC", ylab = "PC", size = 3, colours = palette_pcal,
  percent = round(pca_pcal$eig/sum(pca_pcal$eig)*100, 1)
))

# Visualise PC2 vs. PC3
(pca_pcal2_plt <- scatter_plot(
  dataframe = pca_pcal$li, group_ids = pcal_genind$pop, labels = pcal_genind$pop,
  type = "labels", axes = c(2,3), xlab = "PC", ylab = "PC", size = 3, colours = palette_pcal,
  percent = round(pca_pcal$eig/sum(pca_pcal$eig)*100, 1)
))

# Visualise PC3 vs. PC4
(pca_pcal3_plt <- scatter_plot(
  dataframe = pca_pcal$li, group_ids = pcal_genind$pop, labels = pcal_genind$pop,
  type = "labels", axes = c(3,4), xlab = "PC", ylab = "PC", size = 3, colours = palette_pcal,
  percent = round(pca_pcal$eig/sum(pca_pcal$eig)*100, 1)
))

# Visualise PCX vs. PCY
x <- 4; y <- 5
(pca_pcal4_plt <- scatter_plot(
  dataframe = pca_pcal$li, group_ids = pcal_genind$pop, labels = pcal_genind$pop,
  type = "labels", axes = c(x, y), xlab = "PC", ylab = "PC", size = 3, colours = palette_pcal,
  percent = round(pca_pcal$eig/sum(pca_pcal$eig)*100, 1)
))

# Combine plots
pcal_pcas <- wrap_plots(pca_pcal1_plt, pca_pcal2_plt, guides = "collect")&
  theme(legend.position = "top")&
  guides(fill = guide_legend(nrow = 1, override.aes = aes(label = "")))
pcal_pcas

# Export as supporting figure
ggsave(plot = pcal_pcas, filename = "supporting_figures/Figure_SX_PCA_PCAL.png",
       width = 12, height = 7, units = "in", dpi = 600)

# # Principal component analysis excluding coarse maerl
# pcal_genind_exc_coarse <- popsub(pcal_genind, exclude = c("MawC"))
# pca_pcal_exc_coarse <- dudi.pca(pcal_genind_exc_coarse, scannf = FALSE, nf = 5)
# 
# # Visualise PC1 vs. PC2
# (pca_pcal1_exc_coarse <- scatter_plot(
#   dataframe = pca_pcal_exc_coarse$li, group_ids = pcal_genind_exc_coarse$pop, labels = pcal_genind_exc_coarse$pop,
#   type = "labels", axes = c(1,2), xlab = "PC", ylab = "PC", size = 3, colours = palette_pcal,
#   percent = round(pca_pcal_exc_coarse$eig/sum(pca_pcal_exc_coarse$eig)*100, 1)
# ))
# 
# # Visualise PC2 vs. PC3
# (pca_pcal2_exc_coarse <- scatter_plot(
#   dataframe = pca_pcal_exc_coarse$li, group_ids = pcal_genind_exc_coarse$pop, labels = pcal_genind_exc_coarse$pop,
#   type = "labels", axes = c(2,3), xlab = "PC", ylab = "PC", size = 3, colours = palette_pcal,
#   percent = round(pca_pcal_exc_coarse$eig/sum(pca_pcal_exc_coarse$eig)*100, 1)
# ))

# # Principal component analysis without clones
# MLGs_exc_clones <- read.csv("MLGs_excluding_clones.csv")
# pcal_genind_exc_clones <- pcal_genind[ind = MLGs_exc_clones]
# 
# # Combine plots
# pcal_pcas_mlgs <- wrap_plots(pca_pcal3_plt, pca_pcal4_plt, guides = "collect")
# pcal_pcas_mlgs
# ggsave("supporting_figures/Figure_S4.png", width = 15, height = 7, units = "in", dpi = 600)


# ----------------- #
# Lithothamnion corallioides: PCA ####
# ----------------- #

# Read in Lcorallioides SNPs
lcor_genind <- vcfR2genind(read.vcfR("./outputs/lcorallioides_SNPs.vcf.gz"))

# Merge the following samples into the same group (population)
# • St Mawes maerl from 2015 and 2022
# • All samples from St Austell area (St Austell Bay and Little Gribbin)
lcor_genind$pop <- as.factor(sub("\\_.*", "", indNames(lcor_genind)))
lcor_genind$pop <- as.factor(str_replace(lcor_genind$pop, "Maw15|Maw22", "Maw"))
lcor_genind$pop <- as.factor(str_replace(lcor_genind$pop, "Aus[A-Z]+|Gri", "Aus"))
# lcor_genind$pop <- as.factor(str_replace(lcor_genind$pop, "Mil[1-2]*", "Mil"))
summary(lcor_genind$pop)

# Colour palette
library(randomcoloR)
palette_lcor <- distinctColorPalette(nPop(lcor_genind))

# Principal component analysis
pca_lcor <- dudi.pca(lcor_genind, scannf = FALSE, nf = 10)

# Visualise PC1 vs. PC2
(pca_lcor1_plt <- scatter_plot(
  dataframe = pca_lcor$li, group_ids = lcor_genind$pop, labels = lcor_genind$pop,
  type = "labels", axes = c(1,2), xlab = "PC", ylab = "PC", size = 3, colours = palette_lcor,
  percent = round(pca_lcor$eig/sum(pca_lcor$eig)*100, 1)
))

# Visualise PC2 vs. PC3
(pca_lcor2_plt <- scatter_plot(
  dataframe = pca_lcor$li, group_ids = lcor_genind$pop, labels = lcor_genind$pop,
  type = "labels", axes = c(2,3), xlab = "PC", ylab = "PC", size = 3, colours = palette_lcor,
  percent = round(pca_lcor$eig/sum(pca_lcor$eig)*100, 1)
))

# Visualise PC3 vs. PC4
(pca_lcor3_plt <- scatter_plot(
  dataframe = pca_lcor$li, group_ids = lcor_genind$pop, labels = lcor_genind$pop,
  type = "labels", axes = c(3,4), xlab = "PC", ylab = "PC", size = 3, colours = palette_lcor,
  percent = round(pca_lcor$eig/sum(pca_lcor$eig)*100, 1)
))

# Visualise PCX vs. PCY
x <- 4; y <- 5
(pca_lcor4_plt <- scatter_plot(
  dataframe = pca_lcor$li, group_ids = lcor_genind$pop, labels = lcor_genind$pop,
  type = "labels", axes = c(x, y), xlab = "PC", ylab = "PC", size = 3, colours = palette_lcor,
  percent = round(pca_lcor$eig/sum(pca_lcor$eig)*100, 1)
))

# Combine plots
lcor_pcas <- wrap_plots(pca_lcor1_plt, pca_lcor2_plt, pca_lcor3_plt, guides = "collect")&
  theme(legend.position = "top", legend.text = element_text(size = 14))&
  guides(fill = guide_legend(nrow = 1, override.aes = aes(label = "")))
lcor_pcas

# Export as supporting figure
ggsave("supporting_figures/Figure_SX_PCA_LCOR.png", width = 17, height = 7, units = "in", dpi = 600)

# # Principal component analysis without clones
# MLGs_exc_clones <- read.csv("MLGs_excluding_clones.csv")
# lcor_genind_exc_clones <- lcor_genind[ind = MLGs_exc_clones]
# 
# # Combine both plots
# lcor_pcas_mlgs <- wrap_plots(pca_lcor3_plt, pca_lcor4_plt, guides = "collect")
# lcor_pcas_mlgs
# ggsave("supporting_figures/Figure_S4.png", width = 15, height = 7, units = "in", dpi = 600)

