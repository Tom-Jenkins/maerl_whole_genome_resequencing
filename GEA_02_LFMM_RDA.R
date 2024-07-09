# =========================== #
#
# Maerl Whole Genome Resequencing Project 2024
#
# LFMM and RDA: Outlier SNP Detection
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
library(LEA)
library(readr)
library(data.table)
library(ggplot2)
library(gridExtra)
library(vcfR)
library(ade4)
library(vegan)
library(stringr)
library(mapmixture)
library(dplyr)
library(purrr)
library(R.utils)
library(dartR)
library(ggvegan)

# Read in environmental data per site (remove coordinates from data frame)
env_data <- select(read_csv("./outputs/climate_data.csv", show_col_types = F), !contains(c("lat","lon")))
env_data$site <- str_sub(env_data$site, 1, 3)
env_data

# Test removing env data
# env_data <- dplyr::select(env_data, !contains("so"))
env_data <- dplyr::select(env_data, !contains("sws"))
env_data <- dplyr::select(env_data, !contains("phyc"))
# env_data <- dplyr::select(env_data, !contains("ph_"))

# Expected false discovery rate
alpha <- 0.05

# ----------------- #
# Phymatolithon calcareum ####
# ----------------- #

# Export VCF in geno format
# vcf <- read.vcfR("./outputs/pcalcareum_SNPs.vcf.gz")
# dartR::gl2geno(vcfR::vcfR2genlight(vcf), outfile = "pcalcareum_SNPs", outpath = "./outputs/")

# Read in genotypes in LFMM format
geno <- data.table::fread("./outputs/pcalcareum_SNPs.lfmm")

# Convert to individual format
vcf <- read.vcfR("./outputs/pcalcareum_SNPs.vcf.gz")
geno_ind <- tibble(ind = colnames(vcf@gt)[-1], geno)
geno_ind$site <- str_sub(geno_ind$ind, 1, 3)
head(geno_ind)

# Join data
env_data_ind <- left_join(dplyr::select(geno_ind, ind, site), distinct(env_data), by = "site")
env_data_ind

# Baseline environmental data matrix
remove_future <- c("ssp119","ssp245","ssp585")
env_matrix <- dplyr::select(env_data_ind, !contains(remove_future) & !contains("ind") & !contains("site"))

# Re-test multicollinearity
library(psych)
psych::pairs.panels(env_matrix)

#--------------#
# LFMM2 ####
#--------------#

# Run LEA PCA
pc <- LEA::pca("./outputs/pcalcareum_SNPs.lfmm", scale = TRUE)

# Plot eigenvalues
plot(pc, lwd = 5, col = "blue", cex = 0.7, xlab = "Factors", ylab = "Eigenvalues")

# Number of latent factors (K) should be set 3-9
K <- 9

# Genotypes
genotypes <- dplyr::select(geno_ind, !contains("ind") & !contains("site"))

# Run LFMM2 analysis
lfmm_mod <- lfmm2(genotypes, env_matrix, K = K, lambda = 1e-5)

# Compute P-values adjusted by genomic inflation factor
lfmm_test_full <- lfmm2.test(lfmm_mod, genotypes, env_matrix, linear = T, full = T, genomic.control = T)

# Genomic inflation factor
lfmm_test_full$gif

# Adjusted P-values for full model (one per locus)
lfmm_adj_pval <- p.adjust(lfmm_test_full$pvalues, "BH"); hist(lfmm_adj_pval)

# LFMM GEA outliers
lfmm_loci <- which(lfmm_adj_pval < alpha)
length(lfmm_loci)

#--------------#
# RDA ####
#--------------#

# Run PCA
pca1 <- dudi.pca(genotypes, scale = TRUE, scannf = FALSE, nf = 5)
# scatter_plot(as.data.frame(pca1$li), axes = c(1,2), group_ids = geno_ind$site, centroids = F, segments = F)
# scatter_plot(as.data.frame(pca1$li), axes = c(2,3), group_ids = geno_ind$site, centroids = F, segments = F)

# Run RDA and control for population structure
# https://popgen.nescent.org/2018-03-27_RDA_GEA.html
rda1 <- rda(geno ~ . + Condition(pca1$li$Axis1) + Condition(pca1$li$Axis2), data = env_matrix, scale = T)
rda1
RsquareAdj(rda1)

# Variance Inflation Factors
vif.cca(rda1)

# Variance explained by each canonical axis
summary(eigenvals(rda1, model = "constrained"))

# Screeplot (same number of axes as the number of predictors in the model)
screeplot(rda1)

# Check each canonical axis for significance (**very long run time ~30 mins)
# signif_axis <- anova.cca(rda1, by = "axis", parallel = getOption("mc.cores"))
# signif_axis
# Permutation test for rda under reduced model
# Forward tests for axes
# Permutation: free
# Number of permutations: 999
# 
# Model: rda(formula = geno ~ ph_baseline_2000_2018_depthmean + so_baseline_2000_2019_depthmean + sws_baseline_2000_2019_depthmean + thetao_baseline_2000_2019_depthmean, data = env_matrix, scale = T)
# Df Variance       F Pr(>F)    
# RDA1       1   1947.2 19.9423  0.001 ***
#   RDA2       1    803.0  8.2233  0.001 ***
#   RDA3       1    139.8  1.4317  0.355    
# RDA4       1    136.9  1.4026  0.146    
# Residual 126  12303.1                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Extract SNP loadings on constrained axes
rda1_load <- scores(rda1, choices = c(1:2), display = "species")

# Histogram of the loadings on each RDA axis
hist(rda1_load[,1], main="Loadings on RDA1")
hist(rda1_load[,2], main="Loadings on RDA2")

# Function where x is the vector of loadings and z is the number of standard deviations to use
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)    
  x[x < lims[1] | x > lims[2]]
}

# Identify candidate outlier SNPs 
cand1 <- outliers(rda1_load[,1], 3)
cand2 <- outliers(rda1_load[,2], 3)
cand_loci <- unique(names(c(cand1, cand2)))
rda_loci <- as.integer(str_remove(cand_loci, "V"))
length(rda_loci)

#--------------#
# Export Outliers ####
#--------------#

# Outlier loci common among both methods
outlier_loci <- intersect(lfmm_loci, rda_loci)
length(outlier_loci)

# Extract list of candidate loci from vcf
loci_IDs <- paste(vcf@fix[outlier_loci,][,"CHROM"], vcf@fix[outlier_loci,][,"POS"], sep = "_")
loci_IDs

# Export outlier genotypes in VCF format
write.vcf(vcf[outlier_loci,], "./outputs/pcalcareum_GEA_outlier_loci.vcf.gz")

# Export outlier genotypes in lfmm format
geno_outlier <- dplyr::select(geno, all_of(outlier_loci))
write.table(geno_outlier, file = "./outputs/pcalcareum_GEA_outlier_loci.lfmm", row.names = F, col.names = F)

#--------------#
# Outlier RDA ####
#--------------#

# Env matrix with different column names
colnames(env_matrix) <- str_extract(colnames(env_matrix), "^[^_]+")

# Run RDA with only outlier loci
rda_outlier <- rda(geno_outlier ~ ., data = env_matrix, scale = T)

# Total variance explained
RsquareAdj(rda_outlier)

# Create data.frame for plotting
plot_df <- tibble(
  ind = geno_ind$ind,
  site = geno_ind$site
)
plot_df <- cbind(plot_df, env_matrix)

# Change MawC samples from Maw to MawC
plot_df$site[which(str_detect(plot_df$ind, "Maw11C|Maw22C"))] <- "MawC"

# Change site to factor
site_new_order <- c("Zar","Man","Biz","Aus","Gri","Ger","Her","Maw",
                    "MawC","Swa","Wey","Mor","Tre","Bor","Ons")
plot_df$site <- factor(plot_df$site, levels = site_new_order)

# Vector of colours
set.seed(123)
library(randomcoloR)
library(scales)
sample_cols <- c(Aus = "grey60", Biz = "#FF69B4", Bor = "#E17E68",
                 Ger = "#FCCDE5", Gri = "#FAA0A0", Her = "#D9D9D9",
                 Man = "#F3CFC6", Maw = "#FF00FF", MawC = "#FDB462",
                 Mor = "#BEBADA", Ons = "#FFFFB3", Swa = "grey60", 
                 Tre = "#8DD3C7", Wey = "#B3DE69", Zar = "#80B1D3")
# sample_cols <- distinctColorPalette(n_distinct(plot_df$site), runTsne = FALSE)
scales::show_col(sample_cols)

# Order colours
sample_cols <- sample_cols[levels(plot_df$site)]

# Biplot using ordiplot from vegan
# https://popgen.nescent.org/2018-03-27_RDA_GEA.html
ordiplot(rda_outlier, type = "none", choices = c(1,2), scaling = 3) |>
  # points("species", pch = 20, cex = 1, col = "grey32") |>
  text("biplot", col = "black", cex = 1) |>
  points("sites", pch = 21, cex = 1.8, col = "black", bg = sample_cols[plot_df$site])

# Add legend to triplot
legend("bottomright", legend = levels(plot_df$site), bty = "n", col = "black",
       pch = 21, cex = 1, pt.bg = sample_cols)


# # ----------------- #
# # Lithothamnion corallioides ####
# # ----------------- #
# 
# # Export VCF in geno format
# # vcf <- read.vcfR("./outputs/lcorallioides_SNPs.vcf.gz")
# # dartR::gl2geno(vcfR::vcfR2genlight(vcf), outfile = "lcorallioides_SNPs", outpath = "./outputs/")
# 
# # Read in genotypes in LFMM format
# geno <- data.table::fread("./outputs/lcorallioides_SNPs.lfmm")
# 
# # Convert to individual format
# vcf <- read.vcfR("./outputs/lcorallioides_SNPs.vcf.gz")
# geno_ind <- tibble(ind = colnames(vcf@gt)[-1], geno)
# geno_ind$site <- str_sub(geno_ind$ind, 1, 3)
# head(geno_ind)
# 
# # Convert to individual format
# env_data_ind <- left_join(dplyr::select(geno_ind, ind, site), distinct(env_data), by = "site")
# env_data_ind
# 
# # Baseline environmental data matrix
# remove_future <- c("ssp119","ssp245","ssp585")
# env_matrix <- dplyr::select(env_data_ind, !contains(remove_future) & !contains("ind") & !contains("site"))
# 
# # Re-test multicollinearity
# library(psych)
# psych::pairs.panels(env_matrix)
# env_matrix <- dplyr::select(env_matrix, !contains("so"))
# env_matrix <- dplyr::select(env_matrix, !contains("ph"))
# # env_matrix <- dplyr::select(env_matrix, !contains("sws"))
# # psych::pairs.panels(env_matrix)
# 
# #--------------#
# # LFMM2 ####
# #--------------#
# 
# # Run LEA PCA
# pc <- LEA::pca("./outputs/lcorallioides_SNPs.lfmm", scale = TRUE)
# 
# # Plot eigenvalues
# plot(pc, lwd = 5, col = "blue", cex = 0.7, xlab = "Factors", ylab = "Eigenvalues")
# 
# # Number of latent factors (K) should be set to 3-6
# K <- 6
# 
# # Genotypes
# genotypes <- dplyr::select(geno_ind, !contains("ind") & !contains("site"))
# 
# # Run LFMM2 analysis
# lfmm_mod <- lfmm2(genotypes, env_matrix, K = K, lambda = 1e-5)
# 
# # Compute P-values adjusted by genomic inflation factor
# lfmm_test_full <- lfmm2.test(lfmm_mod, genotypes, env_matrix, linear = T, full = T, genomic.control = T)
# 
# # Genomic inflation factor
# lfmm_test_full$gif
# 
# # Adjusted P-values for full model (one per locus)
# lfmm_adj_pval <- p.adjust(lfmm_test_full$pvalues, "BH"); hist(lfmm_adj_pval)
# 
# # LFMM GEA outliers
# lfmm_loci <- which(lfmm_adj_pval < alpha)
# length(lfmm_loci)
# 
# #--------------#
# # RDA ####
# #--------------#
# 
# # Run PCA
# pca1 <- dudi.pca(genotypes, scale = TRUE, scannf = FALSE, nf = 6)
# scatter_plot(as.data.frame(pca1$li), axes = c(1,2), group_ids = geno_ind$site, centroids = F, segments = F)
# # scatter_plot(as.data.frame(pca1$li), axes = c(2,3), group_ids = geno_ind$site, centroids = F, segments = F)
# 
# # Run RDA and control for population structure
# # https://popgen.nescent.org/2018-03-27_RDA_GEA.html
# rda1 <- rda(geno ~ . + Condition(pca1$li$Axis1) + Condition(pca1$li$Axis2), data = env_matrix, scale = T)
# rda1
# RsquareAdj(rda1)
# 
# # Variance Inflation Factors
# vif.cca(rda1)
# 
# # Variance explained by each canonical axis
# summary(eigenvals(rda1, model = "constrained"))
# 
# # Screeplot (same number of axes as the number of predictors in the model)
# screeplot(rda1)
# 
# # Check each canonical axis for significance (**very long run time ~30 mins)
# # signif_axis <- anova.cca(rda1, by = "axis", parallel = getOption("mc.cores"))
# # signif_axis
# 
# # Extract SNP loadings on constrained axes
# rda1_load <- scores(rda1, choices = c(1:2), display = "species")
# 
# # Histogram of the loadings on each RDA axis
# hist(rda1_load[,1], main="Loadings on RDA1")
# hist(rda1_load[,2], main="Loadings on RDA2")
# 
# # Function where x is the vector of loadings and z is the number of standard deviations to use
# outliers <- function(x,z){
#   lims <- mean(x) + c(-1, 1) * z * sd(x)    
#   x[x < lims[1] | x > lims[2]]
# }
# 
# # Identify candidate outlier SNPs 
# cand1 <- outliers(rda1_load[,1], 3)
# cand2 <- outliers(rda1_load[,2], 3)
# cand_loci <- unique(names(c(cand1, cand2)))
# rda_loci <- as.integer(str_remove(cand_loci, "V"))
# length(rda_loci)
# 
# #--------------#
# # Export Outliers ####
# #--------------#
# 
# # Outlier loci common among both methods
# outlier_loci <- intersect(lfmm_loci, rda_loci)
# length(outlier_loci)
