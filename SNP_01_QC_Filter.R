# =========================== #
#
# Maerl Whole Genome Resequencing Project 2024
#
# SNP QC and Filtering Analysis
#
# Species:
# Phymatolithon calcareum
# Lithothamnion corallioides
#
# SNP data files:
# Y.vcf
# X.vcf
#
# =========================== #

# In RStudio set working directory to the path where this R script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load packages
library(vcfR)
library(SNPfiltR)
library(stringr)

# ----------------- #
# Phymatolithon calcareum ####
# ----------------- #

# Read in VCF file
vcf <- read.vcfR("./data/pcalcareum_variants_qc3_diploid.vcf")

# Filter SNPs for linkage
vcf_ld <- distance_thin(vcf, min.distance = 1000)

# Filter SNPs not biallelic
vcf_ld_bi <- filter_biallelic(vcf_ld)

# Extract total read depth per locus
dp <- extract.gt(vcf_ld_bi, element = "DP", as.numeric = TRUE)
queryMETA(vcf_ld_bi, element = "DP")

# Filter minimum and maximum mean depth (over all samples)
minDP <- 15
maxDP <- 100
dp_pass <- which(rowMeans(dp) >= minDP & rowMeans(dp) < maxDP);
length(dp_pass)
vcf_ld_bi_dp <- vcf_ld_bi[dp_pass, ]

# Check and remove monomorphic loci
vcf_ld_bi_dp_poly <- vcf_ld_bi_dp[which(is.polymorphic(vcf_ld_bi_dp)),]

# Plot loci read depths per sample
new_dp <- extract.gt(vcf_ld_bi_dp_poly, element = "DP", as.numeric = TRUE)
boxplot(new_dp, las = 3, col = c("#C0C0C0", "#808080"), ylab = "Depth")

# Quick visualisation of genetic structure using PCA
# https://lists.r-forge.r-project.org/pipermail/adegenet-forum/2018-February/001752.html
# library(adegenet); library(mapmixture); library(ggplot2)
# pca_geno <- vcfR2genlight(vcf_ld_bi_dp_poly)
# pca_pcal <- glPca(pca_geno, nf = 3)
# (pca_pcal$eig/sum(pca_pcal$eig)*100)[1:10]
# scatter_plot(as.data.frame(pca_pcal$scores), type = "labels", group_ids = indNames(pca_geno), size = 3)+theme(legend.position="none")

# Export filtered genotypes for Phymatolithon calcareum
vcfR::write.vcf(vcf_ld_bi_dp_poly, file = "./outputs/pcalcareum_SNPs.vcf.gz")

# Export filtered genotypes for Phymatolithon calcareum with sample size >= 3 (not AusII or Swa)
remove_less3 <- which(str_detect(colnames(vcf_ld_bi_dp_poly@gt), "AusII+|Swa", negate = TRUE))
remove_less3_vcf <- vcf_ld_bi_dp_poly[, remove_less3]
sort(colnames(remove_less3_vcf@gt))
remove_less3_vcf_poly <- remove_less3_vcf[which(is.polymorphic(remove_less3_vcf)),]
vcfR::write.vcf(remove_less3_vcf_poly, file = "./outputs/pcalcareum_SNPs_3greater_samples.vcf.gz")

# Export filtered genotypes for non-coarse Phymatolithon calcareum only
non_coarse <- which(str_detect(colnames(vcf_ld_bi_dp_poly@gt), "Maw11C+|Maw22C+|Maw22_03", negate = TRUE))
non_coarse_vcf <- vcf_ld_bi_dp_poly[, non_coarse]
colnames(non_coarse_vcf@gt)
non_coarse_vcf_poly <- non_coarse_vcf[which(is.polymorphic(non_coarse_vcf)),]
vcfR::write.vcf(non_coarse_vcf_poly, file = "./outputs/pcalcareum_SNPs_NonCoarse.vcf.gz")

# ----------------- #
# Lithothamnion corallioides ####
# ----------------- #

# Read in VCF file
vcf <- read.vcfR("./data/lcorallioides_variants_qc3_diploid.vcf")

# Filter SNPs for linkage
vcf_ld <- distance_thin(vcf, min.distance = 1000)

# Filter SNPs not biallelic
vcf_ld_bi <- filter_biallelic(vcf_ld)

# Extract total read depth per locus
dp <- extract.gt(vcf_ld_bi, element = "DP", as.numeric = TRUE)
queryMETA(vcf_ld_bi, element = "DP")

# Filter minimum and maximum mean depth (over all samples)
minDP <- 15
maxDP <- 100
dp_pass <- which(rowMeans(dp) >= minDP & rowMeans(dp) < maxDP);
length(dp_pass)
vcf_ld_bi_dp <- vcf_ld_bi[dp_pass, ]

# Check and remove monomorphic loci
vcf_ld_bi_dp_poly <- vcf_ld_bi_dp[which(is.polymorphic(vcf_ld_bi_dp)),]

# Plot loci read depths per sample
new_dp <- extract.gt(vcf_ld_bi_dp_poly, element = "DP", as.numeric = TRUE)
boxplot(new_dp, las = 3, col = c("#C0C0C0", "#808080"), ylab = "Depth")

# Quick visualisation of genetic structure using PCA
# https://lists.r-forge.r-project.org/pipermail/adegenet-forum/2018-February/001752.html
# library(adegenet); library(mapmixture); library(ggplot2)
# pca_geno <- vcfR2genlight(vcf_ld_bi_dp_poly)
# pca_pcal <- glPca(pca_geno, nf = 3)
# (pca_pcal$eig/sum(pca_pcal$eig)*100)[1:10]
# scatter_plot(as.data.frame(pca_pcal$scores), type = "labels", group_ids = indNames(pca_geno), size = 3)+theme(legend.position="none")

# Export filtered genotypes for Lithothamnion corallioides
vcfR::write.vcf(vcf_ld_bi_dp_poly, file = "./outputs/lcorallioides_SNPs.vcf.gz")

# Export filtered genotypes for Phymatolithon calcareum with sample size >= 3 (not Tud)
remove_less3 <- which(str_detect(colnames(vcf_ld_bi_dp_poly@gt), "Tud", negate = TRUE))
remove_less3_vcf <- vcf_ld_bi_dp_poly[, remove_less3]
sort(colnames(remove_less3_vcf@gt))
remove_less3_vcf_poly <- remove_less3_vcf[which(is.polymorphic(remove_less3_vcf)),]
vcfR::write.vcf(remove_less3_vcf_poly, file = "./outputs/lcorallioides_SNPs_3greater_samples.vcf.gz")
