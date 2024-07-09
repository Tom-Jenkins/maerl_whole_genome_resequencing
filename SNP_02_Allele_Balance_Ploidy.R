# =========================== #
#
# Maerl Whole Genome Resequencing Project 2024
#
# SNP Variant QC and Ploidy Analysis
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
library(pinfsc50)
library(adegenet)
library(stringr)


# ----------------- #
# Allele Balance Function
# ----------------- #

plot_allele_balance <- function(ad_matrix, sample, depth_threshold, col = "lightblue") {
  
  # Identify SNPs where the minor allele is greater than the depth threshold
  sample_ad <- ad_matrix[, sample]
  sample_ad_matrix <- apply(str_split(sample_ad, ",", simplify = T), 2, as.integer)
  snp_pass <- apply(sample_ad_matrix, 1, function(i) i[1] & i[2] > depth_threshold)
  
  # Subset SNPs that pass the threshold
  ad_pass <- sample_ad[snp_pass]
  ad_pass <- ad_pass[!is.na(ad_pass)]
  
  # Extract integer counts for allele 1 (highest count)
  allele1 <- unlist(lapply(strsplit(ad_pass, ","), function(i) max(as.integer(i))))
  
  # Extract integer counts for allele 2 (lowest count)
  allele2 <- unlist(lapply(strsplit(ad_pass, ","), function(i) min(as.integer(i))))

  # Calculate allele depth proportions from the counts
  ad1 <- allele1 / (allele1 + allele2)
  ad2 <- allele2 / (allele1 + allele2)
  
  # Plot allele balance proportions
  # hist(ad2, breaks = seq(0,1,by=0.02), col = "#1f78b4", xaxt="n", main = sample)
  # hist(ad1, breaks = seq(0,1,by=0.02), col = "#a6cee3", add = TRUE)
  # axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","1/3","3/4",1))
  
  # Plot allele balance proportions
  df <- data.frame(value = c(ad1,ad2))
  plt <- ggplot(df, aes(x = value)) +
    geom_histogram(binwidth = 0.01, fill = col, color = "black", linewidth = 0.05)+
    coord_cartesian(expand = FALSE, clip = "off")+
    scale_x_continuous(
      limits = c(0,1), breaks = c(0,1/4,1/3,1/2,2/3,1),
      labels = c("0","1/4","1/3","1/2","2/3","1")
    )+
    ylab("Frequency\n")+
    ggtitle(sample)+
    theme(
      panel.background = element_rect(fill = "white"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(colour = "grey90"),
      axis.title.x = element_blank(),
    )
  return(plt)
}


# ----------------- #
# Phymatolithon calcareum ####
# ----------------- #

# Read in VCF file
vcf <- read.vcfR("./outputs/pcalcareum_SNPs.vcf.gz")

# Extract genotypes
gt <- extract.gt(vcf, element = "GT")

# Extract total read depth per locus
dp_pcal <- extract.gt(vcf, element = "DP", as.numeric = TRUE)
sort(apply(dp_pcal, 2, median)) # Median read_pcal depth per sample

# Extract allele depth (number of read_pcals for each allele)
ad_pcal <- extract.gt(vcf, element = "AD")

# Subset heterozygotes by converting homozygotes (FALSE) to NA
is.na( ad_pcal[ !is_het(gt) ] ) <- TRUE

# Minimum allele depth threshold for interpreting allele balance
depth_threshold <- 10
sort(apply(dp_pcal, 2, median))

# Target sample
target_sample <- "Her_04"

# Plot allele balance
plot_allele_balance(ad_pcal, target_sample, depth_threshold)
dp_pcal[, target_sample] |> summary()

# # Subset all homozygous loci in 'godzilla' triploid maerl
# godzilla_idx <- c("FORMAT","Maw11C_01","Maw11C_02","Maw11C_04","Maw11C_06","Maw11C_07","Maw11C_09","Maw11C_10","Maw22_03","Maw22C_06","Maw22C_10","Maw22C_12","Maw22C_14","Maw22C_15","Maw22C_16","Maw22C_P11")
# godzilla_gt <- extract.gt(vcf[, godzilla_idx], element = "GT")
# godzilla_gt_het <- is_het(godzilla_gt)
# godzilla_gt_homo <- apply(godzilla_gt_het, 1, function(x) all(!godzilla_gt_het[x, ]))
# table(godzilla_gt_homo)
# vcf_homozygous <- vcf[which(godzilla_gt_homo), ]
# 
# # Check all loci are homozygous in godzilla maerl (only 0 or 2)
# geno <- vcfR2genlight(vcf_homozygous)
# tab(geno[which(indNames(geno) %in% godzilla_idx),]) |> table()
# 
# # All the above QC but for the triploid freebayes VCF file
# triploid <- read.vcfR("./data/pcalcareum_variants_triploid_qc2.vcf")
# triploid_ld <- distance_thin(triploid, min.distance = 1000)
# triploid_dp <- extract.gt(triploid_ld, element = "DP", as.numeric = TRUE)
# triploid_dp_pass <- which(rowMeans(triploid_dp) >= minDP & rowMeans(triploid_dp) < maxDP);
# length(triploid_dp_pass)
# triploid_ld_dp <- triploid_ld[triploid_dp_pass, ]


# ----------------- #
# Lithothamnion corallioides ####
# ----------------- #

# Read in VCF file
vcf <- read.vcfR("./outputs/lcorallioides_SNPs.vcf.gz")

# Extract genotypes
gt <- extract.gt(vcf, element = "GT")

# Extract total read_lcor depth per locus
dp_lcor <- extract.gt(vcf, element = "DP", as.numeric = TRUE)
sort(apply(dp_lcor, 2, median)) # Median read_lcor depth per sample

# Extract allele depth (number of read_lcors for each allele)
ad_lcor <- extract.gt(vcf, element = "AD")

# Subset heterozygotes by converting homozygotes (FALSE) to NA
is.na( ad_lcor[ !is_het(gt) ] ) <- TRUE

# Minimum allele depth threshold for interpreting allele balance
depth_threshold <- 10
sort(apply(dp_lcor, 2, median))

# Target sample
target_sample <- "Mil1_14"

# Plot allele balance
plot_allele_balance(ad_lcor, target_sample, depth_threshold)
dp_lcor[, target_sample] |> summary()


# ----------------- #
# Figure 3 ####
# ----------------- #

# Custom ggplot2 theme
custom_theme <- theme(
  panel.background = element_rect(fill = "white"),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.y = element_line(colour = "grey90"),
  axis.title.x = element_blank(),
)

# Parameters
n <- 100000
sd <- 0.08
col_diploid <- "#a6cee3"
col_triploid <- "#b2df8a"
set.seed(123)

# Diploid example 
diploid_df <- data.frame(value = rnorm(n, 1/2, sd))
diploid_plt <- ggplot(diploid_df, aes(x = value)) +
  geom_histogram(binwidth = 0.01, fill = col_diploid, color = "black", linewidth = 0.05)+
  coord_cartesian(expand = FALSE, clip = "off")+
  scale_x_continuous(
    limits = c(0,1), breaks = c(0,1/4,1/3,1/2,2/3,1),
    labels = c("0","1/4","1/3","1/2","2/3","1")
  )+
  ylab("Frequency\n")+
  ggtitle("Diploid")+
  custom_theme
diploid_plt

# Plot diploid samples
(Her_06 <- plot_allele_balance(ad_pcal, "Her_06", depth_threshold, col_diploid))
(Bor_11 <- plot_allele_balance(ad_pcal, "Bor_11", depth_threshold, col_diploid))
(Maw15_77 <- plot_allele_balance(ad_pcal, "Maw15_77", depth_threshold, col_diploid))
(Maw15_84 <- plot_allele_balance(ad_pcal, "Maw15_84", depth_threshold, col_diploid))

# Triploid simulated example
triploid_df <- data.frame(value = c(rnorm(n/2, 1/3, sd), rnorm(n/2, 2/3 , sd)))
triploid_plt <- ggplot(triploid_df, aes(x = value)) +
  geom_histogram(binwidth = 0.01, fill = col_triploid, color = "black", linewidth = 0.05)+
  coord_cartesian(expand = FALSE, clip = "off")+
  scale_x_continuous(
    limits = c(0,1), breaks = c(0,1/4,1/3,1/2,2/3,1),
    labels = c("0","1/4","1/3","1/2","2/3","1")
  )+
  ylab("Frequency\n")+
  ggtitle("Triploid")+
  custom_theme
triploid_plt

# Plot triploid samples
(Maw22C_16 <- plot_allele_balance(ad_pcal, "Maw22C_16", depth_threshold, col_triploid))
(Maw22C_06 <- plot_allele_balance(ad_pcal, "Maw22C_06", depth_threshold, col_triploid))
(Maw22C_14 <- plot_allele_balance(ad_pcal, "Maw22C_14", depth_threshold, col_triploid))
(Maw22C_12 <- plot_allele_balance(ad_pcal, "Maw22C_12", depth_threshold, col_triploid))

#--------------#
# Figure: Composer ####
#--------------#

# Load patchwork
library(patchwork)

# Layout design
layout <- "
  ABC
  DEF
"

# Plot layout
plt_list <- list(
  diploid_plt, Maw15_77, Maw15_84,
  triploid_plt, Maw22C_16+ggtitle("Maw22C_16 (Coarse)"), Maw22C_06+ggtitle("Maw22C_16 (Coarse)")
)
wrap_plots(plt_list, design = layout)
# + plot_annotation(tag_levels = "A")
ggsave("Figure_03.png", width = 12, height = 8, units = "in", dpi = 600)
ggsave("Figure_03.pdf", width = 12, height = 8, units = "in")
