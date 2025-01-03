# =========================== #
#
# Maerl Whole Genome Re-sequencing Project 2024
#
# SNP Clonality and Genetic Diversity
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
# Phymatolithon calcareum: clones ####
# ----------------- #

# Read in Pcalcareum SNPs
pcal_genind <- vcfR2genind(read.vcfR("outputs/pcalcareum_SNPs.vcf.gz"))
pcal_genind

# Merge the following samples into the same group (population) for Fst analysis
# • St Mawes coarse maerl from 2011 and 2022
# • St Mawes maerl from 2015 and 2022
# • All samples from St Austell area (St Austell Bay and Little Gribbin)
pcal_genind$pop <- as.factor(sub("\\_.*", "", indNames(pcal_genind)))
pcal_genind$pop <- as.factor(str_replace(pcal_genind$pop, "Maw22C|Maw11C", "MawC"))
pcal_genind$pop <- as.factor(str_replace(pcal_genind$pop, "Maw15|Maw22", "Maw"))
pcal_genind$pop <- as.factor(str_replace(pcal_genind$pop, "Gri|AusII", "Aus"))
pcal_genind$pop <- as.factor(str_replace(pcal_genind$pop, "Her", "Nar"))
summary(pcal_genind$pop)

# Calculate Prevosti's genetic distance
pcal_dist = prevosti.dist(pcal_genind)

# Isolate distance value for two individuals
pcal_dist_mat <- as.matrix(pcal_dist)
ind1 <- "Maw22C_14"
ind2 <- "Maw22C_15"
pcal_dist_mat[ind1, ind2]

# Clone cut-off threshold based on Prevosti's genetic distance
pcal_threshold <- mlg.filter(pcal_genind, distance = pcal_dist, stats = "THRESHOLDS", threshold = 1)
(pcal_clone_threshold <- cutoff_predictor(pcal_threshold))

# Plot histogram of distance values
png("supporting_figures/Figure_SX_PCAL_HIST.png", width = 10, height = 5, units = "in", res = 300)
par(mfrow = c(1, 2))
hist(pcal_dist, breaks = 100, xlab = "Prevosti's genetic distance")
abline(v = pcal_clone_threshold, col = "red")
hist(pcal_dist[pcal_dist < 0.12], breaks = 10, xlab = "Prevosti's genetic distance")
abline(v = pcal_clone_threshold, col = "red")
dev.off()

# Create a table of Multi-Locus Genotypes (MLGs) using clone cut-off threshold
pcal_clones_idx <- mlg.filter(pcal_genind, threshold = pcal_clone_threshold, distance = pcal_dist, algorithm = "nearest_neighbor")
pcal_clones_df <- data.frame(
  IND = indNames(pcal_genind),
  POP = pcal_genind$pop,
  MLG = factor(pcal_clones_idx, labels = 1:n_distinct(pcal_clones_idx))
) |> arrange(.data = _, MLG, POP)
pcal_clones_df
pcal_clones_df |> arrange(.data = _, POP, IND, MLG)

# Number of unique genotypes
# (equal to the number of individuals due to polymorphic SNP filters)
poppr::mlg(pcal_genind)

# Number of MLGs
(nMLGs <- pcal_clones_df$MLG |> table() |> length())

# Create a MLG dataset with clones removed (one individual randomly selected from each MLG and POP)
pcal_clones_df$POP_MLG <- str_c(pcal_clones_df$POP, "_", pcal_clones_df$MLG)
pcal_clones_tibble <- pcal_clones_df |> group_by(POP_MLG) |> summarise(IND_list = list(IND))
inds_no_clones <- unlist(lapply(1:nrow(pcal_clones_tibble), function(i) sample(pcal_clones_tibble[i,]$IND_list[[1]], size = 1)))
pcal_MLG <- pcal_genind[inds_no_clones, ]
pcal_MLG <- pcal_MLG[loc = names(which(isPoly(pcal_MLG)))] # keep only polymorphic SNPs
pcal_MLG
indNames(pcal_MLG)

# Facet labels
facet_labs <- paste0(names(table(pcal_clones_df$POP)), paste0(" (N=", table(pcal_clones_df$POP), ")"))
pcal_clones_df$POP_LABS <- factor(pcal_clones_df$POP, labels = facet_labs)

# Add column that represents whether a clonal lineage is spread over more than one site
# Lineages spread: 19, 22
pcal_clones_df$SPREAD <- ifelse(pcal_clones_df$MLG == 19, "19",
                           ifelse(pcal_clones_df$MLG == 22, "22",
                                  "Site-Specific Lineage"))

# Filter dataset
pcal_clones_df <- pcal_clones_df |>
  filter(POP != "Bor", POP != "Mor", POP != "Ons", POP != "Tre", POP != "Swa")

# Plot MLGs
mlg_pcal <- ggplot(data = pcal_clones_df)+
  geom_bar(aes(x = MLG, fill = SPREAD), show.legend = TRUE, colour = NA, linewidth = 0.1)+
  scale_fill_manual(values = c("#cbd5e8","#fdcdac","#cccccc"))+
  facet_wrap(~ POP_LABS, scales = "free", nrow = 3)+
  scale_y_continuous(
    limits = function(x) c(0, max(x)),
    breaks = function(x) if(max(x) <= 5) {seq(0, max(x), 1)} else {seq(0, max(x), 5)}
  )+
  ylab("Number of Individuals")+
  xlab("Clonal lineage (MLG)")+
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90"),
    strip.background = element_rect(fill = "grey30"),
    strip.text = element_text(colour = "white", face = "bold", size = 10),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
  )
mlg_pcal
  
# Prepare a MLG table for running Pareto Beta
pcal_mlg_table <- count(pcal_clones_df, MLG, POP) |> 
  mutate(MLG = str_c("MLG.", MLG)) |> 
  pivot_wider(names_from = MLG, values_from = n) |> 
  column_to_rownames(var = "POP") |> 
  as.matrix()
pcal_mlg_table[is.na(pcal_mlg_table)] <- 0
pcal_mlg_table

# Copy functions from ?poppr::diversity_stats
library("poweRlaw")
power_law_beta <- function(x){
  xpow <- displ(x[x > 0])                 # Generate the distribution
  xpow$setPars(estimate_pars(xpow))       # Estimate the parameters
  xdat <- plot(xpow, draw = FALSE)        # Extract the data
  xlm <- lm(log(y) ~ log(x), data = xdat) # Run log-log linear model for slope
  return(-coef(xlm)[2])
}
Beta <- function(x){
  x <- drop(as.matrix(x))
  if (length(dim(x)) > 1){
    res <- apply(x, 1, power_law_beta)
  } else {
    res <- power_law_beta(x)
  }
  return(res)
}

# Calculate Pareto Beta using Poppr function (first remove sites with one MLG)
diversity_stats(pcal_mlg_table) # NaN for sites where only one MLG is present
mlg_table_Beta <- pcal_mlg_table[which(!is.nan(diversity_stats(pcal_mlg_table)[,4])), ]
diversity_stats(mlg_table_Beta, B = Beta) |> as.data.frame() |> rownames_to_column() |> arrange(rowname)

# Compute heterozygosity and Fis
pcal_basic_stats <- basic.stats(pcal_genind)
pcal_basic_stats$Ho |> apply(X = _, MARGIN = 2, FUN = mean, na.rm=TRUE) |> round(2)
pcal_basic_stats$Hs |> apply(X = _, MARGIN = 2, FUN = mean, na.rm=TRUE) |> round(2)
pcal_basic_stats$Fis |> apply(X = _, MARGIN = 2, FUN = mean, na.rm=TRUE) |> round(2)

# Calculate number of somatic mutations for coarse St Mawes
coarse_geno <- popsub(pcal_genind, sublist = "MawC") |> genind2df()

somatic_mutations <- function(column) {
  # Find the most common genotype (mode) in the column
  most_common <- names(sort(table(column), decreasing = TRUE))[1]
  
  # Count individuals that have a different genotype to most common
  sum(column != most_common)
}

# Number of somatic mutations across loci
# 1 = one individual has somatic mutation at locus
# 5 = five individuals have somatic mutation at locus
num_mutations <- sapply(coarse_geno[-1], somatic_mutations)
table(num_mutations)

# Percentage somatic mutations
sum(table(num_mutations)[2:length(table(num_mutations))]) / sum(table(num_mutations))


# ----------------- #
# Lithothamnion corallioides: clones ####
# ----------------- #

# Read in Lcorallioides SNPs
lcor_genind <- vcfR2genind(read.vcfR("./outputs/lcorallioides_SNPs.vcf.gz"))

# Merge the following samples into the same group (population)
# • St Mawes maerl from 2015 and 2022
# • All samples from St Austell area (St Austell Bay and Little Gribbin)
lcor_genind$pop <- as.factor(sub("\\_.*", "", indNames(lcor_genind)))
lcor_genind$pop <- as.factor(str_replace(lcor_genind$pop, "Maw15|Maw22", "Maw"))
lcor_genind$pop <- as.factor(str_replace(lcor_genind$pop, "Aus[A-Z]+|Gri", "Aus"))
lcor_genind$pop <- as.factor(str_replace(lcor_genind$pop, "Wey|Swa", "WeySwa"))
lcor_genind$pop <- as.factor(str_replace(lcor_genind$pop, "Mil[1-2]*", "Mil"))
summary(lcor_genind$pop)

# Calculate Prevosti's genetic distance
lcor_dist = prevosti.dist(lcor_genind)

# Isolate distance value for two individuals
lcor_dist_mat <- as.matrix(lcor_dist)
ind1 <- "Mil1_11"
ind2 <- "Mil1_13"
lcor_dist_mat[ind1, ind2]

# Clone cut-off threshold based on Prevosti's genetic distance
lcor_threshold <- mlg.filter(lcor_genind, distance = lcor_dist, stats = "THRESHOLDS", threshold = 1)
(lcor_clone_threshold <- cutoff_predictor(lcor_threshold))

# Plot histogram of distance values
png("supporting_figures/Figure_SX_LCOR_HIST.png", width = 10, height = 5, units = "in", res = 300)
par(mfrow = c(1, 2))
hist(lcor_dist, breaks = 100, xlab = "Prevosti's genetic distance")
abline(v = lcor_clone_threshold, col = "red")
hist(lcor_dist[lcor_dist < 0.25], breaks = 20, xlab = "Prevosti's genetic distance")
abline(v = lcor_clone_threshold, col = "red")
dev.off()

# Create a table of Multi-Locus Genotypes (MLGs) using clone cut-off threshold
lcor_clones_idx <- mlg.filter(lcor_genind, threshold = lcor_clone_threshold, distance = lcor_dist, algorithm = "nearest_neighbor")
lcor_clones_df <- data.frame(
  IND = indNames(lcor_genind),
  POP = lcor_genind$pop,
  MLG = factor(lcor_clones_idx, labels = 1:n_distinct(lcor_clones_idx))
)
lcor_clones_df |> arrange(.data = _, POP, IND, MLG)
lcor_clones_df |> arrange(.data = _, MLG, POP)

# Number of unique genotypes (equal to the number of individuals due to polymorphic SNP filters)
poppr::mlg(lcor_genind)

# Number of MLGs
(nMLGs <- lcor_clones_df$MLG |> table() |> length())

# Create a MLG dataset with clones removed (one individual selected from each MLG and POP)
lcor_clones_df$POP_MLG <- str_c(lcor_clones_df$POP, "_", lcor_clones_df$MLG)
lcor_clones_tibble <- lcor_clones_df |> group_by(POP_MLG) |> summarise(IND_list = list(IND))
lcor_inds_no_clones <- unlist(lapply(1:nrow(lcor_clones_tibble), function(i) sample(lcor_clones_tibble[i,]$IND_list[[1]], size = 1)))
lcor_MLG <- lcor_genind[lcor_inds_no_clones, ]
lcor_MLG <- lcor_MLG[loc = names(which(isPoly(lcor_MLG)))] # keep only polymorphic SNPs
lcor_MLG

# Facet labels
facet_labs <- paste0(names(table(lcor_clones_df$POP)), paste0(" (N=", table(lcor_clones_df$POP), ")"))
lcor_clones_df$POP_LABS <- factor(lcor_clones_df$POP, labels = facet_labs)

# Add column that represents whether a clonal lineage is spread over more than one site
# Lineages spread: None
# clones_df$SPREAD <- ifelse(clones_df$MLG == 38, "38", ifelse(clones_df$MLG == 44, "44", "n/a"))

# Filter dataset
lcor_clones_df <- lcor_clones_df |>
  filter(POP != "Tud")

# Plot MLGs
mlg_lcor <- ggplot(data = lcor_clones_df)+
  geom_bar(aes(x = MLG, fill = NULL), fill = "#cccccc", show.legend = TRUE, colour = NA, linewidth = 0.1)+
  # scale_fill_manual("Multi-Site Lineages", values = c("#f1e2cc","#b3e2cd","#cccccc"))+
  facet_wrap(~ POP_LABS, scales = "free", nrow = 2)+
  scale_y_continuous(
    limits = function(x) c(0, max(x)),
    breaks = function(x) if(max(x) <= 5) {seq(0, max(x), 1)} else {seq(0, max(x), 5)}
  )+
  ylab("Number of Individuals")+
  xlab("\nClonal lineage (MLG)")+
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90"),
    strip.background = element_rect(fill = "grey30"),
    strip.text = element_text(colour = "white", face = "bold", size = 10),
    legend.position = "top",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
  )
mlg_lcor

# Prepare a MLG table for running Pareto Beta
lcor_mlg_table <- count(lcor_clones_df, MLG, POP) |> 
  mutate(MLG = str_c("MLG.", MLG)) |> 
  pivot_wider(names_from = MLG, values_from = n) |> 
  column_to_rownames(var = "POP") |> 
  as.matrix()
lcor_mlg_table[is.na(lcor_mlg_table)] <- 0
lcor_mlg_table

# Copy functions from ?poppr::diversity_stats
library("poweRlaw")
power_law_beta <- function(x){
  xpow <- displ(x[x > 0])                 # Generate the distribution
  xpow$setPars(estimate_pars(xpow))       # Estimate the parameters
  xdat <- plot(xpow, draw = FALSE)        # Extract the data
  xlm <- lm(log(y) ~ log(x), data = xdat) # Run log-log linear model for slope
  return(-coef(xlm)[2])
}
Beta <- function(x){
  x <- drop(as.matrix(x))
  if (length(dim(x)) > 1){
    res <- apply(x, 1, power_law_beta)
  } else {
    res <- power_law_beta(x)
  }
  return(res)
}

# Calculate Pareto Beta using Poppr function (first remove sites with one MLG)
diversity_stats(lcor_mlg_table) # NaN for sites where only one MLG is present
mlg_table_Beta <- lcor_mlg_table[which(!is.nan(diversity_stats(lcor_mlg_table)[,4])), ]
diversity_stats(mlg_table_Beta, B = Beta) |> as.data.frame() |> rownames_to_column() |> arrange(rowname)

# Compute heterozygosity and Fis
lcor_basic_stats <- basic.stats(lcor_genind)
lcor_basic_stats$Ho |> apply(X = _, MARGIN = 2, FUN = mean, na.rm=TRUE) |> round(2)
lcor_basic_stats$Hs |> apply(X = _, MARGIN = 2, FUN = mean, na.rm=TRUE) |> round(2)
lcor_basic_stats$Fis |> apply(X = _, MARGIN = 2, FUN = mean, na.rm=TRUE) |> round(2)


#--------------#
# Figure 4: Composer ####
#--------------#

# Load patchwork
library(patchwork)

# Layout design
layout <- "
  A
  A
  A
  B
  B
"

# Plot layout
plt_list <- list(
  mlg_pcal + theme(axis.title.x = element_blank(), plot.title = element_text(size = 14, face = "bold")) + ggtitle(expression(italic("Phymatolithon calcareum")))+ theme(plot.margin = margin(b = 20)),
  mlg_lcor + theme(axis.title.x = element_text(size = 14), plot.title = element_text(size = 14, face = "bold")) + ggtitle(expression(italic("Lithothamnion corallioides")))
)
wrap_plots(plt_list, design = layout)
ggsave("figures/Figure_05.png", width = 10, height = 10, units = "in", dpi = 600)
ggsave("figures/Figure_05.pdf", width = 10, height = 10, units = "in")


