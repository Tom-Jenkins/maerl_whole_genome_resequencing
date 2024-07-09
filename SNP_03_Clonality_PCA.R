# =========================== #
#
# Maerl Whole Genome Resequencing Project 2024
#
# SNP Genetic Structure Analysis
#
# Species:
# Phymatolithon calcareum
# Lithothamnion corallioides
#
# SNP data files:
# ./outputs/
# ./outputs/
#
# =========================== #

# In RStudio set working directory to the path where this R script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load packages
library(vcfR)
library(adegenet)
library(stringr)
library(dplyr)
library(mapmixture)
library(poppr)
library(hierfstat)
library(mmod)
library(pegas)
library(dartR)
library(RColorBrewer)
library(scales)
library(patchwork)

# Read in Pcalcareum SNPs
pcal_genind <- vcfR2genind(read.vcfR("./outputs/pcalcareum_SNPs.vcf.gz"))
indNames(pcal_genind) |> sort()

# Add and edit population labels
pcal_genind$pop <- as.factor(sub("\\_.*", "", indNames(pcal_genind)))
pcal_genind$pop <- as.factor(str_replace(pcal_genind$pop, "Maw22C|Maw11C", "MawC"))
pcal_genind$pop <- as.factor(str_replace(pcal_genind$pop, "Maw15|Maw22", "Maw"))
summary(pcal_genind$pop)

# Read in Lcorallioides SNPs
lcor_genind <- vcfR2genind(read.vcfR("./outputs/lcorallioides_SNPs.vcf.gz"))

# Add and edit population labels
lcor_genind$pop <- as.factor(sub("\\_.*", "", indNames(lcor_genind)))
lcor_genind$pop <- as.factor(str_replace(lcor_genind$pop, "Maw15|Maw22", "Maw"))
summary(lcor_genind$pop)

# ----------------- #
# Phymatolithon calcareum: clones ####
# ----------------- #

# Calculate Prevosti's genetic distance
ind_dist = prevosti.dist(pcal_genind)

# Plot histogram of distance values
hist(ind_dist, breaks = 100)
hist(ind_dist[ind_dist < 0.10], breaks = 10)

# Isolate distance value for two individuals
ind_dist_mat <- as.matrix(ind_dist)
ind1 <- "Maw15_P1"
ind2 <- "Maw22_04"
ind_dist_mat[ind1, ind2]

# Clone cut-off threshold based on Prevosti's genetic distance
thresholds <- mlg.filter(pcal_genind, distance = ind_dist, stats = "THRESHOLDS", threshold = 1)
(clone_threshold <- cutoff_predictor(thresholds))

# Create a table of Multi-Locus Genotypes (MLGs) using clone cut-off threshold
clones_idx <- mlg.filter(pcal_genind, threshold = clone_threshold, distance = ind_dist, algorithm = "nearest_neighbor")
clones_df <- data.frame(
  IND = indNames(pcal_genind),
  POP = pcal_genind$pop,
  MLG = factor(clones_idx, labels = 1:n_distinct(clones_idx))
)
clones_df |> arrange(.data = _, POP, IND, MLG)
clones_df |> arrange(.data = _, MLG, POP)

# Number of unique genotypes (equal to the number of individuals due to polymorphic SNP filters)
poppr::mlg(pcal_genind)

# Number of MLGs
(nMLGs <- clones_df$MLG |> table() |> length())

# Facet labels
facet_labs <- paste0(names(table(clones_df$POP)), paste0(" (N=", table(clones_df$POP), ")"))
clones_df$POP_LABS <- factor(clones_df$POP, labels = facet_labs)

# Add column that represents whether a clonal lineage is spread over more than one site
# Lineages spread: 38, 41
clones_df$SPREAD <- ifelse(clones_df$MLG == 38, "38", ifelse(clones_df$MLG == 41, "41", "n/a"))

# Plot MLGs
mlg_pcal <- ggplot(data = clones_df)+
  geom_bar(aes(x = MLG, fill = SPREAD), show.legend = TRUE, colour = "black", linewidth = 0.1)+
  scale_fill_manual("Multi-Site Lineages", values = c("#cbd5e8","#fdcdac","#cccccc"))+
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
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
  )
mlg_pcal
  
# Create a MLG dataset with clones removed (one individual randomly selected from each MLG and POP)
clones_df$POP_MLG <- str_c(clones_df$POP, "_", clones_df$MLG)
clones_tibble <- clones_df |> group_by(POP_MLG) |> summarise(IND_list = list(IND))
set.seed(123)
inds_no_clones <- unlist(lapply(1:nrow(clones_tibble), function(i) sample(clones_tibble[i,]$IND_list[[1]], size = 1)))
pcal_MLG <- pcal_genind[inds_no_clones, ]
pcal_MLG <- pcal_MLG[loc = names(which(isPoly(pcal_MLG)))] # keep only polymorphic SNPs
pcal_MLG
indNames(pcal_MLG)


# ----------------- #
# Lithothamnion corallioides: clones ####
# ----------------- #

# Calculate Prevosti's genetic distance
ind_dist = prevosti.dist(lcor_genind)

# Plot histogram of distance values
hist(ind_dist, breaks = 100)
hist(ind_dist[ind_dist < 0.25], breaks = 20)

# Isolate distance value for two individuals
ind_dist_mat <- as.matrix(ind_dist)
ind1 <- "Mil1_11"
ind2 <- "Mil1_13"
ind_dist_mat[ind1, ind2]

# Clone cut-off threshold based on Prevosti's genetic distance
thresholds <- mlg.filter(lcor_genind, distance = ind_dist, stats = "THRESHOLDS", threshold = 1)
(clone_threshold <- cutoff_predictor(thresholds))

# Create a table of Multi-Locus Genotypes (MLGs) using clone cut-off threshold
clones_idx <- mlg.filter(lcor_genind, threshold = clone_threshold, distance = ind_dist, algorithm = "nearest_neighbor")
clones_df <- data.frame(
  IND = indNames(lcor_genind),
  POP = lcor_genind$pop,
  MLG = factor(clones_idx, labels = 1:n_distinct(clones_idx))
)
clones_df |> arrange(.data = _, POP, IND, MLG)
clones_df |> arrange(.data = _, MLG, POP)

# Number of unique genotypes (equal to the number of individuals due to polymorphic SNP filters)
poppr::mlg(lcor_genind)

# Number of MLGs
(nMLGs <- clones_df$MLG |> table() |> length())

# Facet labels
facet_labs <- paste0(names(table(clones_df$POP)), paste0(" (N=", table(clones_df$POP), ")"))
clones_df$POP_LABS <- factor(clones_df$POP, labels = facet_labs)

# Add column that represents whether a clonal lineage is spread over more than one site
# Lineages spread: 38, 44
clones_df$SPREAD <- ifelse(clones_df$MLG == 38, "38", ifelse(clones_df$MLG == 44, "44", "n/a"))

# Plot MLGs
mlg_lcor <- ggplot(data = clones_df)+
  geom_bar(aes(x = MLG, fill = SPREAD), show.legend = TRUE, colour = "black", linewidth = 0.1)+
  scale_fill_manual("Multi-Site Lineages", values = c("#f1e2cc","#b3e2cd","#cccccc"))+
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

# Create a MLG dataset with clones removed (one individual selected from each MLG and POP)
clones_df$POP_MLG <- str_c(clones_df$POP, "_", clones_df$MLG)
clones_tibble <- clones_df |> group_by(POP_MLG) |> summarise(IND_list = list(IND))
set.seed(1234)
inds_no_clones <- unlist(lapply(1:nrow(clones_tibble), function(i) sample(clones_tibble[i,]$IND_list[[1]], size = 1)))
lcor_MLG <- lcor_genind[inds_no_clones, ]
lcor_MLG <- lcor_MLG[loc = names(which(isPoly(lcor_MLG)))] # keep only polymorphic SNPs
lcor_MLG


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
  mlg_pcal + theme(axis.title.x = element_blank(), plot.title = element_text(size = 14, face = "bold")) + ggtitle(expression(italic("Phymatolithon calcareum"))),
  mlg_lcor + theme(axis.title.x = element_text(size = 14), plot.title = element_text(size = 14, face = "bold")) + ggtitle(expression(italic("Lithothamnion corallioides")))
)
wrap_plots(plt_list, design = layout)
ggsave("Figure_04.png", width = 10, height = 10, units = "in", dpi = 600)
ggsave("Figure_04.pdf", width = 10, height = 10, units = "in")


# ----------------- #
# Phymatolithon calcareum: PCA ####
# ----------------- #

# Subset data without coarse maerl
pops <- c("AusII","Biz","Bor","Ger","Gri","Her","Man","Maw","Maw","Mor","Ons","Swa","Tre","Wey","Zar")
pcal_genind_sub <- pcal_genind[pop = pops]
pcal_genind_sub <- pcal_genind_sub[loc = names(which(isPoly(pcal_genind_sub)))]

# Colour palette
library(randomcoloR)
set.seed(123)
palette <- distinctColorPalette(nPop(pcal_genind))

# Principal component analysis: all individuals axis 1 vs. 2
pca_pcal1 <- dudi.pca(pcal_genind, scannf = FALSE, nf = 5)
(pca_pcal1_plt <- scatter_plot(
  dataframe = pca_pcal1$li, group_ids = pcal_genind$pop, labels = pcal_genind$pop,
  type = "labels", axes = c(1,2), size = 3, colours = palette,
  percent = round(pca_pcal1$eig/sum(pca_pcal1$eig)*100, 1)
))

# Principal component analysis: all individuals axis 2 vs. 3
pca_pcal2 <- dudi.pca(pcal_genind, scannf = FALSE, nf = 5)
(pca_pcal2_plt <- scatter_plot(
  dataframe = pca_pcal1$li, group_ids = pcal_genind$pop, labels = pcal_genind$pop,
  type = "labels", axes = c(2,3), size = 3, colours = palette,
  percent = round(pca_pcal1$eig/sum(pca_pcal1$eig)*100, 1)
))

# Combine both plots
pcal_pcas <- wrap_plots(pca_pcal1_plt, pca_pcal2_plt, guides = "collect")

# Principal component analysis: MLGs
pca_pcal3 <- dudi.pca(popsub(pcal_MLG, exclude = c("MawC", "Wey")), scannf = FALSE, nf = 5)
pca_pcal3_plt <- scatter_plot(pca_pcal3$li, group_ids = popsub(pcal_MLG, exclude = c("MawC", "Wey"))$pop, type = "labels", axes = c(1,2), size = 3, colours = palette)
pca_pcal4_plt <- scatter_plot(pca_pcal3$li, group_ids = popsub(pcal_MLG, exclude = c("MawC", "Wey"))$pop, type = "labels", axes = c(2,3), size = 3, colours = palette)

# Combine both plots
pcal_pcas_mlgs <- wrap_plots(pca_pcal3_plt, pca_pcal4_plt, guides = "collect")
pcal_pcas_mlgs
ggsave("supporting_figures/Figure_S4.png", width = 15, height = 7, units = "in", dpi = 600)

# Principal component analysis: without coarse maerl
# pca_pcal2 <- dudi.pca(pcal_genind_sub, scannf = FALSE, nf = 5)
# scatter_plot(pca_pcal2$li, group_ids = pcal_genind_sub$pop, type = "labels")

# ----------------- #
# Lithothamnion corallioides: PCA ####
# ----------------- #

# Principal component analysis: all individuals
lcor_pca1 <- dudi.pca(lcor_genind, scannf = FALSE, nf = 5)
pca3$li[which(pca3$li$Axis1 > 50), ]
# CHECK Hel_57 and Hel_38 balance
# Interpretation:
# Clone could be very old

# Principal component analysis: all individuals axis 1 vs. 2
(pca_lcor1_plt <- scatter_plot(
  dataframe = lcor_pca1$li, group_ids = lcor_genind$pop, labels = lcor_genind$pop,
  type = "labels", axes = c(1,2), size = 3, colours = palette,
  percent = round(lcor_pca1$eig/sum(lcor_pca1$eig)*100, 1)
))

# Principal component analysis: all individuals axis 2 vs. 3
(pca_lcor2_plt <- scatter_plot(
  dataframe = lcor_pca1$li, group_ids = lcor_genind$pop, labels = lcor_genind$pop,
  type = "labels", axes = c(2,3), size = 3, colours = palette,
  percent = round(lcor_pca1$eig/sum(lcor_pca1$eig)*100, 1)
))

wrap_plots(pca_lcor1_plt, pca_lcor2_plt, guides = "collect")
ggsave("supporting_figures/Figure_S6.png", width = 15, height = 7, units = "in", dpi = 600)

# Principal component analysis: MLGs
lcor_pca2 <- dudi.pca(lcor_MLG, scannf = FALSE, nf = 5)

# Principal component analysis: all individuals axis 1 vs. 2
(pca_lcor3_plt <- scatter_plot(
  dataframe = lcor_pca2$li, group_ids = lcor_MLG$pop, labels = lcor_MLG$pop,
  type = "labels", axes = c(1,2), size = 3, colours = palette,
  percent = round(lcor_pca2$eig/sum(lcor_pca2$eig)*100, 1)
))
ggsave("supporting_figures/Figure_Lcor_MLGs_only.png", width = 10, height = 7, units = "in", dpi = 600)
