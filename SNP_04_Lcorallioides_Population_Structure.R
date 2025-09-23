# ========== #
#
# Maerl WGS Analysis 2025
#
# Clonality, Genetic Diversity and Population Structure
#
# Species:
# Lithothamnion corallioides
#
# ========== #

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
library(dartR)
library(LEA)
library(sf)
library(ggplot2)
library(patchwork)
library(rnaturalearth)
library(rnaturalearthhires)
# library(shadowtext)

# ---------- #
# Lithothamion corallioides clonal diversity ####
# ---------- #

# Read in Lcorallioides SNPs
lcor_genind <- vcfR2genind(read.vcfR("./outputs/Lcorallioides_SNPs.vcf.gz"))

# Merge the following samples into the same group (population)
# • St Mawes maerl from 2015 and 2022
# • All samples from St Austell area (St Austell Bay and Little Gribbin)
lcor_genind$pop <- as.factor(sub("\\_.*", "", indNames(lcor_genind)))
lcor_genind$pop <- as.factor(str_replace(lcor_genind$pop, "Maw15|Maw22", "Maw"))
lcor_genind$pop <- as.factor(str_replace(lcor_genind$pop, "Aus[A-Z]+|Gri", "Aus"))
# lcor_genind$pop <- as.factor(str_replace(lcor_genind$pop, "Wey|Swa", "WeySwa"))
# lcor_genind <- popsub(lcor_genind, exclude = c("Tud"))
# lcor_genind$pop <- as.factor(str_replace(lcor_genind$pop, "Mil[1-2]*", "Mil"))
summary(lcor_genind$pop)

# Check and remove any monomorphic SNPS after filtering individuals and sites
lcor_genind <- lcor_genind[loc=names(which(isPoly(lcor_genind)))] 

# Convert genind to snpclone
lcor_snpclone <- as.snpclone(gi2gl(lcor_genind))

# Number of multi-locus genotypes (MLGs)
# Different from the multi-locus lineages (MLL) / clonal lineage / genet
mlg(lcor_genind)

# See Arnaud-Haond et al. 2007
# https://onlinelibrary.wiley.com/doi/10.1111/j.1365-294X.2007.03535.x#ss3-title

# Calculate Hamming genetic distance
lcor_hamming = bitwise.dist(lcor_snpclone, mat = FALSE, percent = TRUE)
hist(lcor_hamming, breaks = 100, xlab = "Hamming genetic distance", main = "Hamming genetic distance")

# Clone cut-off threshold for defining MLL
max_distance <- max(as.vector(lcor_hamming))
lcor_threshold <- mlg.filter(lcor_snpclone, distance = lcor_hamming, stats = "THRESHOLDS", threshold = max_distance, algorithm = "average_neighbor")
(lcor_clone_threshold <- cutoff_predictor(lcor_threshold))

# Plot histogram of Hamming distance values
png("supporting_figures/Lcorallioides_hamming_histogram.png", width = 10, height = 5, units = "in", res = 300)
hist(lcor_hamming, breaks = 100, xlab = "Hamming genetic distance", main = "Hamming genetic distance")
abline(v = lcor_clone_threshold, col = "red")
dev.off()

# Filter snpclone object to get the multi-locus lineages (MLLs)
# Proxy of the number of genets
mlg.filter(lcor_snpclone, distance = lcor_hamming, algorithm = "average_neighbor") <- lcor_clone_threshold
lcor_snpclone
mlg(lcor_snpclone)

# Table of clonal lineages
lcor_snpclone_tab <- mlg.table(lcor_snpclone)

# Prepare data.frame in long format, filter out zeros and order by MLL number
lcor_snpclone_tab_long <- lcor_snpclone_tab |> 
  as.data.frame() |> 
  cbind(Site = rownames(lcor_snpclone_tab)) |> 
  pivot_longer(
    cols = starts_with("MLG"), 
    names_to = "MLG", 
    values_to = "Count",
  ) |> 
  dplyr::filter(Count != 0) |> 
  mutate(MLG = str_replace(MLG, "MLG.", "")) |> 
  mutate(MLG = factor(MLG, levels = sort(unique(as.integer(MLG)))))

# Check whether MLLs cross populations
mlg.crosspop(lcor_snpclone)

# Add column that represents whether a clonal lineage overlaps with more than one site
(overlapping_mll <- str_remove(names(mlg.crosspop(lcor_snpclone)), "MLG."))
lcor_snpclone_tab_long <- lcor_snpclone_tab_long |>
  mutate(
    Overlap = case_when(
      MLG == as.integer(overlapping_mll[1]) ~ overlapping_mll[1],
      MLG == as.integer(overlapping_mll[2]) ~ overlapping_mll[2],
      TRUE      ~ "Site-Specific Lineage"
    )
  )

# Plot MLLs
(lcor_MLLs <- ggplot(data = lcor_snpclone_tab_long)+
  geom_bar(aes(x=MLG, y=Count, fill=Overlap), stat = "identity", linewidth = 0.1)+
  scale_fill_manual(values = c("#90a4ae","#a1887f","grey70"))+
  facet_wrap(~Site, scales = "free", ncol = 4)+
  scale_y_continuous(
    limits = c(0, max(lcor_snpclone_tab_long$Count)),
    breaks = c(0, 1, 2, 3, 7)
  )+
  ylab("Number of clones")+
  xlab("Clonal lineage (MLL)")+
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
    axis.title.x = element_text(vjust = -1),
  )
)

# Subset MawC MLL-19
hel <- popsub(lcor_snpclone, sublist = "Hel")
mlg46 <- hel[hel@mlg == 47, ]
mlg46_different_alleles <- bitwise.dist(mlg46, mat = FALSE, percent = FALSE)
mean(mlg46_different_alleles); median(mlg46_different_alleles); range(mlg46_different_alleles)

# Total differences per individual
total_diff <- rowSums(as.matrix(mlg46_different_alleles))

# Find the "consensus" individual, i.e. genetically smallest distance to all others
consensus_idx <- which.min(total_diff)

# Extract differences from consensus clone
# I.e. candidate somatic mutations
diff_to_consensus <- as.matrix(mlg46_different_alleles)[, consensus_idx]
diff_to_consensus

# In MLG 19, individuals differ from the consensus by X–Y alleles (putative somatic mutations).
# If we knew the age of these maerl samples, we might be able to estimate somatic mutation rate.
paste("Mean:", mean(diff_to_consensus[-1]))
paste("Median:", median(diff_to_consensus[-1]))
paste("Range:", range(diff_to_consensus[-1]))

# Function to calculate Pareto Beta ?poppr::diversity_stats
library(poweRlaw)
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

# Clonal diversity statistics (remove Tud, Swa and Wey sites first)
lcor_snpclone_tab <- lcor_snpclone_tab[which(!rownames(lcor_snpclone_tab) %in% c("Tud","Swa","Wey")), ]
(lcor_snpclone_stats <- diversity_ci(popsub(lcor_snpclone, blacklist = c("Tud","Swa","Wey")), n = 1000L, ci = 95, rarefy = TRUE, n.rare = min(summary(lcor_snpclone$pop)), raw = FALSE))
(lcor_snpclone_ParetoB <- diversity_stats(lcor_snpclone_tab, ParetoB = Beta))

# Create a clone-corrected snpclone object containing genets only
lcor_genets <- clonecorrect(lcor_snpclone)
summary(lcor_genets$pop)

# ---------- #
# Lithothamion corallioides genetic diversity ####
# ---------- #

# Basic stats on ramets (full dataset)
lcor_ramets_genind <- gl2gi(gl.filter.monomorphs(as(lcor_snpclone, "genlight")))
ramets_diversity <- basic.stats(lcor_ramets_genind)
(ramets_Ho <- apply(ramets_diversity$Ho, 2, \(x) round(mean(x), digits=2)))
(ramets_Hs <- apply(ramets_diversity$Hs, 2, \(x) round(mean(x), digits=2)))
(ramets_Fis <- apply(ramets_diversity$Fis, 2, \(x) round(mean(x, na.rm=TRUE), digits=2)))

# Basic stats on genets
lcor_genets_genind <- gl2gi(gl.filter.monomorphs(as(lcor_genets, "genlight")))
genets_diversity <- basic.stats(lcor_genets_genind)
(genets_Ho <- apply(genets_diversity$Ho, 2, \(x) round(mean(x), digits=2)))
(genets_Hs <- apply(genets_diversity$Hs, 2, \(x) round(mean(x), digits=2)))
(genets_Fis <- apply(genets_diversity$Fis, 2, \(x) round(mean(x, na.rm=TRUE), digits=2)))

# Plot genets versus ramets
plot(ramets_Ho, genets_Ho)
plot(ramets_Hs, genets_Hs)
plot(ramets_Fis, genets_Fis)

# ---------- #
# Lithothamion corallioides genetic differentiation ####
# ---------- #

# Compute pairwise Fst on ramets
summary(lcor_ramets_genind$pop)
# ramets_Fst <- hierfstat::genet.dist(lcor_ramets_genind, method = "WC84")
# save(ramets_Fst, file = "outputs/Lcorallioides_ramets_Fst.RData")
load("outputs/Lcorallioides_ramets_Fst.RData")
round(ramets_Fst, 3)

# Compute pairwise Fst on genets (not recommended based on Meirmans 2024 paper)
summary(lcor_genets_genind$pop)
# genets_Fst <- hierfstat::genet.dist(lcor_genets_genind, method = "WC84")
# save(genets_Fst, file = "outputs/Lcorallioides_genets_Fst.RData")
load("outputs/Lcorallioides_genets_Fst.RData")
round(genets_Fst, 3)

# ---------- #
# Lithothamion corallioides population structure ramets ####
# ---------- #

# Export SNPs as VCF (must download plink executable first)
# https://www.cog-genomics.org/plink/
lcor_ramets_genind
summary(lcor_ramets_genind$pop)
plink_exe <- "C:/Users/tj311/plink_win64_20250819"
lcor_ramets_genind |>
  gi2gl() |>
  gl2vcf(plink_path = plink_exe, outfile = "Lcorallioides_SNPs_ramets", outpath = "outputs/")

# Execute PopCluster externally
# First convert VCF to DAT file using PopCluster

# Read in PopCluster statistics file
file_K <- "outputs/PopCluster/Lcorallioides_ramets/Lcorallioides_ramets.K"
stats <- read.table(text = readLines(file_K, n = 9), header = TRUE)

# Shorten best run column
stats$BestRun <- str_extract(stats$BestRun, "[^\\\\]+$")
stats

# Plot PopCluster stats
png("supporting_figures/Lcorallioides_PopClusterStats.png", width = 9, height = 5, units = "in", res = 300)
par(mfrow = c(1, 2))
plot(x = stats$K, y = as.numeric(stats$DLK2), xlab = "K", ylab = "DLK2")
plot(x = stats$K, y = as.numeric(stats$FST.FIS), xlab = "K", ylab = "FST.FIS")
dev.off()

# Read in admixture proportions: best run for each K
source("misc/read_popcluster_results.R")
files <- "outputs/PopCluster/Lcorallioides_ramets/"
admix_ls <- lapply(stats$BestRun, function(i) readPopCluster(files, i, lcor_ramets_genind))
admix_ls <- lapply(1:length(admix_ls), \(x) admix_ls[[x]][which(!admix_ls[[x]]$Site == "Tud"), ])
names(admix_ls) <- paste0("K", stats$K)

# K to plot
K <- 5

# Add admixture results to new data.frame
admix_K <- admix_ls[[K]]

# Order individuals in each site by their highest cluster contributions
admix_K <- admix_K |> 
  rowwise() |> 
  mutate(
    max_cluster = which.max(c_across(3:ncol(admix_K))),
    max_value = max(c_across(3:ncol(admix_K)))
  ) |> 
  ungroup() |> 
  arrange(Site, max_cluster, desc(max_value)) |> 
  select(-max_cluster, -max_value)

# Order site labels
site_order <- c("Mil1","Mil2","Hel","Maw","Aus","Swa","Wey")

# Colour palette
col_lcor <- c("#FDDAEC","#B3CDE3","#FBB4AE","#FFFFCC","#CCEBC5","#FED9A6","#DECBE4","#F2F2F2","#E5D8BD")

# Plot STRUCTURE barplot
(plt_structure_ramets <- structure_plot(
  admixture_df = admix_K,
  cluster_cols = col_lcor,
  site_order = site_order,
  site_labels_size = 4,
  site_ticks = FALSE,
  site_dividers = TRUE, divider_width = 0.5,
  site_labels_y = -0.08,
  ylabel = "Ancestry proportion"
))

# Theme to remove all y axis content
yaxis <- theme(
  # axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
)

# Plot all structure plots
plts <- lapply(2:length(admix_ls), function(i) {

    admix_df <- data.frame(
      Site = admix_ls[[i]]$Site,
      Ind = admix_ls[[i]]$Individual,
      admix_ls[[i]][, 3:ncol(admix_ls[[i]])]
    )
  
    # Order individuals in each site by their highest cluster contributions
    admix_df <- admix_df |> 
      rowwise() |> 
      mutate(
        max_cluster = which.max(c_across(3:ncol(admix_df))),
        max_value = max(c_across(3:ncol(admix_df)))
      ) |> 
      ungroup() |> 
      arrange(Site, max_cluster, desc(max_value)) |> 
      select(-max_cluster, -max_value)
      
    if (i != 7) {
      # Hide x axis site labels
      structure_plot(admixture_df = admix_df, cluster_cols = col_lcor, site_order = site_order, legend = "none", site_ticks = F, display_site_labels = F, divider_width = 0.5, ylabel = paste0("K", i))+ yaxis
    } else {
      # Show x axis site labels on last plot
      structure_plot(admixture_df = admix_df, cluster_cols = col_lcor, site_order = site_order, legend = "none", site_ticks = F, display_site_labels = T, site_labels_y = -0.15, site_labels_size = 2.5, divider_width = 0.5, ylabel = paste0("K", i))+ yaxis
    }
  }
)

# Export as supporting figure
png("supporting_figures/Lcorallioides_PopCluster_K2_K8.png", height = 12, width = 8, units = "in", res = 600)
x <- wrap_plots(plts, ncol = 1) & theme(plot.margin = margin(l = 15, b = 0, t = 0, r = 5))
wrap_elements(x)+ theme(plot.margin = margin(l = 10, b = 15, r = 10, t = 10))
dev.off()

# Read in, filter and prepare L. corallioides coordinates for plotting
coord_lcor <- read.csv("data/site_coordinates3.csv") |> 
  filter(Species == "Lcorallioides") |> 
  dplyr::select(!Species) |> 
  filter(Site != "AusI", Site != "Gri") |> 
  mutate(Site = str_replace(Site, "AusII", "Aus")) |>
  dplyr::select(Site, Latitude, Longitude)
coord_lcor

# Edit Mil1 and Mil2 coordinates for plotting pie charts on map
coord_lcor[coord_lcor$Site == "Mil1", ]$Longitude <- -5.25
coord_lcor[coord_lcor$Site == "Mil2", ]$Longitude <- -5.05
coord_lcor[coord_lcor$Site == "Maw", ]$Latitude <- 50.20
coord_lcor[coord_lcor$Site == "Hel", ]$Latitude <- 50.07

# Create data.frame for labels coordinates
coords_lcor_labs <- coord_lcor |> 
  mutate(Lat_labs = ifelse(Site != "Hel" & Site != "Mil1", Latitude+0.13, Latitude)) |> 
  mutate(Lon_labs = ifelse(Site != "Hel" & Site != "Mil1", Longitude, Longitude-0.25))

# Admixture map: England and Wales
(lcor_admix_map <- mapmixture(
  admixture_df = admix_K, coords_df = coord_lcor,
  basemap = rnaturalearthhires::countries10[, "geometry"],
  boundary = c(xmin = -6.0, xmax = -0.5, ymin = 49.90, ymax = 52.01),
  crs = 4326,
  pie_size = 0.15, cluster_cols = col_lcor[1:K],
  scalebar_position = "br", arrow_position = "br", scalebar_size = 1.5, arrow_size = 1.5,
  pie_border = 0.10, plot_title = "Lithothamnion corallioides"
)+
  geom_label(
    data = coords_lcor_labs,
    aes(x = Lon_labs, y = Lat_labs, label = Site),
    size = 4
  )+
  # annotate("shadowtext", x = textX, y = textY, label = textLab, size = 4, colour = "black", bg.color = "white")+
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold.italic"),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(colour = "#f0f0f0"),
  )
)

# Lcorallioides structure plot
(lcor_structure_plt <- structure_plot(
  admixture_df = admix_K, cluster_cols = col_lcor,
  site_order = site_order, legend = "top",
  site_ticks = F, site_labels_size = 4, site_labels_y = -0.15
)+
  theme(
    axis.title.y = element_text(size = 12, hjust = 0.6),
    legend.text = element_text(size = 12, margin = margin(l = 10, r = 10)),
  )+
  guides(fill = guide_legend(nrow = 1))
)

# Compose figure
layout <- "
  A
  A
  B
  B
  C
"
wrap_plots(
  lcor_MLLs+ labs(tag = "A"),
  lcor_admix_map+ labs(tag = "B"),
  lcor_structure_plt+ labs(tag = "C"),
  design = layout)
ggsave("figures/Figure_04.png", device = png, type = "cairo", width = 11.5, height = 15, units = "in", dpi = 600)
ggsave("figures/Figure_04.pdf", width = 11.5, height = 15, units = "in")

# Colour palette
library(randomcoloR)
palette_lcor <- distinctColorPalette(nPop(lcor_genind))

# Principal component analysis
pca_lcor <- dudi.pca(popsub(lcor_ramets_genind, exclude = "Tud"), scannf = FALSE, nf = 10)

# Visualise PC1 vs. PC2
(pca_lcor1_plt <- scatter_plot(
  dataframe = pca_lcor$li,
  group_ids = popsub(lcor_ramets_genind, exclude = "Tud")$pop,
  labels = popsub(lcor_ramets_genind, exclude = "Tud")$pop,
  type = "labels", axes = c(1,2), xlab = "PC", ylab = "PC", size = 3, colours = palette_lcor,
  percent = round(pca_lcor$eig/sum(pca_lcor$eig)*100, 1)
))

# Visualise PC2 vs. PC3
(pca_lcor2_plt <- scatter_plot(
  dataframe = pca_lcor$li,
  group_ids = popsub(lcor_ramets_genind, exclude = "Tud")$pop,
  labels = popsub(lcor_ramets_genind, exclude = "Tud")$pop,
  type = "labels", axes = c(2,3), xlab = "PC", ylab = "PC", size = 3, colours = palette_lcor,
  percent = round(pca_lcor$eig/sum(pca_lcor$eig)*100, 1)
))

# Visualise PC3 vs. PC4
(pca_lcor3_plt <- scatter_plot(
  dataframe = pca_lcor$li,
  group_ids = popsub(lcor_ramets_genind, exclude = "Tud")$pop,
  labels = popsub(lcor_ramets_genind, exclude = "Tud")$pop,
  type = "labels", axes = c(3,4), xlab = "PC", ylab = "PC", size = 3, colours = palette_lcor,
  percent = round(pca_lcor$eig/sum(pca_lcor$eig)*100, 1)
))

# Combine plots
lcor_pcas <- wrap_plots(pca_lcor1_plt, pca_lcor2_plt, pca_lcor3_plt, ncol = 2, guides = "collect")&
  theme(legend.position = "top", legend.text = element_text(size = 14))&
  guides(fill = guide_legend(nrow = 1, override.aes = aes(label = "")))
lcor_pcas

# Export as supporting figure
ggsave(plot = lcor_pcas, filename = "supporting_figures/Lcorallioides_PCA.png", width = 10, height = 10, units = "in", dpi = 600)
