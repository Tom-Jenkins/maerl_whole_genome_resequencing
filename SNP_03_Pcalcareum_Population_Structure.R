# ========== #
#
# Maerl WGS Analysis 2025
#
# Clonality, Genetic Diversity and Population Structure
#
# Species:
# Phymatolithon calcareum
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
# library(shadowtext)

# ---------- #
# Phymatolithon calcareum clonal diversity ####
# ---------- #

# Read in Pcalcareum SNPs
pcal_genind <- vcfR2genind(read.vcfR("outputs/pcalcareum_SNPs.vcf.gz"))
pcal_genind

# Merge or remove sample sites
# • St Mawes coarse maerl from 2011 and 2022
# • St Mawes maerl from 2015 and 2022
# • All samples from St Austell area (St Austell Bay and Little Gribbin)
# • Remove geographical samples from Europe from Jenkins et al. 2021
pcal_genind$pop <- as.factor(sub("\\_.*", "", indNames(pcal_genind)))
pcal_genind$pop <- as.factor(str_replace(pcal_genind$pop, "Maw22C|Maw11C", "MawC"))
pcal_genind$pop <- as.factor(str_replace(pcal_genind$pop, "Maw15|Maw22", "Maw"))
pcal_genind$pop <- as.factor(str_replace(pcal_genind$pop, "Gri|AusII", "Aus"))
pcal_genind$pop <- as.factor(str_replace(pcal_genind$pop, "Her|Nar", "Ger"))
pcal_genind <- popsub(pcal_genind, exclude = c("Bor","Mor","Ons","Tre","Swa"))
summary(pcal_genind$pop)

# Check and remove any monomorphic SNPS after filtering individuals and sites
pcal_genind <- pcal_genind[loc=names(which(isPoly(pcal_genind)))] 

# Convert genind to snpclone
pcal_snpclone <- as.snpclone(gi2gl(pcal_genind))

# Number of multi-locus genotypes (MLGs)
# Different from the multi-locus lineages (MLL) / clonal lineage / genet
mlg(pcal_genind)

# See Arnaud-Haond et al. 2007
# https://onlinelibrary.wiley.com/doi/10.1111/j.1365-294X.2007.03535.x#ss3-title

# Calculate Hamming genetic distance
pcal_hamming = bitwise.dist(pcal_snpclone, mat = FALSE, percent = TRUE)
hist(pcal_hamming, breaks = 100, xlab = "Hamming genetic distance", main = "Hamming genetic distance")

# Clone cut-off threshold for defining MLL
max_distance <- max(as.vector(pcal_hamming))
pcal_threshold <- mlg.filter(pcal_snpclone, distance = pcal_hamming, stats = "THRESHOLDS", threshold = max_distance, algorithm = "average_neighbor")
(pcal_clone_threshold <- cutoff_predictor(pcal_threshold))

# Plot histogram of Hamming distance values
png("supporting_figures/Pcalcareum_hamming_histogram.png", width = 10, height = 5, units = "in", res = 300)
hist(pcal_hamming, breaks = 100, xlab = "Hamming genetic distance", main = "Hamming genetic distance")
abline(v = pcal_clone_threshold, col = "red")
dev.off()

# Filter snpclone object to get the multi-locus lineages (MLLs)
# Proxy of the number of genets
mlg.filter(pcal_snpclone, distance = pcal_hamming, algorithm = "average_neighbor") <- pcal_clone_threshold
pcal_snpclone
mlg(pcal_snpclone)

# Table of clonal lineages
pcal_snpclone_tab <- mlg.table(pcal_snpclone)

# Prepare data.frame in long format, filter out zeros and order by MLL number
pcal_snpclone_tab_long <- pcal_snpclone_tab |> 
  as.data.frame() |> 
  cbind(Site = rownames(pcal_snpclone_tab)) |> 
  pivot_longer(
    cols = starts_with("MLG"), 
    names_to = "MLG", 
    values_to = "Count",
  ) |> 
  dplyr::filter(Count != 0) |> 
  mutate(MLG = str_replace(MLG, "MLG.", "")) |> 
  mutate(MLG = factor(MLG, levels = sort(unique(as.integer(MLG)))))

# Check whether MLLs cross populations
mlg.crosspop(pcal_snpclone)

# Add column that represents whether a clonal lineage overlaps with more than one site
(overlapping_mll <- str_remove(names(mlg.crosspop(pcal_snpclone)), "MLG."))
pcal_snpclone_tab_long <- pcal_snpclone_tab_long |>
  mutate(
    Overlap = case_when(
      MLG == as.integer(overlapping_mll[1]) ~ overlapping_mll[1],
      MLG == as.integer(overlapping_mll[2]) ~ overlapping_mll[2],
      TRUE      ~ "Site-Specific Lineage"
    )
  )

# Plot MLLs
(pcal_MLLs <- ggplot(data = pcal_snpclone_tab_long)+
  geom_bar(aes(x=MLG, y=Count, fill=Overlap), stat = "identity", linewidth = 0.1)+
  scale_fill_manual(values = c("#90a4ae","#a1887f","grey70"))+
  facet_wrap(~Site, scales = "free", ncol = 4)+
  scale_y_continuous(
    limits = c(0, max(pcal_snpclone_tab_long$Count)),
    breaks = c(0, 5, 10, 13)
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
mawc <- popsub(pcal_snpclone, sublist = "MawC")
mlg19 <- mawc[mawc@mlg == 70, ]
mlg19_different_alleles <- bitwise.dist(mlg19, mat = FALSE, percent = FALSE)
mean(mlg19_different_alleles); median(mlg19_different_alleles); range(mlg19_different_alleles)

# Total differences per individual
total_diff <- rowSums(as.matrix(mlg19_different_alleles))

# Find the "consensus" individual, i.e. genetically smallest distance to all others
consensus_idx <- which.min(total_diff)

# Extract differences from consensus clone
# I.e. candidate somatic mutations
diff_to_consensus <- as.matrix(mlg19_different_alleles)[, consensus_idx]
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

# Clonal diversity statistics
(pcal_snpclone_stats <- diversity_ci(pcal_snpclone_tab, n = 1000L, ci = 95, rarefy = TRUE, n.rare = min(summary(pcal_snpclone$pop)), raw = FALSE))
(pcal_snpclone_ParetoB <- diversity_stats(pcal_snpclone_tab, ParetoB = Beta))

# Create a clone-corrected snpclone object containing genets only
pcal_genets <- clonecorrect(pcal_snpclone)
summary(pcal_genets$pop)

# ---------- #
# Phymatolithon calcareum genetic diversity ####
# ---------- #

# Basic stats on ramets (full dataset)
ramets_diversity <- basic.stats(pcal_genind)
apply(ramets_diversity$Ho, 2, \(x) round(mean(x), digits=2))
apply(ramets_diversity$Hs, 2, \(x) round(mean(x), digits=2))
apply(ramets_diversity$Fis, 2, \(x) round(mean(x, na.rm=TRUE), digits=2))

# Basic stats on genets
pcal_genets_genind <- gl2gi(gl.filter.monomorphs(as(pcal_genets, "genlight")))
genets_diversity <- basic.stats(pcal_genets_genind)
apply(genets_diversity$Ho, 2, \(x) round(mean(x), digits=2))
apply(genets_diversity$Hs, 2, \(x) round(mean(x), digits=2))
apply(genets_diversity$Fis, 2, \(x) round(mean(x, na.rm=TRUE), digits=2))

# ---------- #
# Phymatolithon calcareum genetic differentiation ####
# ---------- #

# Compute pairwise Fst on ramets
summary(pcal_genind$pop)
# ramets_Fst <- hierfstat::genet.dist(pcal_genind, method = "WC84")
# save(ramets_Fst, file = "outputs/Pcalcareum_ramets_Fst.RData")
load("outputs/Pcalcareum_ramets_Fst.RData")
round(ramets_Fst, 3)

# Compute pairwise Fst on genets (not recommended based on Meirmans 2024 paper)
summary(pcal_genets_genind$pop)
# genets_Fst <- hierfstat::genet.dist(pcal_genets_genind, method = "WC84")
# save(genets_Fst, file = "outputs/Pcalcareum_genets_Fst.RData")
load("outputs/Pcalcareum_genets_Fst.RData")
round(genets_Fst, 3)

# ---------- #
# Phymatolithon calcareum population structure ramets ####
# ---------- #

# Export SNPs as VCF for PopCluster input (must download plink executable first)
# https://www.cog-genomics.org/plink/
pcal_genind
summary(pcal_genind$pop)
plink_exe <- "C:/Users/tj311/plink_win64_20250819"
pcal_genind |>
  gi2gl() |>
  gl2vcf(plink_path = plink_exe, outfile = "Pcalcareum_SNPs_ramets", outpath = "outputs/")

# Execute PopCluster externally
# First convert VCF to DAT file using PopCluster

# Read in PopCluster statistics file
file_K <- "outputs/PopCluster/Pcalcareum_ramets/Pcalcareum_ramets.K"
stats <- read.table(text = readLines(file_K, n = 9), header = TRUE)

# Shorten best run column
stats$BestRun <- str_extract(stats$BestRun, "[^\\\\]+$")
stats

# Plot PopCluster stats
png("supporting_figures/Pcalcareum_PopClusterStats.png", width = 9, height = 5, units = "in", res = 300)
par(mfrow = c(1, 2))
plot(x = stats$K, y = as.numeric(stats$DLK2), xlab = "K", ylab = "DLK2")
plot(x = stats$K, y = as.numeric(stats$FST.FIS), xlab = "K", ylab = "FST.FIS")
dev.off()

# Read in admixture proportions: best run for each K
source("misc/read_popcluster_results.R")
files <- "outputs/PopCluster/Pcalcareum_ramets/"
admix_ls <- lapply(stats$BestRun, function(i) readPopCluster(files, i, pcal_genind))
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
site_order <- c("Zar","Man","Maw","MawC","Biz","Ger","Aus","Wey")

# Cluster colours
col_pcal <- c("#80B1D3","#FDB462","#FCCDE5","#B3DE69","#D9D9D9","#E17E68","#8DD3C7","#BEBADA","#FFFFB3")

# Plot STRUCTURE barplot
(plt_structure_ramets <- structure_plot(
  admixture_df = admix_K,
  cluster_cols = col_pcal,
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
      structure_plot(admixture_df = admix_df, cluster_cols = col_pcal, site_order = site_order, legend = "none", site_ticks = F, display_site_labels = F, divider_width = 0.5, ylabel = paste0("K", i))+ yaxis
    } else {
      # Show x axis site labels on last plot
      structure_plot(admixture_df = admix_df, cluster_cols = col_pcal, site_order = site_order, legend = "none", site_ticks = F, display_site_labels = T, site_labels_y = -0.15, site_labels_size = 2.5, divider_width = 0.5, ylabel = paste0("K", i))+ yaxis
    }
  }
)

# Export as supporting figure
png("supporting_figures/Pcalcareum_PopCluster_K2_K8.png", height = 12, width = 8, units = "in", res = 600)
x <- wrap_plots(plts, ncol = 1) & theme(plot.margin = margin(l = 15, b = 0, t = 0, r = 5))
wrap_elements(x)+ theme(plot.margin = margin(l = 10, b = 15, r = 10, t = 10))
dev.off()

# Read in high resolution England coastline geopackage
england <- st_read("data/england_ol_2001.gpkg")

# Read in P. calcareum site coordinates
coord_pcal <- subset(read.csv("data/site_coordinates2.csv"), Species == "Pcalcareum")

# Read in, filter and prepare P. calcareum coordinates for plotting
coord_pcal <- read.csv("data/site_coordinates2.csv") |> 
  filter(Species == "Pcalcareum") |> 
  filter(Site != "Gri", Site != "Nar", Site != "Mor", Site != "Tre", Site != "Bor", Site != "Ons", Site != "Swa", Site != "Bem") |> 
  mutate(Site = str_replace(Site, "AusII", "Aus")) |> 
  dplyr::select(!Species) |> 
  dplyr::select(Site, Latitude, Longitude)
coord_pcal

# Create data.frame for labels coordinates
coords_pcal_labs <- coord_pcal |> 
  mutate(Lat_labs = ifelse(Site != "Ger" & Site != "MawC", Latitude, Latitude+0.03)) |> 
  mutate(Lon_labs = ifelse(Site != "Ger" & Site != "MawC", Longitude-0.065, Longitude))

# Admixture map: south-west England
SWEngland <- c(xmin = -5.75, xmax = -4.3, ymin = 49.95, ymax = 50.36)
pcal_admix_map <- mapmixture(
  admixture_df = admix_K, coords_df = coord_pcal,
  basemap = england[, "geom"], basemap_border_lwd = 0.05,
  boundary = SWEngland, crs = 4326, plot_title = "Phymatolithon calcareum",
  pie_size = 0.03, cluster_cols = col_pcal[1:K],
  pie_border = 0.10, scalebar_size = 1.9, arrow_size = 1.8
)+
  geom_label(
    data = coords_pcal_labs,
    aes(x = Lon_labs, y = Lat_labs, label = Site),
    size = 4.5
  )+
  annotate("text", x = -5.27, y = 50.16, label = "Cornwall", size = 5)+
  # annotate("shadowtext", x = textX, y = textY, label = textLab, size = 4, colour = "black", bg.color = "white")+
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold.italic"),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(colour = "#f0f0f0"),
  )

# Admixture map: Britain
admix_K_inset <- filter(admix_K, Site %in% c("Wey","Zar"))
coord_pcal_inset <- filter(coord_pcal, Site %in% c("Wey","Zar"))
pcal_admix_map_inset <- mapmixture(
  admixture_df = admix_K_inset, coords_df = coord_pcal_inset,
  basemap = rnaturalearthhires::countries10[, "geometry"], basemap_border_lwd = 0.05,
  boundary = c(xmin = -11.50, xmax = 2.50, ymin = 49.00, ymax = 55.10),
  # crs = st_crs(england),
  crs = 4326,
  pie_size = 0.7, cluster_cols = col_pcal[1:K],
  scalebar_position = "bl", arrow = F, scalebar_size = 1.5,
  pie_border = 0.10,
)+
  geom_label(
    data = coord_pcal_inset,
    aes(x = ifelse(Site == "Zar", Longitude+1.45, Longitude+1.7), y = Latitude, label = Site),
    size = 4.5
  )+
  # annotate("shadowtext", x = insettextX, y = insettextY, label = insettextLab, size = 3.5, colour = "black", bg.color = "white")+
  annotate("rect", xmin = SWEngland["xmin"], xmax = SWEngland["xmax"], ymin = SWEngland["ymin"], ymax = SWEngland["ymax"], colour = "black", alpha = 0.50, lwd = 0.2)+
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    # panel.grid.major = element_line(colour = "#f0f0f0"),
    plot.background = element_blank(),
    plot.margin = margin(0,0,0,0, "pt"),
  )

# Inset plot
inset_plot <- pcal_admix_map+
  inset_element(pcal_admix_map_inset, left = 0.65, bottom = 0.05, right = 0.98, top = 0.7, align_to = "full")
# inset_plot

# Pcalcarem K5 structure plot
pcal_admix_structure <- structure_plot(
  admixture_df = admix_K, cluster_cols = col_pcal[1:5],
  site_order = site_order, legend = "top",
  site_ticks = F, site_labels_size = 4, site_labels_y = -0.15
)+
  theme(
    axis.title.y = element_text(size = 12, hjust = 0.6),
    legend.text = element_text(size = 12, margin = margin(l = 10, r = 10)),
  )+
  guides(fill = guide_legend(nrow = 1))

# Combine admixture map and structure plot
layout <- "
  A
  A
  B
  B
  C
"
wrap_plots(
  pcal_MLLs+ labs(tag = "A"),
  wrap_elements(inset_plot)+ labs(tag = "B"),
  pcal_admix_structure+ labs(tag = "C"),
  design = layout)
ggsave("figures/Figure_03.png", device = png, type = "cairo", width = 11.5, height = 15, units = "in", dpi = 600)
ggsave("figures/Figure_03.pdf", width = 11.5, height = 15, units = "in")

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

# Combine plots
pcal_pcas <- wrap_plots(pca_pcal1_plt, pca_pcal2_plt, guides = "collect")&
  theme(legend.position = "top")&
  guides(fill = guide_legend(nrow = 1, override.aes = aes(label = "")))
pcal_pcas

# Export as supporting figure
ggsave(plot = pcal_pcas, filename = "supporting_figures/Pcalcareum_PCA.png",
       width = 12, height = 7, units = "in", dpi = 600)
