# =========================== #
#
# Maerl Whole Genome Re-sequencing Project 2024
#
# SNP Population Structure Analysis
#
# Species:
# Phymatolithon calcareum
# Lithothamnion corallioides
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
library(RColorBrewer)
library(scales)
library(patchwork)
library(sf)
library(ggplot2)
library(shadowtext)

# Diverging colour palette
diverging_palette <- brewer.pal(9, "Paired")
diverging_palette <- brewer.pal(9, "Set3")
scales::show_col(diverging_palette)

# ----------------- #
# Prepare genotypes for PopCluster software ####
# ----------------- #

# Read in VCF
pcal_vcf <- read.vcfR("outputs/pcalcareum_SNPs.vcf.gz")

# Export only sites from Britain (including Zara Shoal)
britain <- which(str_detect(colnames(pcal_vcf@gt), "Tre|Mor|Bor|Ons", negate = TRUE))
britain_vcf <- pcal_vcf[, britain]
sort(colnames(britain_vcf@gt))
britain_vcf_poly <- britain_vcf[which(is.polymorphic(britain_vcf)),]
vcfR::write.vcf(britain_vcf_poly, file = "outputs/pcalcareum_SNPs_Britain.vcf.gz")

# # Export filtered genotypes for Phymatolithon calcareum with sample size >= 3 (not AusII or Swa)
# remove_less3 <- which(str_detect(colnames(pcal_vcf@gt), "AusII+|Swa", negate = TRUE))
# remove_less3_vcf <- pcal_vcf[, remove_less3]
# sort(colnames(remove_less3_vcf@gt))
# remove_less3_vcf_poly <- remove_less3_vcf[which(is.polymorphic(remove_less3_vcf)),]
# vcfR::write.vcf(remove_less3_vcf_poly, file = "outputs/pcalcareum_SNPs_3greater_samples.vcf.gz")
#
# # Export filtered genotypes for non-coarse Phymatolithon calcareum only
# non_coarse <- which(str_detect(colnames(pcal_vcf@gt), "Maw11C+|Maw22C+", negate = TRUE))
# non_coarse_vcf <- pcal_vcf[, non_coarse]
# colnames(non_coarse_vcf@gt)
# non_coarse_vcf_poly <- non_coarse_vcf[which(is.polymorphic(non_coarse_vcf)),]
# vcfR::write.vcf(non_coarse_vcf_poly, file = "outputs/pcalcareum_SNPs_NonCoarse.vcf.gz")

# # Read in VCF
# lcor_vcf <- read.vcfR("outputs/lcorallioides_SNPs.vcf.gz")

# # Export filtered genotypes for Lithothamnion corallioides with sample size >= 3 (not Tud)
# remove_less3 <- which(str_detect(colnames(vcf_ld_bi_dp_poly@gt), "Tud", negate = TRUE))
# remove_less3_vcf <- vcf_ld_bi_dp_poly[, remove_less3]
# sort(colnames(remove_less3_vcf@gt))
# remove_less3_vcf_poly <- remove_less3_vcf[which(is.polymorphic(remove_less3_vcf)),]
# vcfR::write.vcf(remove_less3_vcf_poly, file = "./outputs/lcorallioides_SNPs_3greater_samples.vcf.gz")


# ----------------- #
# Function to read in PopCluster results
# ----------------- #

readPopCluster <- function(path, file, genind) {
  
  # Read the lines of the text file
  lines <- readLines(paste0(path, file))
  
  # Find the index of the marker that indicates the start of the table
  start_index <- grep("Inferred ancestry of individuals", lines)
  
  # Find the index of the first empty line after the start of the table
  end_index <- start_index + which(lines[start_index:length(lines)] == "")[1] - 1
  
  # Skip lines before the start of the table and read the table into a data frame
  table_data <- read.table(text = lines[(start_index + 2):end_index], header = FALSE, fill = TRUE)
  
  # Subset columns
  admix_data <- subset(table_data, select = c(3,7:ncol(table_data)))
  colnames(admix_data) <- c("Individual", paste(rep("Cluster"), 1:(ncol(admix_data)-1)))
  
  # Add a column for site labels
  admix_data$Site <- as.character(genind$pop)
  
  admix_data <- dplyr::select(admix_data, Site, Individual, everything())
  
  # Return table
  return(admix_data)
}


# ----------------- #
# Phymatolithon calcareum: PopCluster ####
# ----------------- #

# Read in Pcalcareum SNPs
pcal_genind <- vcfR2genind(read.vcfR("outputs/pcalcareum_SNPs_Britain.vcf.gz"))

# Merge the following samples into the same group (population)
# • St Mawes coarse maerl from 2011 and 2022
# • St Mawes maerl from 2015 and 2022
# • All samples from St Austell area (St Austell Bay and Little Gribbin)
pcal_genind$pop <- as.factor(sub("\\_.*", "", indNames(pcal_genind)))
pcal_genind$pop <- as.factor(str_replace(pcal_genind$pop, "Maw22C|Maw11C", "MawC"))
pcal_genind$pop <- as.factor(str_replace(pcal_genind$pop, "Maw15|Maw22", "Maw"))
pcal_genind$pop <- as.factor(str_replace(pcal_genind$pop, "Gri|AusII", "Aus"))
summary(pcal_genind$pop)

# Read in statistics file
file <- "outputs/Pcalcareum_PopCluster_Britain_Admixture_Unequal_MediumScaling/"
file_K <- "Unequal_Medium.K"
stats <- read.table(text = readLines(paste0(file, file_K), n = 9), header = TRUE)

# Shorten best run column
stats$BestRun <- str_extract(stats$BestRun, "[^\\\\]+$")
stats

# Plot stats
png("supporting_figures/PopClusterStats_PCAL.png", width = 600, height = 350)
par(mfrow = c(1, 2))
plot(x = stats$K, y = as.numeric(stats$DLK2), xlab = "K", ylab = "DLK2")
plot(x = stats$K, y = as.numeric(stats$FST.FIS), xlab = "K", ylab = "FST.FIS")
dev.off()

# Read in admixture proportions: best run for each K
admix_df <- lapply(stats$BestRun, function(i) readPopCluster(file, i, pcal_genind))
names(admix_df) <- paste0("K", stats$K)

# Order site labels
site_order <- c("Zar","Man","Maw","MawC","Biz","Ger","Nar","Aus","Swa","Wey")

# Colour palettes
col_britain <- c(cluster1="#80B1D3", cluster2="#FDB462", cluster3="#FCCDE5", cluster4="#B3DE69", cluster5="#D9D9D9")
colours <- c("#E17E68","#FDB462","#80B1D3","#B3DE69","#8DD3C7","#BEBADA","#D9D9D9","#FFFFB3","#FCCDE5")

# Plot single structure plot
structure_plot(admix_df$K5, cluster_cols = col_britain, site_order = site_order, legend = "right")
structure_plot(admix_df$K5, labels = "individual", cluster_cols = col_britain, site_order = site_order, legend = "right")

# Theme to remove all y axis content
yaxis <- theme(
  # axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
)

# Plot all structure plots
plts <- lapply(1:length(admix_df[2:8]), function(i) {
      
    if (i != 7) {
      # Hide x axis site labels
      structure_plot(admixture_df = admix_df[2:8][[i]], cluster_cols = colours, site_order = site_order, legend = "none", site_ticks = F, display_site_labels = F, divider_width = 0.5, ylabel = paste0("K", i+1))+ yaxis
    } else {
      # Show x axis site labels on last plot
      structure_plot(admixture_df = admix_df[2:8][[i]], cluster_cols = colours, site_order = site_order, legend = "none", site_ticks = F, display_site_labels = T, site_labels_y = -0.15, site_labels_size = 2.5, divider_width = 0.5, ylabel = paste0("K", i+1))+ yaxis
    }
  }
)

# Export as supporting figure
png("supporting_figures/Pcalcareum_Britain_structure_K2_K8.png", height = 12, width = 8, units = "in", res = 600)
x <- wrap_plots(plts, ncol = 1) & theme(plot.margin = margin(l = 15, b = 0, t = 0, r = 5))
wrap_elements(x)+ theme(plot.margin = margin(l = 10, b = 15, r = 10, t = 10))
dev.off()

# Read in P. calcareum coordinates
coord_pcal <- subset(read.csv("data/site_coordinates2.csv"), Species == "Pcalcareum")

# Read in, filter and prepare P. calcareum coordinates for plotting
coord_pcal <- read.csv("data/site_coordinates2.csv") |> 
  filter(Species == "Pcalcareum") |> 
  filter(Site != "Gri", Site != "Mor", Site != "Tre", Site != "Bor", Site != "Ons", Site != "Swa", Site != "Bem") |> 
  mutate(Site = str_replace(Site, "AusII", "Aus")) |> 
  dplyr::select(!Species) |> 
  dplyr::select(Site, Latitude, Longitude)
coord_pcal

# Read in high resolution England coastline geopackage
england <- st_read("data/england_ol_2001.gpkg")

# Text labels on map
textX <- c(-4.95, -5.0, -4.67)
textY <- c(50.04, 50.11, 50.29)
textLab <- c("The Manacles", "Falmouth", "St Austell")

# Subset admixture data frame
admix_df_plotting <- admix_df$K5 |> filter(Site != "Swa")

# Admixture map: England
britain <- c(xmin = -5.75, xmax = -4.3, ymin = 49.95, ymax = 50.35)
Fig5A <- mapmixture(
  admixture_df = admix_df_plotting, coords_df = coord_pcal,
  basemap = england[, "geom"], basemap_border_lwd = 0.05,
  boundary = britain, crs = 4326, plot_title = "Phymatolithon calcareum",
  pie_size = 0.025, cluster_cols = col_britain,
  pie_border = 0.10, scalebar_size = 1.9, arrow_size = 1.8
)+
  annotate("text", x = -5.27, y = 50.16, label = "Cornwall", size = 5)+
  annotate("shadowtext", x = textX, y = textY, label = textLab, size = 4, colour = "black", bg.color = "white")+
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold.italic"),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(colour = "#f0f0f0"),
  )
# Fig5A

# Inset text labels on map
insettextX <- c(-0.60, -2.90)
insettextY <- c(51.2, 54.0)
insettextLab <- c("Weymouth", "Strangford Lough")

# Admixture map: Britain
admix_df_inset <- filter(admix_df_plotting, Site %in% c("Wey","Zar"))
coords_df_inset <- filter(coord_pcal, Site %in% c("Wey","Zar"))
# inset_box <- transform_bbox(britain, st_crs(england))
inset_box <- britain
Fig5A_inset <- mapmixture(
  admixture_df = admix_df_inset, coords_df = coords_df_inset,
  basemap = rnaturalearthhires::countries10[, "geometry"], basemap_border_lwd = 0.05,
  boundary = c(xmin = -11.50, xmax = 2.50, ymin = 49.00, ymax = 55.10),
  # crs = st_crs(england),
  crs = 4326,
  pie_size = 0.7, cluster_cols = col_britain,
  scalebar_position = "bl", arrow = F, scalebar_size = 1.5,
  pie_border = 0.10,
)+
  annotate("shadowtext", x = insettextX, y = insettextY, label = insettextLab, size = 3.5, colour = "black", bg.color = "white")+
  annotate("rect", xmin = inset_box["xmin"], xmax = inset_box["xmax"], ymin = inset_box["ymin"], ymax = inset_box["ymax"], colour = "black", alpha = 0.50, lwd = 0.2)+
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
# Fig5A_inset

# Inset plot
Fig5A_final <- Fig5A+
  inset_element(Fig5A_inset, left = 0.65, bottom = 0.05, right = 0.98, top = 0.7, align_to = "full")
# Fig5A_final

# Pcalcarem K5 structure plot
Fig5B <- structure_plot(
  admixture_df = admix_df$K5, cluster_cols = col_britain,
  site_order = site_order, legend = "top",
  site_ticks = F, site_labels_size = 3.5, site_labels_y = -0.15
)+
  theme(
    axis.title.y = element_text(size = 8, hjust = 0.6),
  )+
  guides(fill = guide_legend(nrow = 1))
# Fig5B

# Combine admixture map and structure plot
layout <- "
  A
  A
  A
  B
"
wrap_plots(wrap_elements(Fig5A_final), Fig5B, design = layout)
ggsave("figures/Figure_03.png", device = png, type = "cairo", width = 12, height = 8, units = "in", dpi = 600)
ggsave("figures/Figure_03.pdf", width = 12, height = 8, units = "in")


# ----------------- #
# Lithothamnion corallioides: PopCluster ####
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

# Read in statistics file
file <- "outputs/Lcorallioides_PopCluster_Britain_Admixture_Unequal_MediumScaling//"
file_K <- "Unequal_Medium.K"
stats <- read.table(text = readLines(paste0(file, file_K), n = 9), header = TRUE)

# Shorten best run column
stats$BestRun <- str_extract(stats$BestRun, "[^\\\\]+$")
stats

# Plot stats
png("supporting_figures/PopClusterStats_LCOR.png", width = 600, height = 350)
par(mfrow = c(1, 2))
plot(x = stats$K, y = as.numeric(stats$DLK2), xlab = "K", ylab = "DLK2")
plot(x = stats$K, y = as.numeric(stats$FST.FIS), xlab = "K", ylab = "FST.FIS")
dev.off()

# Read in admixture proportions: best run for each K
admix_df_lcor <- lapply(stats$BestRun, function(i) readPopCluster(file, i, lcor_genind))
names(admix_df_lcor) <- paste0("K", stats$K)

# Order site labels
site_order <- c("Mil1","Mil2","Tud","Hel","Maw","Aus","Swa","Wey")

# Colour palette
scales::show_col(brewer.pal(9, "Pastel1"))
# Cluster order: Maw&Hel, Hel, Mil, MilIso, AusMix, Swa&Wey
col_lcor <- c("#FDDAEC","#B3CDE3","#FBB4AE","#FFFFCC","#CCEBC5")
col_lcor_all <- c("#FDDAEC","#B3CDE3","#FBB4AE","#FFFFCC","#CCEBC5","#FED9A6","#DECBE4","#F2F2F2","#E5D8BD")

# Plot single structure plot
structure_plot(admix_df_lcor$K5, cluster_cols = col_lcor, labels = "individual", site_order = site_order, legend = "right")

# Theme to remove all y axis content
yaxis <- theme(
  # axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
)

# Plot all structure plots
plts <- lapply(1:length(admix_df_lcor[2:8]), function(i) {
  
    if (i != 7) {
      # Hide x axis site labels
      structure_plot(admix_df_lcor[2:8][[i]], cluster_cols = col_lcor_all, site_order = site_order, legend = "none", site_ticks = F, display_site_labels = F, divider_width = 0.5, ylabel = paste0("K", i+1))+ yaxis
    } else {
      # Show x axis site labels on last plot
      structure_plot(admix_df_lcor[2:8][[i]], cluster_cols = col_lcor_all, site_order = site_order, legend = "none", site_ticks = F, display_site_labels = T, site_labels_y = -0.15, site_labels_size = 2.5, divider_width = 0.5, ylabel = paste0("K", i+1))+ yaxis
    }
  }
)

# Export as supporting figure
png("supporting_figures/Lcorallioides_structure_K2_K9.png", height = 12, width = 8, units = "in", res = 600)
x <- wrap_plots(plts, ncol = 1) & theme(plot.margin = margin(l = 15, b = 0, t = 0, r = 5))
wrap_elements(x)+ theme(plot.margin = margin(l = 10, b = 10, r = 10, t = 10))
dev.off()

# Lcorallioides K6 structure plot
Fig6B <- structure_plot(
  admixture_df = admix_df_lcor$K5, cluster_cols = col_lcor,
  site_order = site_order, legend = "top",
  site_ticks = F, site_labels_size = 3.5, site_labels_y = -0.15
)+
  theme(
    axis.title.y = element_text(size = 8, hjust = 0.6),
  )+
  guides(fill = guide_legend(nrow = 1))
Fig6B

# Read in, filter and prepare L. corallioides coordinates for plotting
coord_lcor <- read.csv("data/site_coordinates3.csv") |> 
  filter(Species == "Lcorallioides") |> 
  dplyr::select(!Species) |> 
  filter(Site != "AusI", Site != "Gri") |> 
  mutate(Site = str_replace(Site, "AusII", "Aus")) |>
  filter(Site != "Mil2") |> 
  mutate(Site = str_replace(Site, "Mil1", "Mil")) |>
  dplyr::select(Site, Latitude, Longitude)
coord_lcor

# Text labels on map
textX <- c(-4.69, -4.45, -5.40, -2.4, -1.8)
textY <- c(50.09, 50.21, 51.6, 50.48, 50.53)
textLab <- c("Falmouth", "St Austell", "Milford Haven", "Weymouth", "Swanage")

# Subset admixture data frame
admix_df_plotting <- admix_df_lcor$K5 |>
  mutate(Site = str_replace(Site, "Mil2|Mil1", "Mil")) |> 
  filter(Site != "Tud")

# Admixture map: England and Wales
Fig6A <- mapmixture(
  admixture_df = admix_df_plotting, coords_df = coord_lcor,
  basemap = rnaturalearthhires::countries10[, "geometry"],
  boundary = c(xmin = -6.0, xmax = -0.5, ymin = 49.90, ymax = 52.05),
  crs = 4326,
  pie_size = 0.15, cluster_cols = col_lcor,
  scalebar_position = "br", arrow_position = "br", scalebar_size = 1.5, arrow_size = 1.5,
  pie_border = 0.10, plot_title = "Lithothamnion corallioides"
)+
  annotate("shadowtext", x = textX, y = textY, label = textLab, size = 4, colour = "black", bg.color = "white")+
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold.italic"),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(colour = "#f0f0f0"),
  )
# Fig6A

# Combine admixture map and structure plot
layout <- "
  A
  A
  A
  B
"
wrap_plots(Fig6A, Fig6B, design = layout)
ggsave("figures/Figure_04.png", device = png, type = "cairo", width = 12, height = 8, units = "in", dpi = 600)
ggsave("figures/Figure_04.pdf", width = 12, height = 8, units = "in")

# ----------------- #
# Figure 5: Composite ####
# ----------------- #

# # Layout
# layout <- "
#   A
#   A
#   A
#   B
#   C
#   C
#   C
#   D
# "
# 
# # Plot
# wrap_plots(wrap_elements(Fig5A_final), Fig5B, Fig5C, Fig5D, design = layout)
# ggsave("Figure_05.png", device = png, type = "cairo", width = 15, height = 8, units = "in", dpi = 600)
