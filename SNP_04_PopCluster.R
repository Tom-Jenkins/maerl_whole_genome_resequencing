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

# Read in Pcalcareum SNPs
pcal_genind <- vcfR2genind(read.vcfR("./outputs/pcalcareum_SNPs_3greater_samples.vcf.gz"))
indNames(pcal_genind) |> sort()

# Add and edit population labels
pcal_genind$pop <- as.factor(sub("\\_.*", "", indNames(pcal_genind)))
pcal_genind$pop <- as.factor(str_replace(pcal_genind$pop, "Maw22C|Maw11C", "MawC"))
pcal_genind$pop <- as.factor(str_replace(pcal_genind$pop, "Maw15|Maw22", "Maw"))
summary(pcal_genind$pop)

# Read in Lcorallioides SNPs
lcor_genind <- vcfR2genind(read.vcfR("./outputs/lcorallioides_SNPs_3greater_samples.vcf.gz"))

# Add and edit population labels
lcor_genind$pop <- as.factor(sub("\\_.*", "", indNames(lcor_genind)))
lcor_genind$pop <- as.factor(str_replace(lcor_genind$pop, "Maw15|Maw22", "Maw"))
summary(lcor_genind$pop)

# Function to read in PopCluster results
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

# Read in statistics file
file <- "outputs/Pcalcareum_PopCluster_MedScaling_Min3/"
file_K <- "Pcalcareum_PopCluster_MedScaling_Min3.K"
stats <- read.table(text = readLines(paste0(file, file_K), n = 10), header = TRUE)

# Shorten best run column
stats$BestRun <- str_extract(stats$BestRun, "(Pcalcareum_PopCluster_MedScaling_Min3_).*$")
stats

# Plot stats
plot(x = stats$K, y = as.numeric(stats$DLK2))
plot(x = stats$K, y = as.numeric(stats$FST.FIS))

# Read in admixture proportions: best run for each K
admix_df <- lapply(stats$BestRun, function(i) readPopCluster(file, i, pcal_genind))
names(admix_df) <- paste0("K", stats$K)

# Order site labels
site_order <- c("Zar","Man","Biz","Gri","Ger","Her","Maw","MawC","Wey","Mor","Tre","Bor","Ons")

# Colour palette
# Cluster order: Bor, Mawc, Zar, Wey, Tre, Mor, Her, Cornwall mix, Ons
colours <- c("#E17E68","#FDB462","#80B1D3","#B3DE69","#8DD3C7","#BEBADA","#D9D9D9","#FFFFB3","#FCCDE5")
colours2 <- c("#E17E68","#FDB462","#80B1D3","#B3DE69","#8DD3C7","#BEBADA","#D9D9D9","#FCCDE5","#FFFFB3")

# Plot single structure plot
# structure_plot(admix_df$K2, cluster_cols = colours, site_order = site_order, legend = "right")

# Theme to remove all y axis content
yaxis <- theme(
  axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
)

# Plot all structure plots
plts <- lapply(1:length(admix_df[1:8]), function(i) {
    
    if (i != 7 && i != 8) {
      structure_plot(admixture_df = admix_df[1:8][[i]], cluster_cols = colours, site_order = site_order, legend = "none", site_ticks = F, display_site_labels = F, divider_width = 0.5)+ yaxis
    } else if (i == 7) {
      structure_plot(admixture_df = admix_df[1:8][[i]], cluster_cols = colours2, site_order = site_order, legend = "none", site_ticks = F, display_site_labels = F, divider_width = 0.5)+ yaxis
    } else if (i == 8) {
      structure_plot(admixture_df = admix_df[1:8][[i]], cluster_cols = colours2, site_order = site_order, legend = "none", site_ticks = F, display_site_labels = F, divider_width = 0.5)+ yaxis
    }
  }
)
png("supporting_figures/Pcalareum_structure_K2_K9.png", height = 12, width = 8, units = "in", res = 600)
x <- wrap_plots(plts, ncol = 1) & theme(plot.margin = margin(l = 15, b = 0, t = 0, r = 5))
wrap_elements(x)+ theme(plot.margin = margin(l = 10, b = 10, r = 10, t = 10))
dev.off()

# Read in P. calcareum coordinates
coord_pcal <- subset(read.csv("data/site_coordinates3.csv"), Species == "Pcalcareum")

# Reformat for mapmixture
coord_pcal <- dplyr::select(dplyr::select(coord_pcal, !Species), Site, Latitude, Longitude)

# Read in high resolution England coastline geopackage
england <- st_read("data/england_ol_2001.gpkg")

# Text labels on map
textX <- c(-4.95, -5.0, -4.67)
textY <- c(50.04, 50.11, 50.29)
textLab <- c("The Manacles", "Falmouth", "St Austell")

# Admixture map: Britain
britain <- c(xmin = -5.75, xmax = -4.3, ymin = 49.95, ymax = 50.35)
Fig5A <- mapmixture(
  admixture_df = admix_df$K9, coords_df = coord_pcal,
  basemap = england[, "geom"],
  boundary = britain, crs = 4326, plot_title = "Phymatolithon calcareum",
  pie_size = 0.03, cluster_cols = colours2,
  pie_border = 0.10, scalebar_size = 1.9, arrow_size = 1.8
)+
  annotate("text", x = -5.27, y = 50.16, label = "Cornwall", size = 5, fontface = "bold")+
  annotate("shadowtext", x = textX, y = textY, label = textLab, size = 4, colour = "black", bg.color = "white")+
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold.italic"),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(colour = "#f0f0f0"),
  )
# Fig5A

# Inset text labels on map
insettextX <- c(-9.6, -6.50, -0.15, -2.55)
insettextY <- c(43.95, 48.2, 51.5, 53.50)
insettextLab <- c("Galicia", "Brittany", "Weymouth", "Strangford Lough")

# Admixture map: Western Europe 
admix_df_inset <- filter(admix_df$K9, Site %in% c("Ons","Bor","Tre","Mor","Wey","Zar"))
coords_df_inset <- filter(coord_pcal, Site %in% c("Ons","Bor","Tre","Mor","Wey","Zar"))
# inset_box <- transform_bbox(britain, st_crs(england))
inset_box <- britain
Fig5A_inset <- mapmixture(
  admixture_df = admix_df_inset, coords_df = coords_df_inset,
  basemap = rnaturalearthhires::countries10[, "geometry"],
  boundary = c(xmin = -11.50, xmax = 2.50, ymin = 39.00, ymax = 55.10),
  # crs = st_crs(england),
  crs = 4326,
  pie_size = 0.9, cluster_cols = colours2,
  scalebar_position = "br", arrow = F, scalebar_size = 1.9,
  pie_border = 0.10,
)+
  annotate("shadowtext", x = insettextX, y = insettextY, label = insettextLab, size = 3, colour = "black", bg.color = "white")+
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
  inset_element(Fig5A_inset, left = 0.6, bottom = 0.1, right = 1.15, top = 0.8, align_to = "full")

# Pcalcarem K9 structure plot
Fig5B <- structure_plot(
  admixture_df = admix_df$K9, cluster_cols = colours2,
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
ggsave("Figure_05.png", device = png, type = "cairo", width = 12, height = 8, units = "in", dpi = 600)
ggsave("Figure_06.pdf", width = 12, height = 8, units = "in")

# ----------------- #
# Lithothamnion corallioides: PopCluster ####
# ----------------- #

# Read in statistics file
file <- "outputs/Lcorallioides_PopCluster_MedScaling_Min3/"
file_K <- "Lcorallioides_PopCluster_MedScaling_Min3.K"
stats <- read.table(text = readLines(paste0(file, file_K), n = 10), header = TRUE)

# Shorten best run column
stats$BestRun <- str_extract(stats$BestRun, "(Lcorallioides_PopCluster_MedScaling_Min3_).*$")
stats

# Plot stats
plot(x = stats$K, y = as.numeric(stats$DLK2))
plot(x = stats$K, y = as.numeric(stats$FST.FIS))

# Read in admixture proportions: best run for each K
admix_df <- lapply(stats$BestRun, function(i) readPopCluster(file, i, lcor_genind))
names(admix_df) <- paste0("K", stats$K)

# Order site labels
site_order <- c("Mil1","Mil2","Hel","Maw","AusI","AusII","Gri","Swa","Wey")

# Colour palette
scales::show_col(brewer.pal(9, "Pastel1"))
# Cluster order: Maw&Hel, Hel, Mil, MilIso, AusMix, Swa&Wey
colours <- c("#FDDAEC","#B3CDE3","#FBB4AE","#FFFFCC","#FED9A6","#CCEBC5","#DECBE4","#F2F2F2","#E5D8BD")

# Plot single structure plot
structure_plot(admix_df$K6, cluster_cols = colours, site_order = site_order, legend = "right")

# Theme to remove all y axis content
yaxis <- theme(
  axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
)

# Plot all structure plots
plts <- lapply(1:length(admix_df[1:8]), function(i) {
    structure_plot(admixture_df = admix_df[1:8][[i]], cluster_cols = colours, site_order = site_order, legend = "none", site_ticks = F, display_site_labels = F, divider_width = 0.5)+ yaxis
  }
)
png("supporting_figures/Lcorallioides_structure_K2_K9.png", height = 12, width = 8, units = "in", res = 600)
x <- wrap_plots(plts, ncol = 1) & theme(plot.margin = margin(l = 15, b = 0, t = 0, r = 5))
wrap_elements(x)+ theme(plot.margin = margin(l = 10, b = 10, r = 10, t = 10))
dev.off()

# Lcorallioides K6 structure plot
Fig6B <- structure_plot(
  admixture_df = admix_df$K6, cluster_cols = colours,
  site_order = site_order, legend = "top",
  site_ticks = F, site_labels_size = 3.5, site_labels_y = -0.15
)+
  theme(
    axis.title.y = element_text(size = 8, hjust = 0.6),
  )+
  guides(fill = guide_legend(nrow = 1))
Fig6B

# Read in L. corallioides coordinates
coord_lcor <- subset(read.csv("data/site_coordinates3.csv"), Species == "Lcorallioides")

# Reformat for mapmixture
coord_lcor <- dplyr::select(dplyr::select(coord_lcor, !Species), Site, Latitude, Longitude)

# Text labels on map
textX <- c(-4.69, -4.45, -5.40, -2.4, -1.8)
textY <- c(50.09, 50.21, 51.6, 50.48, 50.53)
textLab <- c("Falmouth", "St Austell", "Milford Haven", "Weymouth", "Swanage")

# Admixture map: England and Wales
Fig6A <- mapmixture(
  admixture_df = admix_df$K6, coords_df = coord_lcor,
  basemap = rnaturalearthhires::countries10[, "geometry"],
  boundary = c(xmin = -6.0, xmax = -0.5, ymin = 49.90, ymax = 52.05),
  crs = 4326,
  pie_size = 0.15, cluster_cols = colours[1:6],
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
ggsave("Figure_06.png", device = png, type = "cairo", width = 12, height = 8, units = "in", dpi = 600)
ggsave("Figure_06.pdf", width = 12, height = 8, units = "in")

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
