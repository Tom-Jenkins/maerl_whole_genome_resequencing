# =========================== #
#
# Maerl Whole Genome Resequencing Project 2024
#
# Figure 1
#
# =========================== #

# In RStudio set working directory to the path where this R script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
set.seed(123)

# Load packages
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthhires)
library(ggspatial)
library(sf)
library(mapmixture)
library(patchwork)
library(jpeg)
library(grid)
library(ggrepel)
library(ggsflabel)

# Function to create a grob with a circle split into two colors vertically
create_split_circle <- function(x, y, radius = 0.1, color1 = "red", color2 = "blue") {
  angle_seq <- seq(-pi/2, pi/2, length.out = 100)
  
  # Coordinates for the top half-circle
  top_half <- data.frame(
    x = c(x, x + radius * cos(angle_seq), x),
    y = c(y, y + radius * sin(angle_seq), y)
  )
  
  # Coordinates for the bottom half-circle
  bottom_half <- data.frame(
    x = c(x, x + radius * cos(angle_seq + pi), x),
    y = c(y, y + radius * sin(angle_seq + pi), y)
  )
  
  grobTree(
    # Top half-circle
    polygonGrob(
      x = unit(top_half$x, "npc"), 
      y = unit(top_half$y, "npc"), 
      gp = gpar(fill = color1, col = NA)
    ),
    # Bottom half-circle
    polygonGrob(
      x = unit(bottom_half$x, "npc"), 
      y = unit(bottom_half$y, "npc"), 
      gp = gpar(fill = color2, col = NA)
    )
  )
}

#--------------#
# Figure 1A ####
#--------------#

# CRS
CRS = 3035

# Europe basemap
basemap <- ne_countries(scale = "large", continent = "Europe")[, c("admin")]
basemap <- st_transform(basemap, crs = CRS)

# Boundary
bbox <- mapmixture::transform_bbox(c(xmin = -17, xmax = 4, ymin = 40, ymax = 58), CRS)

# Read in coordinates csv as sf object
coords <- read_csv("data/site_coordinates.csv") |>
  st_as_sf(x = _, coords = c("Longitude","Latitude"), crs = 4326) |>
  st_transform(x = _, crs = CRS)
coords

# Convert species to factor
coords$Species <- factor(coords$Species, levels = c("Pcalcareum","Lcorallioides","Mixed"))

# Theme
gg_theme <- theme(
  panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5),
  panel.background = element_rect(fill = "white"),
  panel.grid.major = element_line(colour = "#f0f0f0"),
  plot.title = element_text(size = 10)
)

# Subset coords to remove Cornwall sites
coords_basemap <- coords |>
  dplyr::filter(Site == "Bem" | Site == "Wey" | Site == "Swa" | Site == "Mil" | Site == "Tud" | Site == "Ons" | Site == "Bor" | Site == "Zar" | Site == "Mor" | Site == "Tre")

# Rect boundary box
rect_bbox <- mapmixture::transform_bbox(
  bbox = c(xmin = -5.20, xmax = -4.50, ymin = 49.95, ymax = 50.50),
  CRS = st_crs(basemap)
)

# Plot sampling basemap
map1 <- ggplot()+
  geom_sf(data = basemap)+
  # geom_sf(data = coords_basemap, aes(fill = Species), size = 4, shape = 21, colour = "black")+
  geom_sf(data = coords_basemap, size = 2, colour = "black")+
  scale_fill_manual("Species Identified:",
                    values = c("#f00c93","pink","#dfc5fe"),
                    labels = expression(italic("P. calcareum"), italic("L. corallioides"), "Both")
                    )+
  annotation_north_arrow(
    data = basemap, location = "tl", which_north = "true",
    height = unit(0.6, "cm"), width = unit(0.6, "cm"),
    pad_y = unit(0.8, "cm"),
    style = north_arrow_orienteering(text_size = 5)
  )+
  annotation_scale(
    data = basemap, location = "tl",
    width_hint = 0.2, bar_cols = c("black","white"),
    # height = unit(0.35, "cm"),
    text_cex = 1
  )+
  annotate("rect", xmin = rect_bbox["xmin"], xmax = rect_bbox["xmax"], ymin = rect_bbox["ymin"], ymax = rect_bbox["ymax"], colour = "black", alpha = 0.50)+
  annotate(
    geom = "segment", arrow = arrow(length = unit(0.10, "inches"), type = "closed"),
    x = rect_bbox["xmin"], xend = rect_bbox["xmin"]*0.95,
    y = rect_bbox["ymin"]*1.01, yend = rect_bbox["ymin"]*1.01
  )+
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"]))+
  xlab("Longitude")+
  ylab("Latitude")+
  gg_theme+
  theme(
    legend.position = "top",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 13),
  )+
  guides(fill = guide_legend(override.aes = aes(label = " ", size = 5)))+
  # Labels
  # geom_sf_label_repel(
  #   data = coords_basemap, aes(label = Site),
  #   min.segment.length = 0, force = 5, size = 4,  label.padding = unit(0.08, "cm")
  # )
  geom_sf_label_repel(
    data = coords_basemap, aes(label = Site, fill = Species),
    min.segment.length = 0, force = 5, size = 4,  label.padding = unit(0.08, "cm")
  )
# map1

# Create data frame for half circles
mix_df <- filter(coords_basemap, Species == "Mixed")
half_circles_df <- data.frame(
  Site = mix_df$Site,
  X = st_coordinates(mix_df)[, "X"],
  Y = st_coordinates(mix_df)[, "Y"]
)

# Loop through the data frame and add custom markers
# for (i in 1:nrow(half_circles_df)) {
#   map1 <- map1 + annotation_custom(
#     grob = create_split_circle(half_circles_df$X[i], half_circles_df$Y[i], radius = 5, color1 = "darkred", color2 = "pink"), 
#     xmin = half_circles_df$X[i], xmax = half_circles_df$X[i], 
#     ymin = half_circles_df$Y[i], ymax = half_circles_df$Y[i]
#   )
# }
# map1

# Read in high resolution England coastline geopackage
england <- st_read("data/england_ol_2001.gpkg")

# Subset coords to keep Cornwall sites
coords_inset <- coords |>
  dplyr::filter(Site == "Biz" | Site == "Maw" | Site == "Ger" | Site == "Hel" | Site == "Man" | Site == "Her" | Site == "AusI" | Site == "AusII" | Site == "Gri")

# Inset boundary box
inset_bbox <- mapmixture::transform_bbox(
  bbox = c(xmin = -5.20, xmax = -4.50, ymin = 49.95, ymax = 50.50),
  CRS = st_crs(england)
)

# Plot inset map
inset_map <- ggplot()+
  geom_sf(data = england)+
  # geom_sf(data = coords_inset, aes(fill = Species), size = 4, shape = 21, alpha = 0.95, colour = "black", show.legend = F)+
  geom_sf(data = coords_inset, size = 2, colour = "black")+
  scale_fill_manual(values = c("#f00c93","pink","#dfc5fe"))+
  annotation_scale(
    data = basemap, location = "br",
    width_hint = 0.10, bar_cols = c("black","white"),
    height = unit(0.35, "cm"),
    text_cex = 1
  )+
  coord_sf(xlim = c(inset_bbox["xmin"], inset_bbox["xmax"]), ylim = c(inset_bbox["ymin"], inset_bbox["ymax"]))+
  xlab("Longitude")+
  ylab("Latitude")+
  gg_theme+
  theme(
    panel.border = element_rect(linewidth = 2, colour = "black"),
  )+
  geom_sf_label_repel(
    data = coords_inset, aes(label = Site, fill = Species),
    show.legend = FALSE, point.size = 3,
    min.segment.length = 0, force = 20, size = 4, label.padding = unit(0.08, "cm")
  )
# inset_map

# New theme for inset map
inset_theme <- theme(
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
  plot.background = element_blank(),
  plot.margin = grid::unit(c(0,0,0,0), "mm")
)

# Add inset to main map
Fig1A <- map1+ inset_element(
  inset_map+ inset_theme,
  left = 0.02, bottom = 0.25, right =  0.53, top = 0.90, align_to = "panel"
)
# Fig1A

#--------------#
# Figure 1B ####
#--------------#

# Read in jpeg and convert to ggplot object
maerl_jpeg <- readJPEG("data/maerl_Matt_Slater_Cornwall_Wildlife_Trust_2021_compressed.jpg") |>
  rasterGrob(image = _, interpolate = TRUE)
grid::textGrob("egege")

# Plot jpeg
Fig1B <- ggplot()+
  annotation_custom(maerl_jpeg, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_void()
# Fig1B

#--------------#
# Figure 1C ####
#--------------#

# Load ggtree
library(ggtree)
library(rphylopic)

# Read in species tree
tree <- read.tree("outputs/busco_genes_ASTRAL_spp1.tree")
tree$edge.length[which(is.na(tree$edge.length))] <- 0.1

# Change tip labels
tree$tip.label[1] <- expression(paste(italic("Porphyra umbilicalis"), phantom()^"DNA"))
tree$tip.label[2] <- expression(paste(italic("Chondrus crispus"), phantom()^"DNA"))
tree$tip.label[3] <- expression(paste(italic("Gracilaria domingensis"), phantom()^"DNA"))
tree$tip.label[4] <- expression(paste(italic("Gracilaria gracilis"), phantom()^"DNA"))
tree$tip.label[5] <- expression(paste(italic("Sporolithon durum"), phantom()^"RNA"))
tree$tip.label[6] <- expression(paste(italic(bold("Phymatolithon calcareum")), phantom()^"DNA"))
tree$tip.label[7] <- expression(paste(italic("Lithothamnion proliferum"), phantom()^"RNA"))
tree$tip.label[8] <- expression(paste(italic(bold("Lithothamnion corallioides")), phantom()^"DNA"))
tree$tip.label[9] <- expression(paste(italic("Calliarthron tuberculosum"), phantom()^"RNA"))
tree$tip.label[10] <- expression(paste(italic("Lithophyllum insipidum"), phantom()^"RNA"))
tree$tip.label[11] <- expression(paste(italic("Porolithon onkodes"), phantom()^"DNA"))
tree$tip.label[12] <- expression(paste(italic("Porolithon onkodes"), phantom()^"RNA"))
tree$tip.label

# Build a ggplot with a geom_tree
plt_tree <- ggtree(tree)+
  geom_treescale()+
  geom_tiplab(parse = TRUE)+
  xlim(0, 10)+
  # Corallinophycidae: calcifying red algae
  # geom_cladelabel(node = 17, label = "Corallinophycidae", color = "purple", offset = 2)+
  geom_hilight(node = 17, fill = "purple", extend = 2.5, alpha = 0.05) +
  # Florideophyceae
  geom_cladelabel(node = 14, label = "Florideophyceae", angle = 270, color = "purple", offset = 3, extend = 0.47, barsize = 1.5, offset.text = 0.25, fontsize = 5, hjust = 0.5)+
  # Bangiophyceae
  geom_cladelabel(node = 1, label = "", color = "red3", offset = 3, align = T, extend = 0.35, barsize = 1.5, offset.text = -1.5, fontsize = 4.5)+
  # Phylopics
  add_phylopic(uuid = "795cf765-7b4a-47a3-8bc3-4c4120fe2393", ysize = 0.8, x = 9.0, y = 3.5, alpha = 0.25)+
  add_phylopic(uuid = "7721d472-7b68-411a-baeb-cab565e4bf89", ysize = 0.8, x = 5.1, y = 2, alpha = 0.25)+
  add_phylopic(uuid = "bdb624c7-049d-45b6-ae41-b1b74186dc69", ysize = 0.8, x = 2.7, y = 1, alpha = 0.25)
# plt_tree

#--------------#
# Figure: Composer ####
#--------------#

# Load patchwork
library(patchwork)

# Layout design
layout <- "
  AABB
  AACC
  AACC
"

# Plot layout
plt_list <- list(
  wrap_elements(Fig1A),
  Fig1B,
  plt_tree
)
fig <- wrap_plots(plt_list, design = layout) + plot_annotation(tag_levels = "A")
# fig

# Export figure
# ggsave(plot = fig, filename = "Figure_01.jpg", width = 12, height = 6, units = "in", dpi = 600)

# Export figure with extra annotations of coarse maerl
# tiff(filename = "Figure_01.tif", width = 17, height = 10, units = "in", res = 600)
jpeg(filename = "Figure_01.jpg", width = 17, height = 10, units = "in", res = 900)
fig
grid.text(
  expression(paste("Coarse ", italic("P. calcareum"))),
  x = 0.84, y = 0.92,
  gp = gpar(fontsize = 12, col = "white", fontface = "bold")
)
grid.segments(
  x0 = 0.81, y0 = 0.90, x1 = 0.77, y1 = 0.85,
  gp = gpar(col = "white", lwd = 1, fill = "black"),
  arrow = arrow(type = "closed", angle = 30, length = unit(0.10, "inches"))
)
dev.off()

# BUSCO scores
# busco <- tibble(
#   Species = c("P. calcareum","L. corallioides"),
#   Complete = c(70, 40),
#   Duplicated = c(5, 5),
#   Fragmented = c(10, 30),
#   Missing = c(10, 20)
# )
# busco
# 
# # Transform data frame
# busco_long <- pivot_longer(busco, cols = 2:5, names_to = "BUSCO", values_to = "Number")
# busco_long
# 
# # Stacked barplot
# ggplot(data = busco_long)+ 
#   geom_col(aes(x = Species, y = Number, fill = BUSCO), width = 0.95)+ 
#   scale_fill_manual(values = c("#2CBBEF", "#0099CF", "#F3E600", "#FF343E"))+
#   coord_flip()+
#   theme_void()
