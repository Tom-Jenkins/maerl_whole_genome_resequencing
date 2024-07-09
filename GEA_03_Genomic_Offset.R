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

# Load packages
library(data.table)
library(readr)
library(stringr)
library(vcfR)
library(dplyr)
library(LEA)
library(purrr)
library(vegan)
library(ggvegan)
library(ggplot2)
library(ggtext)
library(rnaturalearthhires)
library(sf)
library(ggsflabel)
library(patchwork)

# Read in genotypes in LFMM format
geno_outlier <- data.table::fread("outputs/pcalcareum_GEA_outlier_loci.lfmm")

# Read in environmental data per site
env_data <- read_csv("./outputs/climate_data.csv", show_col_types = F)
env_data$site <- str_sub(env_data$site, 1, 3)
env_data

# Convert to individual format
vcf <- read.vcfR("./outputs/pcalcareum_GEA_outlier_loci.vcf.gz")
geno_ind <- tibble(ind = colnames(vcf@gt)[-1], geno_outlier)
geno_ind$site <- str_sub(geno_ind$ind, 1, 3)
head(geno_ind)

# Join data
env_data_ind <- left_join(dplyr::select(geno_ind, ind, site), distinct(env_data, site, .keep_all = T), by = "site")
env_data_ind

# Baseline environmental data
baseline <- select(env_data_ind, contains(c("tao_baseline","so_baseline","ph_baseline")))

# SSP119 'Sustainability'
ssp119 <- select(env_data_ind, contains(c("tao_ssp119","so_ssp119","ph_ssp119")))

# SSP245 'Middle of the Road'
ssp245 <- select(env_data_ind, contains(c("tao_ssp245","so_ssp245","ph_ssp245")))

# SSP585 'Fossil-Fueled'
ssp585 <- select(env_data_ind, contains(c("tao_ssp585","so_ssp585","ph_ssp585")))

# Genetic offset SSP119
offset119 <- genetic.offset(geno_outlier, env = baseline, pred.env = ssp119, K = 9, scale = TRUE)
sort(set_names(round(offset119$offset, digit = 2), env_data_ind$ind))

# Genetic offset SSP245
offset245 <- genetic.offset(geno_outlier, env = baseline, pred.env = ssp245, K = 9, scale = TRUE)
sort(set_names(round(offset245$offset, digit = 2), env_data_ind$ind))

# Genetic offset SSP585
offset585 <- genetic.offset(geno_outlier, env = baseline, pred.env = ssp585, K = 9, scale = TRUE)
sort(set_names(round(offset585$offset, digit = 2), env_data_ind$ind))


#--------------#
# Figure 7A ####
#--------------#

# Run RDA with only outlier loci
rda_outlier <- rda(geno_outlier ~ ., data = baseline, scale = T)

# Total variance explained
RsquareAdj(rda_outlier)

# Compute Variance Inflation Factors
vegan::vif.cca(rda_outlier)

# Create data.frame for plotting
plot_df <- tibble(
  ind = geno_ind$ind,
  site = geno_ind$site
)
plot_df <- cbind(plot_df, baseline)

# Change MawC samples from Maw to MawC
plot_df$site[which(str_detect(plot_df$ind, "Maw11C|Maw22C"))] <- "MawC"

# Change site to factor
site_new_order <- c("Zar","Man","Biz","Aus","Gri","Ger","Her","Maw",
                    "MawC","Swa","Wey","Mor","Tre","Bor","Ons")
plot_df$site <- factor(plot_df$site, levels = site_new_order)

# Vector of colours
library(scales)
sample_cols <- c(Aus = "grey60", Biz = "#FF69B4", Bor = "#E17E68",
                 Ger = "#FCCDE5", Gri = "#FAA0A0", Her = "#D9D9D9",
                 Man = "#F3CFC6", Maw = "#FF00FF", MawC = "#FDB462",
                 Mor = "#BEBADA", Ons = "#FFFFB3", Swa = "grey60", 
                 Tre = "#8DD3C7", Wey = "#B3DE69", Zar = "#80B1D3")
scales::show_col(sample_cols)

# Order colours
sample_cols <- sample_cols[levels(plot_df$site)]

# RDA axes proportion of variance explained
# https://stackoverflow.com/questions/62542609/extracting-proportion-of-variance-explained-from-summaryrda-for-axis-labels
RDA1_percent <- round(summary(rda_outlier)$cont$importance[2,"RDA1"]*100, digits = 1)
RDA2_percent <- round(summary(rda_outlier)$cont$importance[2,"RDA2"]*100, digits = 1)

# Scaling to use
scaling = 3

# Prepare sample points
sites_scores <- fortify(rda_outlier, display = "wa", axes = 1:3, scaling = scaling)
sites_scores$location <- plot_df$site
head(sites_scores)

# Extract data used to plot environmental predictor arrows
arrow1 <- ggplot_build(autoplot(rda_outlier, layers = c("sites","biplot"), scaling = scaling))$data[[2]][1,]
arrow2 <- ggplot_build(autoplot(rda_outlier, layers = c("sites","biplot"), scaling = scaling))$data[[2]][2,]
arrow3 <- ggplot_build(autoplot(rda_outlier, layers = c("sites","biplot"), scaling = scaling))$data[[2]][3,]

# RDA
rda_pcal <- ggplot(data = sites_scores, aes(x = RDA1, y = RDA2, fill = location))+
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = 3)+
  geom_vline(xintercept = 0, linewidth = 0.5, linetype = 3)+
  # Arrow 1
  annotate("segment", x = 0, xend = arrow1$xend, y = 0, yend = arrow1$yend,
           arrow = arrow(length = unit(0.15, "inches")), colour = "black", linewidth = 0.5)+
  annotate("text", x = arrow1$xend, y = arrow1$yend+0.15, label = "Temperature", col = "black", size = 5)+
  # Arrow 2
  annotate("segment", x = 0, xend = arrow2$xend, y = 0, yend = arrow2$yend-0.2,
           arrow = arrow(length = unit(0.15, "inches")), colour = "black", linewidth = 0.5)+
  annotate("text", x = arrow2$xend+0.15, y = arrow2$yend-0.05, label = "Salinity", col = "black", size = 5)+
  # Arrow 3
  annotate("segment", x = 0, xend = arrow3$xend, y = 0, yend = arrow3$yend-0.2,
           arrow = arrow(length = unit(0.15, "inches")), colour = "black", linewidth = 0.5)+
  annotate("text", x = arrow3$xend, y = arrow3$yend-0.30, label = "pH", col = "black", size = 5)+
  scale_fill_manual(name = "Location", values = sample_cols)+
  # scale_x_continuous(limits = c(-2.5,1.5))+
  scale_y_continuous(position = "right")+
  # Sample points
  geom_point(shape = 21, colour = "black", size = 4)+
  xlab(paste0("RDA1 (", RDA1_percent, "%)"))+
  ylab(paste0("RDA2 (", RDA2_percent, "%)"))+
  ggtitle("Outlier Redundancy Analysis")+
  # Text
  # annotate("text", x = -2, y = -2, label = "Warmer & Drier", col = "grey70", size = 4.5)+
  # annotate("text", x = -2, y = 1.2, label = "Warmer & Wetter", col = "grey70", size = 4.5)+
  # annotate("text", x = 1, y = 1.2, label = "Cooler & Wetter", col = "grey70", size = 4.5)+
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA, colour = "black"),
    plot.title = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 11),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13),
    legend.position = "left",
  )
rda_pcal
# ggsave("Figure_07.png", height = 8, width = 10, dpi = 600, units = "in")


#--------------#
# Figure 7B ####
#--------------#

# Temperature change
env_data_ind |> 
  select(site, contains(c("tao_")) & contains(c("baseline","ssp245"))) |> 
  mutate(temp_change = thetao_ssp245_2020_2100_depthmean - thetao_baseline_2000_2019_depthmean) |> 
  distinct() |> 
  pull(thetao_baseline_2000_2019_depthmean) |> range()

# Salinity change
env_data_ind |> 
  select(site, contains(c("so_")) & contains(c("baseline","ssp245"))) |> 
  mutate(salinity_change = so_ssp245_2020_2100_depthmean - so_baseline_2000_2019_depthmean) |> 
  distinct() |> 
  pull(salinity_change) |> range()

# pH change
env_data_ind |> 
  select(site, contains(c("ph_")) & contains(c("baseline","ssp245"))) |> 
  mutate(pH_change = ph_ssp245_2020_2100_depthmean - ph_baseline_2000_2018_depthmean) |> 
  distinct() |> 
  pull(pH_change) |> range()

# Basemap
basemap <- rnaturalearthhires::countries10[, c("geometry")]

# Boundary
boundary = c(xmin = -11.50, xmax = 2.50, ymin = 39.00, ymax = 55.10)

# Data frame of genomic offset values
offset_df <- data.frame(
  ind = env_data_ind$ind,
  site = env_data_ind$site,
  lat = env_data_ind$lat,
  lon = env_data_ind$lon,
  offset119 = round(offset119$offset, digit = 2),
  offset245 = round(offset245$offset, digit = 2),
  offset585 = round(offset585$offset, digit = 2)
)
head(offset_df)

# Mean average offsets per site and convert to sf object
offset_sf <- offset_df |> 
  group_by(site, lat, lon) |> 
  summarise(
    offset119_mean = mean(offset119),
    offset245_mean = mean(offset245),
    offset585_mean = mean(offset585)
  ) |> 
  st_as_sf(coords = c("lon","lat"), crs = 4326)
offset_sf

# Theme
gg_theme <- theme(
  panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5),
  panel.background = element_rect(fill = "white"),
  panel.grid.major = element_line(colour = "#f0f0f0"),
  plot.title = element_text(size = 12, hjust = 0.5),
  axis.title = element_text(size = 11),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 10),
  legend.position = "right",
)

# Label coordinates on map
labels_df <- tribble(
  ~site, ~lat, ~lon,
  "Zar", 54.50, -4.50,
  "Wey", 51.20, -1.30,
  "Biz", 50.09, -4.00,
  "Gri", 50.50, -6.50,
  "Ger", 50.50, -7.50,
  "Her", 50.09, -6.50,
  "Maw", 50.09, -7.50,
  "Man", 50.09, -8.50,
  "Mor", 49.00, -6.50,
  "Tre", 48.00, -6.50,
  "Bor", 42.00, -9.80,
  "Ons", 41.50, -9.80,
)
labels_df

# Text colour vector
text_colour <- c("white","white","black","black","black","black","black","black","white","black","black","white")

# Join to main data frame
offset_sf_join <- left_join(offset_sf, labels_df, by = "site")

# Remove sites which have low sample sizes < 3 (Aus and Swa)
offset_sf_join <- filter(offset_sf_join, site != "Aus", site != "Swa")

# Plot Biz on top of Cornwall sites
Biz <- filter(offset_sf_join, site == "Biz")

# Map
Fig7B <- ggplot()+
  geom_sf(data = basemap)+
  geom_sf(data = offset_sf_join, aes(fill = offset245_mean), size = 2.5, shape = 24, colour = "black", stroke = 0.5)+
  geom_sf(data = Biz, aes(fill = offset245_mean), size = 2.5, shape = 24, colour = "black", stroke = 0.5)+
  # geom_sf(data = offset_sf_join, size = 2.5, shape = 24, fill = "black", colour = "white")+
  # geom_label(
  #   data = offset_sf_join,
  #   aes(label = site, x = lon, y = lat, fill = offset245_mean),
  #   size = 4, colour = text_colour, label.padding = unit(0.10, "cm")
  # )+
  geom_sf_label_repel(
    data = offset_sf_join, aes(label = site, fill = offset245_mean),
    show.legend = FALSE, point.size = 3, max.overlaps = 20,
    min.segment.length = 0, force = 5, size = 4, label.padding = unit(0.08, "cm")
  )+
  coord_sf(xlim = c(boundary["xmin"], boundary["xmax"]), ylim = c(boundary["ymin"], boundary["ymax"]))+
  scale_fill_gradient2(
    name = "Offset", low = "#878787", mid = "#fddbc7", high = "red",
    midpoint = median(offset_sf_join$offset245_mean)
  )+
  # scale_fill_viridis_c(option = "turbo")+
  # scale_fill_distiller(name = "Offset", palette = "RdGy")+
  xlab("Longitude")+
  ylab("Latitude")+
  ggtitle("Genomic Offset – 2050 SSP245")+
  gg_theme
# Fig7B

#--------------#
# Figure: Composer ####
#--------------#

# Layout design
plt_list = list(rda_pcal, Fig7B)
fig <- wrap_plots(plt_list, ncol = 2)+ plot_annotation(tag_levels = "A")
# fig

# Export figure
ggsave(plot = fig, filename = "Figure_07.png", width = 12, height = 7, units = "in", dpi = 600)
ggsave(plot = fig, filename = "Figure_07.pdf", width = 12, height = 7, units = "in")
