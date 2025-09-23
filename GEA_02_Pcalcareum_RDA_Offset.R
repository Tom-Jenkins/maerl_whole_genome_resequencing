# ========== #
#
# Maerl WGS Analysis 2025
#
# Geonotype-Environment Association Analysis
#
# Species:
# Phymatolithon calcareum
#
# ========== #

# In RStudio set working directory to the path where this R script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
set.seed(123)

# Load packages
library(LEA)
library(readr)
library(data.table)
library(ggplot2)
library(gridExtra)
library(vcfR)
library(ade4)
library(adespatial)
library(vegan)
library(stringr)
library(mapmixture)
library(dplyr)
library(purrr)
library(R.utils)
library(dartR)
library(ggvegan)
library(randomcoloR)
library(scales)
library(sf)
library(ggsflabel)
library(rnaturalearth)
library(rnaturalearthhires)
library(patchwork)

# Read in environmental data per site (remove coordinates from data frame)
env_data <- select(read_csv("outputs/climate_data.csv", show_col_types = F), !contains(c("lat","lon")))
env_data

# Read in SNP genotypes for all samples including European maerl beds
vcf_europe <- read.vcfR("outputs/Pcalcareum_SNPs.vcf.gz")
# colnames(vcf_europe@gt)

# Export SNP genotypes in LFMM format
# dartR::gl2geno(vcfR::vcfR2genlight(vcf_europe), outfile = "Pcalcareum_SNPs", outpath = "outputs/")

# Read in genotypes in LFMM format
geno_europe <- data.table::fread("outputs/pcalcareum_SNPs.lfmm")

# Convert environmental data.frame to individual format
geno_ind <- tibble(ind = colnames(vcf_europe@gt)[-1], geno_europe)

# Add column for site names and edit names
geno_ind$site <- sub("_.*", "", geno_ind$ind)
geno_ind$site <- as.factor(str_replace(geno_ind$site, "Maw22C|Maw11C", "MawC"))
geno_ind$site <- as.factor(str_replace(geno_ind$site, "Maw15|Maw22", "Maw"))
geno_ind$site <- as.factor(str_replace(geno_ind$site, "Gri|AusII", "Aus"))
geno_ind$site <- as.factor(str_replace(geno_ind$site, "Her|Nar", "Ger"))
geno_ind <- droplevels(geno_ind[which(!geno_ind$site == "Swa"), ])
geno_ind <- dplyr::select(geno_ind, ind, site, everything())
geno_ind

# Join data
env_data_ind <- left_join(dplyr::select(geno_ind, ind, site), env_data, by = "site")

# Baseline environmental data matrix
remove_future <- c("ssp245")
env_matrix <- dplyr::select(env_data_ind, !contains(remove_future) & !contains("ind") & !contains("site"))

# Re-test multicollinearity
library(psych)
psych::pairs.panels(env_matrix)

# Variable differences among locations
range(env_matrix$thetao_baseline_2000_2019_depthmean)
range(env_matrix$so_baseline_2000_2019_depthmean)
range(env_matrix$sws_baseline_2000_2019_depthmean)

# ----------------- #
# Spatial autocorrelation ####
# ----------------- #

# Read in environmental data coordinates
coords <- select(read_csv("outputs/climate_data.csv", show_col_types = F), contains(c("site","lat","lon")))
coords

# Adjust coordinates so that they fit within the transition object
coords_trans <- coords
coords_trans[1,]$lon <- -5.46
coords_trans[2,]$lon <- -4.95
coords_trans[5,]$lat <- 51.5
coords_trans[6,]$lat <- 51.5
coords_trans[7,]$lat <- 50.2
coords_trans[9,]$lat <- 50.0
coords_trans[10,]$lat <- 50.0
coords_trans[11,]$lat <- 50.0
coords_trans[16,]$lat <- 42.710615
coords_trans[16,]$lon <- -9.070688

# Calculate least-cost geographic distances between sites
library(marmap)
coords_ind <- env_data_ind |> 
  left_join(coords_trans) |> 
  select(site, lon, lat)
# bathydata = getNOAA.bathy(lon1 = -15, lon2 = 30, lat1 = 35, lat2 = 65, resolution = 2)
# trans1 <- trans.mat(bathydata, min.depth = -1)
# save(trans1, file = "outputs/transition_object.RData")
load("outputs/transition_object.RData")
lc_dist_ind = lc.dist(trans1, coords_ind[, c("lon","lat")], res = "dist")

# Compute dbMEMs using adespatial
library(adespatial)
dbmems <- dbmem(lc_dist_ind, MEM.autocor = "positive")

# Note: control for spatial autocorrelation by detrending environmental variables prior to modelling

# Fit linear model to regress temperature against dbmems and extract residuals
temp_regress_data <- cbind(temperature = env_matrix$thetao_baseline_2000_2019_depthmean, dbmems)
temp_regress_mod <- lm(temperature ~ ., data = temp_regress_data)
summary(temp_regress_mod)
temp_residuals <- residuals(temp_regress_mod)

# Fit linear model to regress salinity against dbmems and extract residuals
sal_regress_data <- cbind(salinity = env_matrix$so_baseline_2000_2019_depthmean, dbmems)
sal_regress_mod <- lm(salinity ~ ., data = sal_regress_data)
summary(sal_regress_mod)
sal_residuals <- residuals(sal_regress_mod)

# Fit linear model to regress currents against dbmems and extract residuals
sws_regress_data <- cbind(sws = env_matrix$sws_baseline_2000_2019_depthmean, dbmems)
sws_regress_mod <- lm(sws ~ ., data = sws_regress_data)
summary(sws_regress_mod)
sws_residuals <- residuals(sws_regress_mod)

# Replace environmental values with residuals (controlling for spatial autocorrelation)
env_matrix$thetao_baseline_2000_2019_depthmean <- temp_residuals
env_matrix$so_baseline_2000_2019_depthmean <- sal_residuals
env_matrix$sws_baseline_2000_2019_depthmean <- sws_residuals

# ---------- #
# Redundancy analysis (RDA) ####
# ---------- #

# Expected false discovery rate
alpha <- 0.01

# Individual genotypes
genotypes <- select(geno_ind, -ind, -site)

# Run PCA (use PCs as condition for RDA)
pca1 <- dudi.pca(genotypes, scale = TRUE, scannf = FALSE, nf = 5)
# scatter_plot(as.data.frame(pca1$li), axes = c(1,2), group_ids = geno_ind$site, centroids = F, segments = F)
# scatter_plot(as.data.frame(pca1$li), axes = c(2,3), group_ids = geno_ind$site, centroids = F, segments = F)

# Run partial RDA accounting for population structure and isolation-by-distance (IBD)
# https://popgen.nescent.org/2018-03-27_RDA_GEA.html
rda1 <- rda(genotypes ~ . + Condition(pca1$li$Axis1 + pca1$li$Axis2 + dbmems$MEM1 + dbmems$MEM2), data = env_matrix)
RsquareAdj(rda1)
plot(rda1)

# Variance Inflation Factors
vif.cca(rda1)

# Variance explained by each canonical axis
summary(eigenvals(rda1, model = "constrained"))

# Screeplot (same number of axes as the number of predictors in the model)
screeplot(rda1)

# Check each canonical axis for significance (**very long run time ~30 mins)
# signif_axis <- anova.cca(rda1, by = "axis", parallel = 4)
# signif_axis

# Function from Capblancq and Forester 2021 paper which returns q-values from RDA
# https://github.com/Capblancq/RDA-landscape-genomics/blob/main/src/rdadapt.R
# https://doi.org/10.1111/2041-210X.13722
rdadapt <- function(rda, K)
{
  zscores <- rda$CCA$v[, 1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- robust::covRob(resscale, distance = TRUE, na.action = na.omit, estim = "pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5, df=K)
  reschi2test <- pchisq(resmaha/lambda, K, lower.tail=FALSE)
  qval <- qvalue::qvalue(reschi2test)
  q.values_rdadapt <- qval$qvalues
  return(data.frame(p.values = reschi2test, q.values = q.values_rdadapt))
}
rda_pvals <- rdadapt(rda1, K = 3)

# P-values threshold after Bonferroni correction
thres_env <- alpha/length(rda_pvals$p.values)

# Identify outlier SNPs
rda_outliers_idx <- which(rda_pvals$p.values < thres_env)
length(rda_outliers_idx)

# ---------- #
# Outlier RDA ####
# ---------- #

# Subset genotypes for outlier SNPs
geno_outlier_SNPs <- as.data.frame(genotypes)[, rda_outliers_idx]

# Env matrix with different column names
colnames(env_matrix) <- str_extract(colnames(env_matrix), "^[^_]+")

# Run RDA with only outlier SNPs
rda_outlier <- rda(geno_outlier_SNPs ~ ., data = env_matrix)
# rda_outlier <- rda(geno_outlier_SNPs ~ . + dbmems$MEM1, data = env_matrix)
plot(rda_outlier)

# Total variance explained
RsquareAdj(rda_outlier)

# PCA and visualisation
pca2 <- dudi.pca(geno_outlier_SNPs, scale = TRUE, scannf = FALSE, nf = 5)
sample_cols <- distinctColorPalette(n_distinct(geno_ind$site), runTsne = FALSE)
scatter_plot(as.data.frame(pca2$li), axes = c(1,2), group_ids = geno_ind$site, centroids = F, segments = F, colours = sample_cols)
# scatter_plot(as.data.frame(pca2$li), axes = c(2,3), group_ids = geno_ind$site, centroids = F, segments = F)

# Create data.frame for plotting
plot_df <- tibble(
  ind = geno_ind$ind,
  site = geno_ind$site
)
plot_df <- cbind(plot_df, env_matrix)

# Change site to factor
unique(plot_df$site)
site_new_order <- c("Zar","Man","Biz","Aus","Ger","Maw","MawC","Wey","Mor","Tre","Bor","Ons")
plot_df$site <- factor(plot_df$site, levels = site_new_order)

# Vector of colours
sample_cols <- c(Aus = "#FAA0A0", Biz = "#FF69B4", Bor = "#E17E68",
                 Ger = "#D9D9D9", 
                 Man = "#F3CFC6", Maw = "#FF00FF", MawC = "#FDB462",
                 Mor = "#BEBADA", Ons = "#FFFFB3", 
                 Tre = "#8DD3C7", Wey = "#B3DE69", Zar = "#80B1D3")
# sample_cols <- distinctColorPalette(n_distinct(plot_df$site), runTsne = FALSE)
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
arrow4 <- ggplot_build(autoplot(rda_outlier, layers = c("sites","biplot"), scaling = scaling))$data[[2]][4,]
arrow5 <- ggplot_build(autoplot(rda_outlier, layers = c("sites","biplot"), scaling = scaling))$data[[2]][5,]

# RDA
rda_pcal <- ggplot(data = sites_scores, aes(x = RDA1, y = RDA2, fill = location))+
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = 3)+
  geom_vline(xintercept = 0, linewidth = 0.5, linetype = 3)+
  # Arrow 1
  annotate("segment", x = 0, xend = arrow1$xend, y = 0, yend = arrow1$yend,
           arrow = arrow(length = unit(0.15, "inches")), colour = "black", linewidth = 0.5)+
  annotate("text", x = arrow1$xend+0.40, y = arrow1$yend, label = "Salinity", col = "black", size = 5)+
  # Arrow 2
  annotate("segment", x = 0, xend = arrow2$xend, y = 0, yend = arrow2$yend-0.2,
           arrow = arrow(length = unit(0.15, "inches")), colour = "black", linewidth = 0.5)+
  annotate("text", x = arrow2$xend+0.35, y = arrow2$yend-0.35, label = "Sea water velocity", col = "black", size = 5)+
  # Arrow 3
  annotate("segment", x = 0, xend = arrow3$xend, y = 0, yend = arrow3$yend-0.2,
           arrow = arrow(length = unit(0.15, "inches")), colour = "black", linewidth = 0.5)+
  annotate("text", x = arrow3$xend, y = arrow3$yend, label = "Temperature", col = "black", size = 5)+
  # # Arrow 4
  # annotate("segment", x = 0, xend = arrow4$xend, y = 0, yend = arrow4$yend-0.2,
  #          arrow = arrow(length = unit(0.15, "inches")), colour = "black", linewidth = 0.5)+
  # annotate("text", x = arrow4$xend, y = arrow4$yend, label = "Temperature", col = "black", size = 5)+
  # # Arrow 5
  # annotate("segment", x = 0, xend = arrow5$xend, y = 0, yend = arrow5$yend-0.2,
  #          arrow = arrow(length = unit(0.15, "inches")), colour = "black", linewidth = 0.5)+
  # annotate("text", x = arrow5$xend, y = arrow5$yend, label = "Temperature", col = "black", size = 5)+
  scale_fill_manual(
    name = "Site", values = sample_cols)+
  # scale_x_continuous(limits = c(-2.5,1.5))+
  scale_y_continuous(position = "right")+
  # Sample points
  geom_point(shape = 21, colour = "black", size = 4)+
  xlab(paste0("RDA1 (", RDA1_percent, "%)"))+
  ylab(paste0("RDA2 (", RDA2_percent, "%)"))+
  ggtitle("RDA with outlier SNP genotypes")+
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
    axis.title.x = element_text(vjust = -1),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13),
    legend.position = "left",
    legend.key.height = unit(1.5, "lines"),
    legend.spacing.y = unit(0.8, "cm"),
  )
rda_pcal
# ggsave("Figure_07.png", height = 8, width = 10, dpi = 600, units = "in")

# ---------- #
# Genomic offset ####
# ---------- #

# Baseline environmental data
baseline <- select(env_data_ind, contains(c("tao_baseline","so_baseline","sws_baseline")))

# SSP245 'Middle of the Road'
ssp245 <- select(env_data_ind, contains(c("tao_ssp245","so_ssp245","sws_baseline")))

# Genetic offset SSP245
genomic_offset <- genetic.offset(geno_outlier_SNPs, env = baseline, pred.env = ssp245, K = 5, scale = TRUE)
sort(set_names(round(genomic_offset$offset, digit = 3), env_data_ind$site))
sort(set_names(round(genomic_offset$offset, digit = 3), env_data_ind$ind))

# Temperature change
env_data_ind |> 
  select(site, contains(c("tao_")) & contains(c("baseline","ssp245"))) |> 
  mutate(temp_change = thetao_ssp245_2020_2100_depthmean - thetao_baseline_2000_2019_depthmean) |> 
  distinct() |> 
  pull(temp_change) |> range()

# Salinity change
env_data_ind |> 
  select(site, contains(c("so_")) & contains(c("baseline","ssp245"))) |> 
  mutate(salinity_change = so_ssp245_2020_2100_depthmean - so_baseline_2000_2019_depthmean) |> 
  distinct() |> 
  pull(salinity_change) |> range()

# SWS change
env_data_ind |> 
  select(site, contains(c("sws_")) & contains(c("baseline","ssp245"))) |> 
  mutate(sws_change = sws_ssp245_2020_2100_depthmean - sws_baseline_2000_2019_depthmean) |> 
  distinct() |> 
  pull(sws_change) |> range()

# Basemap
basemap <- rnaturalearthhires::countries10[, c("geometry")]

# Boundary
boundary = c(xmin = -11.50, xmax = 2.50, ymin = 39.00, ymax = 55.10)

# Data frame of genomic offset values
offset_df <- data.frame(
  ind = env_data_ind$ind,
  site = env_data_ind$site,
  offset245 = round(genomic_offset$offset, digit = 3)
)
head(offset_df)

# Mean average offsets per site and convert to sf object
offset_sf <- offset_df |> 
  left_join(coords, by = "site") |> 
  group_by(site, lat, lon) |> 
  summarise(offset245_mean = mean(offset245)) |> 
  st_as_sf(coords = c("lon","lat"), crs = 4326) |> 
  distinct(site, .keep_all = TRUE)
offset_sf

# Label coordinates on map
labels_df <- tribble(
  ~site, ~lat, ~lon,
  "Zar", 54.50, -4.50,
  "Wey", 51.20, -1.30,
  "Biz", 50.09, -4.00,
  "Aus", 50.50, -6.50,
  "Ger", 50.50, -7.50,
  # "Nar", 50.09, -6.50,
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
offset_sf_join <- filter(offset_sf_join, site != "Swa")

# Plot Biz on top of Cornwall sites
Biz <- filter(offset_sf_join, site == "Biz")

# Theme
gg_theme <- theme(
  panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5),
  panel.background = element_rect(fill = "white"),
  panel.grid.major = element_line(colour = "#f0f0f0"),
  plot.title = element_text(size = 12, hjust = 0.5),
  axis.title = element_text(size = 11),
  axis.title.x = element_text(vjust = -1),
  legend.title = element_text(size = 14, margin = margin(b = 15)),
  legend.text = element_text(size = 10),
  legend.position = "right",
)

# Map
offset_map <- ggplot()+
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
    min.segment.length = 0, force = 3, size = 4, label.padding = unit(0.08, "cm")
  )+
  coord_sf(xlim = c(boundary["xmin"], boundary["xmax"]), ylim = c(boundary["ymin"], boundary["ymax"]))+
  scale_fill_gradient2(
    name = "Offset",
    low = "#878787",
    mid = "#fddbc7",
    high = "#fc8d59",
    midpoint = mean(offset_sf_join$offset245_mean),
    # breaks = c(0.02,0.03,0.04)
    # midpoint = 0.03
  )+
  # scale_fill_viridis_c(option = "turbo")+
  # scale_fill_distiller(name = "Offset", palette = "RdGy")+
  xlab("Longitude")+
  ylab("Latitude")+
  ggtitle("Genomic offset – 2050 SSP245")+
  gg_theme
offset_map

# Figure layout design
plt_list = list(rda_pcal, offset_map)
fig <- wrap_plots(plt_list, ncol = 2)+ plot_annotation(tag_levels = "A")
# fig

# Export figure
ggsave(plot = fig, filename = "figures/Figure_05.png", width = 12, height = 7, units = "in", dpi = 600)
ggsave(plot = fig, filename = "figures/Figure_05.pdf", width = 12, height = 7, units = "in")
