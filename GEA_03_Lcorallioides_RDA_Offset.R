# ========== #
#
# Maerl WGS Analysis 2025
#
# Geonotype-Environment Association Analysis
#
# Species:
# Lithothamnion corallioides
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
vcf <- read.vcfR("outputs/Lcorallioides_SNPs_ramets.vcf")

# Export SNP genotypes in LFMM format
# dartR::gl2geno(vcfR::vcfR2genlight(vcf), outfile = "Lcorallioides_SNPs_ramets", outpath = "outputs/")

# Read in genotypes in LFMM format
geno <- data.table::fread("outputs/Lcorallioides_SNPs_ramets.lfmm")

# Convert environmental data.frame to individual format
geno_ind <- tibble(ind = colnames(vcf@gt)[-1], geno)

# Add column for site names and edit names
geno_ind$site <- sub("_.*", "", geno_ind$ind)
geno_ind <- geno_ind[which(!geno_ind$site == "Tud"), ]
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

# Re-test multicollinearity
library(psych)
psych::pairs.panels(env_matrix)

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
rda1 <- rda(genotypes ~ thetao_baseline_2000_2019_depthmean + so_baseline_2000_2019_depthmean + Condition(pca1$li$Axis1 + dbmems$MEM1 + dbmems$MEM2), data = env_matrix)
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
rda_pvals <- rdadapt(rda1, K = 2)

# P-values threshold after Bonferroni correction
thres_env <- alpha/length(rda_pvals$p.values)

# Identify outlier SNPs
rda_outliers_idx <- which(rda_pvals$p.values < thres_env)
length(rda_outliers_idx)
