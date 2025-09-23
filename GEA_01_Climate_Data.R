# ========== #
#
# Maerl WGS Analysis 2025
#
# Extract Climate Data for GEA Analysis
#
# Phymatolithon calcareum
# Lithothamnion corallioides
#
# ========== #

# In RStudio set working directory to the path where this R script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load packages
library(terra)
library(tibble)
library(dplyr)
library(stringr)

# Function to read in world netCDF file and output a cropped tif file of study area
netcdf2tifcrop <- function(netcdf, outpath) {
  path <- netcdf
  nc <- terra::rast(netcdf)
  extent <- terra::ext(c(-15, 10, 35, 60))
  cropped_nc <- terra::crop(nc, extent)
  terra::writeRaster(cropped_nc, filename = paste0(sub("\\.nc$", "", path), "_crop.tif"), overwrite = TRUE)
}

# Run function on all *.nc files
runConversion <- FALSE
if (runConversion) {
  nc_files <- list.files("../BioOracle_NetCDFs/", pattern = ".nc", full.names = TRUE)
  lapply(nc_files, netcdf2tifcrop)
}

# Copy *.tif files to data/ directory
# BioOracle https://www.bio-oracle.org/index.php

# Read in all tif files as rasters
tifs <- list.files(path = "data/environmental_data/", pattern = "\\.tif$", full.names = TRUE)
rasters <- lapply(tifs, terra::rast)

# Rename raster list elements
names(rasters) <- str_extract(tifs, "[a-z|0-9]*_[a-z|0-9]*_[0-9]*_[0-9]*_[a-z]*")
names(rasters) <- str_replace(names(rasters), "mean_", "kdpar_")
names(rasters)

# Combine SpatRasters into a single SpatRaster
allRasters <- rast(rasters)
allRasters[[1]]
plot(allRasters[[1]])

# Create a data frame of longitudes and latitudes
coords <- tribble(
  ~site, ~lon, ~lat,
  "Biz",   -4.986650, 50.136483,
  # "Maw15", -5.030833, 50.155694,
  "Maw", -5.031100, 50.165750,
  "MawC", -5.031100, 50.165750,
  "Ger",   -4.937717, 50.198667,
  "Hel",   -5.100000, 50.097694, # adjusted lon to get value
  "Man",   -5.060200, 50.039417,
  # "Nar",   -4.900000, 50.200000, # adjusted lon to get value
  "Aus",  -4.729850, 50.319417,
  # "AusII", -4.727850, 50.327283,	
  # "Gri",   -4.600000, 50.319150, # adjusted lon to get value
  "Wey",   -2.319622, 50.615993,	
  "Swa",   -1.910432, 50.649958,	
  "Mil1",  -5.130000, 51.690000, # adjusted lon/lat to get value
  "Mil2",  -5.130000, 51.690000, # adjusted lonlat to get value
  # "Tud",   -4.456283, 52.809333,
  "Ons",   -8.915,    42.395,
  "Bor",   -9.020,    42.789,
  "Zar",   -5.500,    54.379, # adjusted lon to get value
  "Mor",   -3.951,    48.711,
  "Tre",   -3.887,    47.795
)
coords

# Convert coordinates to SpatVect object
coords_vect <- vect(coords, geom = c("lon","lat"))
plot(coords_vect, add = TRUE)

# Data frame of climate values
climate_df <- data.frame(
  site = coords_vect$site,
  lat = geom(coords_vect, df = TRUE)[, c("y")],
  lon = geom(coords_vect, df = TRUE)[, c("x")],
  terra::extract(allRasters, coords_vect, ID = FALSE)
)

# Order by temperature and set decimal places
climate_df <- climate_df |> 
  arrange(thetao_baseline_2000_2019_depthmean) |> 
  mutate(across(4:ncol(climate_df), ~ round(.x, digits = 5)))
head(climate_df)

# Test multicollinearity
baseline_datasets <- dplyr::select(climate_df, !contains("lat") & !contains("lon") & !contains("ssp245") & !contains("site"))
head(baseline_datasets)
library(psych)
psych::pairs.panels(baseline_datasets, scale = TRUE)

# Calculate variance inflation factor
library(usdm)
usdm::vif(baseline_datasets)

# Only keep key variables that are not correlated (>0.70)
remove_vars <- c("o2","kdpar")
baseline_datasets_noCor <- dplyr::select(baseline_datasets, !contains(remove_vars))
psych::pairs.panels(baseline_datasets_noCor, scale = TRUE)

# Only keep key variables that are not correlated (>0.70)
remove_vars <- c(remove_vars, "pH")
baseline_datasets_noCor <- dplyr::select(baseline_datasets_noCor, !contains(remove_vars))
psych::pairs.panels(baseline_datasets_noCor, scale = TRUE)
usdm::vif(baseline_datasets_noCor)

# Export climate data that are not correlated
climate_df_noCor <- dplyr::select(climate_df, !contains(remove_vars))

# Export data as CSV file
write.csv(climate_df_noCor, file = "outputs/climate_data.csv", row.names = FALSE)
