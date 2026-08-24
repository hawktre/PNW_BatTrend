## ---------------------------
## Purpose of script: Collect and format NABat grid covariates
##
## Author: Trent VanHawkins
## ---------------------------

## load up the packages we will need:  (uncomment as required)

library(tidyverse)
library(here)
library(sf)
library(terra)
library(elevatr)
library(exactextractr)
library()

# Read in data ------------------------------------------------------------
covars <- readRDS(here("data/processed/occurrence/nw_grid.rds"))

# 1. Download a DEM for your study area (or load one you already have)
# (Setting z = 9 gives roughly 30m resolution)
dem <- get_elev_raster(covars, z = 9) 
dem <- rast(dem) # convert to terra SpatRaster

# 4. Extract the proportion of cliff habitat per 10km grid
covars$dem_mean <- exact_extract(dem, covars, 'mean')
covars$topo_rough <- exact_extract(dem, covars, 'stdev')


saveRDS(covars, here("data/processed/occurrence/nw_grid.rds"))
