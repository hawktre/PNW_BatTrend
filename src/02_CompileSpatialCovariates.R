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

# Read in data ------------------------------------------------------------

## nabat covariates already aggregated for Udell et al., 2022
# nabat_covars <- st_read(
#   dsn = here("data/raw/covariates/NABat_grid_covariates/"),
#   layer = "NABat_grid_covariates"
# )

nabat_covars <- read.csv(here("data/raw/covariates/NABat_grid_covariates/Conus_10km_covs.csv"))
deployments <- readRDS(here("data/processed/detections/deployments_to2025.rds"))

## Fix the sample unit id
deployments$sample_unit_id[which(deployments$sample_unit_id == 96444)] <- 95444

## conus_grts key
### Trent and Laura had diff file names.. check that they're the same, or 
### add flexibility to load
possible_conus_files <- c("complete_conus_mastersample_10km_attributed.shp", 
                "conus_mastersample_10km_attributed.shp")

file_to_use <- possible_conus_files[file.exists(
  here::here("data/raw/batgrid",
             possible_conus_files))][1]

conus10k <- read_sf(here("data/raw/batgrid", file_to_use))

##Landfire gap cover
landfire_or <- terra::rast(here(
  "data/raw/covariates/LandFire/LF2022_OR/LC22_EVT_230.tif"
))
landfire_wa <- rast(here(
  "data/raw/covariates/LandFire/LF2022_WA/LC22_EVT_230.tif"
))
landfire_id <- rast(here(
  "data/raw/covariates/LandFire/LF2022_ID/LC22_EVT_230.tif"
))

# subset and plot data ---------------------------------------------------------------
## Define pnw
pnw <- c("Oregon", "Washington", "Idaho")

## subset and plot
nabat_pnw <- nabat_covars %>%
  filter(admin1 %in% pnw | CONUS_10KM %in% unique(deployments$sample_unit_id))

deployments_spat <- deployments %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# Join to get CONUS -------------------------------------------------------
conus_grts_key <- conus10k %>%
  select(CONUS_10KM, GRTS_ID)

## Join with conus and do some formatting
conus_pnw_covars <- nabat_pnw %>%
  left_join(
    as.data.frame(conus_grts_key),
    by = "CONUS_10KM") %>% 
  st_as_sf()  %>%  
  st_transform(crs = 4326)

## rename to make it easier to call
covars <- conus_pnw_covars
plot(covars["CONUS_10KM"], reset = FALSE)

## There is an issue with SU 114921. Should be 114291.
site <- deployments_spat[
  deployments_spat$sample_unit_id == 114921,
] |> 
  st_transform(st_crs(covars))

plot(
  st_geometry(site),
  add = TRUE,
  col = "red",
  pch = 16,
  cex = 1,
  lwd = 0.25
)

## Check again
all(deployments$sample_unit_id %in% unique(covars$CONUS_10KM))

# Get Cliff_Canyon ---------------------------------------------------------------
## Create a single layer
landfire <- terra::merge(landfire_or, landfire_wa, landfire_id, first = T)

## subset all cliffs and canyons categories
lf_vals <- as.data.frame(unique(landfire)) %>%
  filter(str_detect(EVT_NAME, "Cliff") | str_detect(EVT_NAME, "Canyon"))

## Keep cliffs and canyons only
cliff_canyon <- landfire %in% lf_vals$EVT_NAME

## project everything to equal area projection
covars_ea <- covars %>%
  st_transform(st_crs("EPSG:5070"))

cliff_canyon_ea <- cliff_canyon %>%
  terra::project("EPSG:5070", method = "near")

## Convert to terra spatial vector for next operation
covars_ea_sv <- vect(covars_ea)

## resample cliff_canyon (takes a long time)
cliffcanyon_st <- terra::extract(cliff_canyon_ea, covars_ea_sv, mean)

## Join results of extract with sf object
covars <- covars %>%
  rowid_to_column("ID") %>%
  left_join(cliffcanyon_st, by = "ID")


# Write out the results ---------------------------------------------------
# changed to .gpkg to keep column names from reformatting
write_sf(covars, here("data/processed/occurrence/batgrid_covars.gpkg"))
