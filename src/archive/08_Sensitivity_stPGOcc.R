renv::activate()
## ---------------------------
## Purpose of script: Sensitivity analysis for spatio-temporal dynamic
##                    occupancy model using spOccupancy (stPGOcc) excluding Idaho.
##                    Outputs bundled data objects and fits the stPGOcc model.
##
## Author: Trent VanHawkins
## ---------------------------

# Load Packages -----------------------------------------------------------
library(tidyverse)
library(here)
library(sf)
library(janitor)
library(spOccupancy)

# Load Data ---------------------------------------------------------------
covars <- readRDS(here("data/processed/occurrence/nw_grid.rds")) 
dets   <- readRDS(here("data/processed/detections/nw_nights.rds")) 

# Species to model --------------------------------------------------------
possible_bats <- c(
  "laci",
  "lano",
  "myev",
  "epfu",
  "myyu",
  "myth",
  "myci",
  "myvo",
  "anpa",
  "pahe",
  "euma",
  "myca",
  "mylu",
  "coto"
)

# Cliff-associated species requiring additional cliff cover covariate
cliff_spp <- c("anpa", "euma", "myci", "pahe")

# Prepare Occurrence Covariates -------------------------------------------

## Drop geometry and transform/scale occupancy covariates
covars_join <- covars %>%
  filter(precip >= 0) |> 
  select(-c(lat, long, riverlake)) %>%
  rename(cliff_cover = cliff_canyon) |> 
  mutate(
    p_forest = log(p_forest + 1),
    cliff_cover = log(cliff_cover * 100 + 1),
    across(karst:cliff_cover, ~ scale(.x)[, 1])
  ) %>%
  st_drop_geometry() |> 
  clean_names()

# Prepare Detection Data --------------------------------------------------

## Join detection histories with occurrence covariates
## Filter out Idaho for sensitivity analysis
## Assign replicate IDs within each cell-year combination
## Drop rows missing key occupancy covariates
## IMPORTANT: Sort by site (sample_unit_id), year, and replicate_id for indexing
dets_clean <- dets %>%
  filter(!is.na(state)) %>% # remove records outside study area
  filter(state != "Idaho") %>% # drop Idaho for sensitivity analysis
  select(-state) %>% # drop state column after filtering
  left_join(covars_join, by = "sample_unit_id") %>%
  group_by(sample_unit_id, year) %>%
  mutate(
    replicate_id = as.numeric(factor(
      location_name,
      levels = unique(location_name)
    ))
  ) %>%
  ungroup() %>%
  drop_na(dem_max, p_forest, precip, cliff_cover) %>%
  arrange(sample_unit_id, year, replicate_id)

# Indexing & Spatial Setup ------------------------------------------------

## Dimensions and sorted unique site/year keys
site_ids <- sort(unique(dets_clean$sample_unit_id))
year_ids <- sort(unique(dets_clean$year))

n_sites <- length(site_ids)
n_years <- length(year_ids)

year_scaled <- as.vector(scale(year_ids))
year_mat <- matrix(year_scaled, nrow = n_sites, ncol = n_years, byrow = TRUE)

site_idx <- match(dets_clean$sample_unit_id, site_ids)
year_idx <- match(dets_clean$year, year_ids)
rep_idx  <- dets_clean$replicate_id
n_visits_max <- max(rep_idx)

idx_mat <- cbind(site_idx, year_idx, rep_idx)

## Spatial coordinates for sampled non-Idaho grid cells (projected EPSG:5070 in meters)
coords <- covars %>%
  filter(sample_unit_id %in% site_ids) %>%
  arrange(sample_unit_id) %>%
  st_centroid() %>%
  st_transform(crs = 5070) %>%
  st_coordinates() %>%
  as.matrix()

## Distance matrix for spatial decay prior / initialization
dist_mat <- dist(coords)

# Build Detection & Occurrence Covariates for spOccupancy -----------------

## Sampled cell occupancy covariates
covars_sampled <- covars_join %>%
  filter(sample_unit_id %in% site_ids) %>%
  arrange(sample_unit_id)

occ.covs.std <- list(
  log_fc = covars_sampled$p_forest,
  precip = covars_sampled$precip,
  dem_max = covars_sampled$dem_max,
  year_eff = year_mat
)

occ.covs.cliff <- list(
  log_fc = covars_sampled$p_forest,
  precip = covars_sampled$precip,
  dem_max = covars_sampled$dem_max,
  log_cliff = covars_sampled$cliff_cover,
  year_eff = year_mat
)

## Build 3D detection covariate arrays (n_sites x n_years x n_visits_max)
clutter_num <- as.numeric(as.character(dets_clean$clutter_percent))
tmin_scaled <- as.vector(scale(dets_clean$tmin))
dayl_scaled <- as.vector(scale(dets_clean$dayl))
water_ind   <- dets_clean$water_ind

clutter_arr <- array(NA_real_, dim = c(n_sites, n_years, n_visits_max))
tmin_arr    <- array(NA_real_, dim = c(n_sites, n_years, n_visits_max))
dayl_arr    <- array(NA_real_, dim = c(n_sites, n_years, n_visits_max))
water_arr   <- array(NA_real_, dim = c(n_sites, n_years, n_visits_max))

clutter_arr[idx_mat] <- clutter_num
tmin_arr[idx_mat]    <- tmin_scaled
dayl_arr[idx_mat]    <- dayl_scaled
water_arr[idx_mat]   <- water_ind

det.covs <- list(
  clutter_percent = clutter_arr,
  tmin = tmin_arr,
  dayl = dayl_arr,
  water_ind = water_arr
)

# Model Formulas ----------------------------------------------------------
occ.formula.std   <- ~ log_fc + precip + dem_max + year_eff
occ.formula.cliff <- ~ log_fc + precip + dem_max + log_cliff + year_eff
det.formula       <- ~ clutter_percent + tmin + dayl + water_ind

# Prior Configuration -----------------------------------------------------
# Set PRIOR_TYPE to "weakly_informative" (default) or "informative".
# If "informative", specify the path to your CSV file of priors.
PRIOR_TYPE  <- "weakly_informative" # Options: "weakly_informative", "informative"
PRIORS_FILE <- "data/raw/Posteriors-for-UpdatingPriors-BatModels.csv" 

# Prior Helpers -----------------------------------------------------------
get_vague_priors_spocc <- function(n_occ_params, n_det_params = 5, default_var = 10, dist_mat) {
  list(
    beta.normal = list(mean = rep(0, n_occ_params), var = rep(default_var, n_occ_params)),
    alpha.normal = list(mean = rep(0, n_det_params), var = rep(default_var, n_det_params)),
    sigma.sq.ig = c(2, 1),
    phi.unif = c(3 / max(dist_mat), 3 / min(dist_mat[dist_mat > 0])),
    sigma.sq.t.ig = c(2, 0.5),
    rho.unif = c(-1, 1)
  )
}

load_informative_priors_spocc <- function(file_path, species_name, is_cliff, n_det_params = 5, default_var = 10, dist_mat) {
  n_occ_params <- if (is_cliff) 6 else 5
  priors <- get_vague_priors_spocc(n_occ_params, n_det_params, default_var, dist_mat)
  
  if (!file.exists(file_path)) {
    warning(sprintf("Informative priors file '%s' not found. Falling back to weakly informative priors for species '%s'.", file_path, species_name))
    return(priors)
  }
  
  df <- tryCatch({
    read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)
  }, error = function(e) {
    warning(sprintf("Error reading priors file '%s': %s. Falling back to weakly informative priors.", file_path, e$message))
    NULL
  })
  
  if (is.null(df) || nrow(df) == 0) return(priors)
  
  # Standardize species name to match uppercase Spp column (e.g., "myev" -> "MYEV")
  spp_upper <- toupper(species_name)
  row_idx <- which(toupper(df$Spp) == spp_upper)
  
  if (length(row_idx) == 0) {
    warning(sprintf("No priors found for species '%s' (searched for '%s') in '%s'. Falling back to weakly informative priors.", species_name, spp_upper, file_path))
    return(priors)
  }
  
  spp_data <- df[row_idx[1], ]
  
  # Helper to parse values safely, handling NA representation
  get_val <- function(col_name) {
    if (!col_name %in% colnames(spp_data)) return(NA)
    val <- spp_data[[col_name]]
    if (is.null(val) || is.na(val) || val == "NA" || val == "") return(NA)
    as.numeric(val)
  }
  
  means <- c(get_val("Intercept"), get_val("Forest"), get_val("Precip"), get_val("Elevation"))
  sds   <- c(get_val("InterceptSD"), get_val("ForestSD"), get_val("PrecipSD"), get_val("ElevSD"))
  
  if (is_cliff) {
    means <- c(means, get_val("Cliffs"))
    sds   <- c(sds, get_val("CliffsSD"))
  }
  
  for (k in seq_along(means)) {
    if (!is.na(means[k]) && !is.na(sds[k]) && sds[k] > 0) {
      priors$beta.normal$mean[k] <- means[k]
      priors$beta.normal$var[k]  <- sds[k]^2
    }
  }
  
  return(priors)
}

# Bundle Data for spOccupancy ---------------------------------------------
spocc_data <- list()

for (spp in possible_bats) {
  is_cliff <- spp %in% cliff_spp
  occ.covs <- if (is_cliff) occ.covs.cliff else occ.covs.std
  occ.formula <- if (is_cliff) occ.formula.cliff else occ.formula.std
  
  # Construct 3D detection array y (n_sites x n_years x n_visits_max)
  y_spp <- array(NA_real_, dim = c(n_sites, n_years, n_visits_max))
  y_spp[idx_mat] <- dets_clean[[spp]]
  
  # Initial latent occupancy state z
  z.inits <- apply(y_spp, c(1, 2), function(a) {
    if (all(is.na(a))) return(0)
    as.numeric(sum(a, na.rm = TRUE) > 0)
  })
  
  inits <- list(
    beta = 0,
    alpha = 0,
    sigma.sq = 1,
    phi = 3 / mean(dist_mat),
    z = z.inits,
    rho = 0,
    sigma.sq.t = 0.5
  )
  
  priors <- if (PRIOR_TYPE == "informative") {
    load_informative_priors_spocc(PRIORS_FILE, spp, is_cliff, dist_mat = dist_mat)
  } else {
    n_occ_params <- if (is_cliff) 6 else 5
    get_vague_priors_spocc(n_occ_params, dist_mat = dist_mat)
  }
  
  data.list <- list(
    y = y_spp,
    occ.covs = occ.covs,
    det.covs = det.covs,
    coords = coords
  )
  
  spocc_data[[spp]] <- list(
    data = data.list,
    inits = inits,
    priors = priors,
    occ.formula = occ.formula,
    det.formula = det.formula
  )
}

# Save Data Objects & Metadata --------------------------------------------
output_dir <- here("data/processed/results/stPGOcc/sensitivity/fits")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

## Bundled spOccupancy data
saveRDS(spocc_data, here(paste0("data/processed/results/stPGOcc/spocc_data_", PRIOR_TYPE, "_sensitivity.rds")))

## Site and year index keys
saveRDS(
  list(
    site_ids = site_ids,
    year_ids = year_ids,
    coords = coords,
    cliff_spp = cliff_spp,
    possible_bats = possible_bats
  ),
  here("data/processed/results/stPGOcc/index_keys_sensitivity.rds")
)

# MCMC Sampler & Parallel Settings ----------------------------------------
n.batch      <- 800
batch.length <- 25     # Total samples per chain = 20,000
n.burn       <- 5000
n.thin       <- 5
n.chains     <- 4
n.neighbors  <- 15

# Fit Model for Each Species ----------------------------------------------
# Species to fit: names(spocc_data) for all species, or subset as needed.
species_to_fit <- names(spocc_data)

fit_single_species <- function(spp) {
  cat(sprintf("[%s] Starting stPGOcc fit (Sensitivity - No Idaho) for %s (%s)...\n", format(Sys.time(), "%H:%M:%S"), spp, PRIOR_TYPE))
  spp_obj <- spocc_data[[spp]]
  
  fit <- stPGOcc(
    occ.formula  = spp_obj$occ.formula,
    det.formula  = spp_obj$det.formula,
    data         = spp_obj$data,
    inits        = spp_obj$inits,
    priors       = spp_obj$priors,
    cov.model    = "exponential",
    NNGP         = TRUE,
    n.neighbors  = n.neighbors,
    n.batch      = n.batch,
    batch.length = batch.length,
    n.burn       = n.burn,
    n.thin       = n.thin,
    n.chains     = n.chains,
    ar1          = TRUE,
    n.report     = 100
  )
  
  save_path <- file.path(output_dir, paste0(spp, "_stPGOccfit_", PRIOR_TYPE, ".rds"))
  saveRDS(fit, file = save_path)
  
  cat(sprintf("[%s] Finished %s. Saved fit to %s\n", format(Sys.time(), "%H:%M:%S"), spp, save_path))
  rm(fit)
  gc()
  return(save_path)
}

if (N_CORES > 1 && .Platform$OS.type == "unix") {
  cat(sprintf("\nFitting %d species in parallel across %d cores...\n", length(species_to_fit), N_CORES))
  invisible(parallel::mclapply(species_to_fit, fit_single_species, mc.cores = N_CORES))
} else {
  for (spp in species_to_fit) {
    fit_single_species(spp)
  }
}
cat("\nAll stPGOcc sensitivity model fits complete.\n")
