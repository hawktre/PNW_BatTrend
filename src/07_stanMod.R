
## ---------------------------
## Purpose of script: Prepare detection and covariate data for dynamic
##                    occupancy model in Stan (TV_dynocc_forward.stan).
##                    Outputs a species-level list of data objects and
##                    fits the Stan model for each species.
##
## Author: Trent VanHawkins
## ---------------------------

# Load Packages -----------------------------------------------------------
library(tidyverse)
library(here)
library(sf)
library(janitor)
library(cmdstanr)

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
## Assign replicate IDs within each cell-year combination
## Drop rows missing key occupancy covariates
## Filter to OR, WA, ID only then drop state column
## IMPORTANT: Sort by site (sample_unit_id), year, and replicate_id for Stan indexing
dets <- dets %>%
  filter(!is.na(state)) %>% # remove records outside study area
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

# Indexing Setup ----------------------------------------------------------

## Dimensions and sorted unique site/year keys
site_ids <- sort(unique(dets$sample_unit_id))
year_ids <- sort(unique(dets$year))
n_sites <- length(site_ids)
n_years <- length(year_ids)
n_site_years <- n_sites * n_years

## Create site_yr_idx mapping each observation to its site-year index (1..n_site_years)
site_idx <- match(dets$sample_unit_id, site_ids)
year_idx <- match(dets$year, year_ids)
site_yr_idx <- (site_idx - 1) * n_years + year_idx

# Build Detection Design Matrix -------------------------------------------

## Detection formula
det_formula <- ~ clutter_percent + scale(tmin) + scale(dayl) + water_ind

## Build design matrix from detection data and drop intercept (beta0 is explicit in Stan)
design_matrix_det <- model.matrix(det_formula, data = dets)
vmat_all <- design_matrix_det[, -1, drop = FALSE]
n_det_covs <- ncol(vmat_all)

# Build Occupancy Design Matrices -----------------------------------------

## Separate sampled and unsampled cells for model fitting vs prediction
nw_grida <- covars_join %>% filter(samp_all == 1) %>% arrange(sample_unit_id)
nw_gridb <- covars_join %>% filter(samp_all == 0) %>% arrange(sample_unit_id)

## Recombine with sampled cells first for consistent indexing
nw_grid_all <- rbind(nw_grida, nw_gridb)

## Standard occupancy design matrix: log forest cover, precip, max elevation
xmat_all <- nw_grid_all %>%
  rename(log_fc = p_forest) %>%
  select(log_fc, precip, dem_max) %>%
  as.matrix()

## Cliff species occupancy design matrix: adds log cliff cover
xmat_cliff <- nw_grid_all %>%
  rename(log_fc = p_forest, log_cliff = cliff_cover) %>%
  select(log_fc, precip, dem_max, log_cliff) %>%
  as.matrix()

## Subset to sampled cells for model fitting
xmata <- xmat_all[which(nw_grid_all$samp_all == 1), ] # standard species, sampled
xmatb <- xmat_all[which(nw_grid_all$samp_all == 0), ] # standard species, unsampled (prediction)
xmatc <- xmat_cliff[which(nw_grid_all$samp_all == 1), ] # cliff species, sampled
xmatd <- xmat_cliff[which(nw_grid_all$samp_all == 0), ] # cliff species, unsampled (prediction)

# Prior Configuration -----------------------------------------------------
# Set PRIOR_TYPE to "weakly_informative" (default) or "informative".
# If "informative", specify the path to your CSV file of priors.
PRIOR_TYPE <- "weakly_informative" # Options: "weakly_informative", "informative"
PRIORS_FILE <- "data/raw/Posteriors-for-UpdatingPriors-BatModels.csv" 

# Prior Helpers -----------------------------------------------------------
get_vague_priors <- function(n_occ_covs, default_sd = 2.5) {
  list(
    phi_mu       = 0,
    phi_sigma    = default_sd,
    gamma_mu     = 0,
    gamma_sigma  = default_sd,
    alpha0_mu    = 0,
    alpha0_sigma = default_sd,
    alphas_mu    = rep(0, n_occ_covs),
    alphas_sigma = rep(default_sd, n_occ_covs)
  )
}

load_informative_priors <- function(file_path, species_name, n_occ_covs, default_sd = 3.162278) {
  # Default to weakly informative
  priors <- get_vague_priors(n_occ_covs, default_sd = default_sd)
  
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
    if (!col_name %in% colnames(spp_data)) {
      return(NA)
    }
    val <- spp_data[[col_name]]
    if (is.null(val) || is.na(val) || val == "NA" || val == "") {
      return(NA)
    }
    return(as.numeric(val))
  }
  
  # Helper to set mean and sd
  set_prior <- function(mean_val, sd_val, default_mean = 0, default_sd_val = default_sd) {
    if (is.na(mean_val) || is.na(sd_val) || sd_val <= 0) {
      list(mean = default_mean, sd = default_sd_val)
    } else {
      list(mean = mean_val, sd = sd_val)
    }
  }
  
  # 1. Intercept (alpha0)
  p_intercept <- set_prior(get_val("Intercept"), get_val("InterceptSD"))
  priors$alpha0_mu    <- p_intercept$mean
  priors$alpha0_sigma <- p_intercept$sd
  
  # 2. Gamma (colonization)
  p_gamma <- set_prior(get_val("Gamma"), get_val("GammaSD"))
  priors$gamma_mu    <- p_gamma$mean
  priors$gamma_sigma <- p_gamma$sd
  
  # 3. Phi (persistence)
  p_phi <- set_prior(get_val("Phi"), get_val("PhiSD"))
  priors$phi_mu    <- p_phi$mean
  priors$phi_sigma <- p_phi$sd
  
  # 4. Forest (alphas[1])
  p_forest <- set_prior(get_val("Forest"), get_val("ForestSD"))
  priors$alphas_mu[1]    <- p_forest$mean
  priors$alphas_sigma[1] <- p_forest$sd
  
  # 5. Precip (alphas[2])
  p_precip <- set_prior(get_val("Precip"), get_val("PrecipSD"))
  priors$alphas_mu[2]    <- p_precip$mean
  priors$alphas_sigma[2] <- p_precip$sd
  
  # 6. Elevation (alphas[3])
  p_elev <- set_prior(get_val("Elevation"), get_val("ElevSD"))
  priors$alphas_mu[3]    <- p_elev$mean
  priors$alphas_sigma[3] <- p_elev$sd
  
  # 7. Cliffs (alphas[4]) - only if xmat has at least 4 columns (cliff-associated species)
  if (n_occ_covs >= 4) {
    p_cliffs <- set_prior(get_val("Cliffs"), get_val("CliffsSD"))
    priors$alphas_mu[4]    <- p_cliffs$mean
    priors$alphas_sigma[4] <- p_cliffs$sd
  }
  
  return(priors)
}

# Bundle Data for Stan ----------------------------------------------------

## Build species-level list of data objects for TV_dynocc_forward.stan
stan_data <- list()

for (spp in possible_bats) {
  xmat_use <- if (spp %in% cliff_spp) xmatc else xmata
  n_occ_covs <- ncol(xmat_use)
  
  # Extract species detection history and filter out NA observations if any exist
  dets_spp_raw <- dets[[spp]]
  valid_obs <- !is.na(dets_spp_raw)
  
  dets_vec <- as.integer(dets_spp_raw[valid_obs])
  site_yr_idx_spp <- site_yr_idx[valid_obs]
  vmat_spp <- vmat_all[valid_obs, , drop = FALSE]
  
  # Load prior variables
  priors <- if (PRIOR_TYPE == "informative") {
    load_informative_priors(PRIORS_FILE, spp, n_occ_covs)
  } else {
    get_vague_priors(n_occ_covs)
  }
  
  stan_data[[spp]] <- c(
    list(
      n_sites      = n_sites,
      n_years      = n_years,
      n_obs        = length(dets_vec),
      n_site_years = n_site_years,
      dets         = dets_vec,
      site_yr_idx  = site_yr_idx_spp,
      n_occ_covs   = n_occ_covs,
      xmat         = xmat_use,
      n_det_covs   = ncol(vmat_spp),
      vmat         = vmat_spp
    ),
    priors
  )
}

# Save Data Objects & Metadata --------------------------------------------

output_dir <- here("data/processed/results/stan/full/fits")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

## Bundled Stan data
saveRDS(stan_data, here("data/processed/results/stan/stan_data.rds"))

## Design matrices for parameter naming / prediction downstream
saveRDS(
  list(
    xmata = xmata,
    xmatb = xmatb,
    xmatc = xmatc,
    xmatd = xmatd,
    design_matrix_det = design_matrix_det,
    vmat_all = vmat_all
  ),
  here("data/processed/results/stan/design_matrices.rds")
)

## Site and year index keys for matching model indices to real IDs
saveRDS(
  list(
    site_ids = site_ids,
    year_ids = year_ids
  ),
  here("data/processed/results/stan/index_keys.rds")
)

# Compile Stan Model ------------------------------------------------------

stan_model_file <- here("src/TV_dynocc_forward.stan")
mod <- cmdstan_model(stan_model_file)

# Fit Model for Each Species ----------------------------------------------

# Change index range (e.g. 1:length(possible_bats)) to run all species
for (i in c(1, 2, 6, 12, 13)) {
  spp <- names(stan_data)[i]
  cat("\n--------------------------------------------------\n")
  cat("Fitting Stan model for:", spp, "\n")
  cat("Start time:", format(Sys.time()), "\n")
  
  data_list <- stan_data[[spp]]
  
  fit <- mod$sample(
    data            = data_list,
    chains          = 4,
    parallel_chains = 4,
    iter_warmup     = 1000,
    iter_sampling   = 1000,
    refresh         = 100,
    adapt_delta = 0.95
  )
  
  save_path <- file.path(output_dir, paste0(spp, "_stanfit_", PRIOR_TYPE, ".rds"))
  fit$save_object(file = save_path)
  
  cat("Finished:", spp, "at", format(Sys.time()), "\n")
  cat("Saved fit to:", save_path, "\n")

  rm(fit)
}

