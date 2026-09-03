## ---------------------------
## Purpose of script: Prepare detection and covariate data for dynamic
##                    occupancy model in Stan (dynocc_autologistic.stan).
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
  "laci", "lano", "myev", "epfu", "myyu", "myth", "myci", 
  "myvo", "anpa", "pahe", "euma", "myca", "mylu", "coto"
)

# Cliff-associated species requiring additional cliff cover covariate
cliff_spp <- c("anpa", "euma", "myci", "pahe")

# Prepare Occurrence Covariates -------------------------------------------
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
dets <- dets %>%
  filter(!is.na(state)) %>% 
  select(-state) %>% 
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
site_ids <- sort(unique(dets$sample_unit_id))
year_ids <- sort(unique(dets$year))
n_sites <- length(site_ids)
n_years <- length(year_ids)

# Keep the base vectors for subsetting inside the species loop
site_idx <- match(dets$sample_unit_id, site_ids)
year_idx <- match(dets$year, year_ids)

# Build Detection Design Matrix -------------------------------------------
det_formula <- ~ clutter_percent + scale(tmin) + scale(dayl) + water_ind
design_matrix_det <- model.matrix(det_formula, data = dets)
vmat_all <- design_matrix_det[, -1, drop = FALSE]

# Build Occupancy Design Matrices -----------------------------------------
nw_grida <- covars_join %>% filter(samp_all == 1) %>% arrange(sample_unit_id)
nw_gridb <- covars_join %>% filter(samp_all == 0) %>% arrange(sample_unit_id)
nw_grid_all <- rbind(nw_grida, nw_gridb)

xmat_all <- nw_grid_all %>%
  column_to_rownames("sample_unit_id") |> 
  rename(log_fc = p_forest) %>%
  select(log_fc, precip, dem_max) %>%
  as.matrix()

xmat_cliff <- nw_grid_all %>%
  column_to_rownames("sample_unit_id") |> 
  rename(log_fc = p_forest, log_cliff = cliff_cover) %>%
  select(log_fc, precip, dem_max, log_cliff) %>%
  as.matrix()

xmata <- xmat_all[which(nw_grid_all$samp_all == 1), ] 
xmatb <- xmat_all[which(nw_grid_all$samp_all == 0), ] 
xmatc <- xmat_cliff[which(nw_grid_all$samp_all == 1), ] 
xmatd <- xmat_cliff[which(nw_grid_all$samp_all == 0), ] 

# Prior Configuration -----------------------------------------------------
PRIOR_TYPE <- "weakly_informative" 
PRIORS_FILE <- "data/raw/Posteriors-for-UpdatingPriors-BatModels.csv" 

# Prior Helpers -----------------------------------------------------------
get_vague_priors <- function(n_occ_covs, n_det_covs, default_sd = sqrt(10)) {
  list(
    a_mu         = 0,
    a_sigma      = default_sd,
    b_mu         = 0,
    b_sigma      = default_sd,
    alpha0_mu    = 0,
    alpha0_sigma = default_sd,
    alphas_mu    = rep(0, n_occ_covs),
    alphas_sigma = rep(default_sd, n_occ_covs),
    beta0_mu     = 0,
    beta0_sigma  = default_sd,
    betas_mu     = rep(0, n_det_covs),
    betas_sigma  = rep(default_sd, n_det_covs)
  )
}

load_informative_priors <- function(file_path, species_name, n_occ_covs, default_sd = sqrt(10)) {
  # We need n_det_covs here; passing it directly is safer, but assuming it exists globally for now.
  # Better practice: determine it from vmat_all globally.
  n_det_covs_local <- ncol(vmat_all)
  priors <- get_vague_priors(n_occ_covs, n_det_covs_local, default_sd = default_sd)
  
  if (!file.exists(file_path)) {
    warning(sprintf("Informative priors file not found. Falling back to weakly informative for '%s'.", species_name))
    return(priors)
  }
  
  df <- tryCatch({
    read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(df) || nrow(df) == 0) return(priors)
  
  spp_upper <- toupper(species_name)
  row_idx <- which(toupper(df$Spp) == spp_upper)
  
  if (length(row_idx) == 0) return(priors)
  
  spp_data <- df[row_idx[1], ]
  
  get_val <- function(col_name) {
    if (!col_name %in% colnames(spp_data)) return(NA)
    val <- spp_data[[col_name]]
    if (is.null(val) || is.na(val) || val == "NA" || val == "") return(NA)
    return(as.numeric(val))
  }
  
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
  
  # 2. 'a' parameter (colonization intercept, mapped from Gamma)
  p_a <- set_prior(get_val("Gamma"), get_val("GammaSD"))
  priors$a_mu    <- p_a$mean
  priors$a_sigma <- p_a$sd
  
  # 3. 'b' parameter (survival effect, mapped from Phi)
  p_b <- set_prior(get_val("Phi"), get_val("PhiSD"))
  priors$b_mu    <- p_b$mean
  priors$b_sigma <- p_b$sd
  
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
  
  # 7. Cliffs (alphas[4])
  if (n_occ_covs >= 4) {
    p_cliffs <- set_prior(get_val("Cliffs"), get_val("CliffsSD"))
    priors$alphas_mu[4]    <- p_cliffs$mean
    priors$alphas_sigma[4] <- p_cliffs$sd
  }
  
  return(priors)
}

# Bundle Data for Stan ----------------------------------------------------
stan_data <- list()

for (spp in possible_bats) {
  xmat_use <- if (spp %in% cliff_spp) xmatc else xmata
  n_occ_covs <- ncol(xmat_use)
  
  # Filter valid observations for this specific species
  dets_spp_raw <- dets[[spp]]
  valid_obs <- !is.na(dets_spp_raw)
  
  dets_vec <- as.integer(dets_spp_raw[valid_obs])
  vmat_spp <- vmat_all[valid_obs, , drop = FALSE]
  
  # Extract valid site and year indices
  site_idx_valid <- site_idx[valid_obs]
  year_idx_valid <- year_idx[valid_obs]
  
  # Build Start and End Index Matrices
  start_idx <- matrix(0, nrow = n_sites, ncol = n_years)
  end_idx <- matrix(0, nrow = n_sites, ncol = n_years)
  
  for (k in seq_along(dets_vec)) {
    s <- site_idx_valid[k]
    y <- year_idx_valid[k]
    
    # First occurrence marks the start index
    if (start_idx[s, y] == 0) start_idx[s, y] <- k
    
    # Continually update the end index
    end_idx[s, y] <- k
  }
  
  # Load prior variables
  priors <- if (PRIOR_TYPE == "informative") {
    load_informative_priors(PRIORS_FILE, spp, n_occ_covs)
  } else {
    get_vague_priors(n_occ_covs, ncol(vmat_spp))
  }
  
  stan_data[[spp]] <- c(
    list(
      n_sites    = n_sites,
      n_years    = n_years,
      n_obs      = length(dets_vec),
      dets       = dets_vec,
      start_idx  = start_idx,
      end_idx    = end_idx,
      n_occ_covs = n_occ_covs,
      xmat       = xmat_use,
      n_det_covs = ncol(vmat_spp),
      vmat       = vmat_spp
    ),
    priors
  )
}

# Save Data Objects & Metadata --------------------------------------------
output_dir <- here("data/processed/results/stan/fits")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

saveRDS(stan_data, here(paste0("data/processed/results/stan/stan_data_", PRIOR_TYPE, ".rds")))

saveRDS(
  list(
    xmat_all = xmat_all, xmat_cliff = xmat_cliff,
    xmata = xmata, xmatb = xmatb, xmatc = xmatc, xmatd = xmatd,
    design_matrix_det = design_matrix_det, vmat_all = vmat_all
  ),
  here("data/processed/results/stan/design_matrices.rds")
)

saveRDS(
  list(
    site_ids = site_ids, year_ids = year_ids,
    cliff_spp = cliff_spp, possible_bats = possible_bats
  ),
  here("data/processed/results/stan/index_keys.rds")
)

# Compile and Fit Model ---------------------------------------------------
stan_model_file <- here("src/dynocc_autologistic.stan")
mod <- cmdstan_model(stan_model_file)

species_to_run <- seq_along(possible_bats)

for (i in species_to_run) {
  spp <- names(stan_data)[i]
  cat("\n--------------------------------------------------\n")
  cat("Fitting Stan model for:", spp, "(", i, "/", length(possible_bats), ")\n")
  cat("Prior type:", PRIOR_TYPE, "\n")
  cat("Start time:", format(Sys.time()), "\n")
  
  fit <- mod$sample(
    data            = stan_data[[spp]],
    chains          = 4,
    parallel_chains = 4,
    iter_warmup     = 1000,
    iter_sampling   = 1000,
    refresh         = 0,
    adapt_delta     = 0.95
  )
  
  save_path <- file.path(output_dir, paste0(spp, "_stanfit_", PRIOR_TYPE, ".rds"))
  fit$save_object(file = save_path)
  
  cat("Finished:", spp, "at", format(Sys.time()), "\n")
  cat("Saved fit to:", save_path, "\n")
  rm(fit)
}

