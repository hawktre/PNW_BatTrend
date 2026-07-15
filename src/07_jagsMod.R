## ---------------------------
## Purpose of script: Prepare detection and covariate data for dynamic
##                    occupancy model in JAGS. Outputs a species-level
##                    list of data objects and runs the model for each species.
##
## Author: Trent VanHawkins
## ---------------------------

# Load Packages -----------------------------------------------------------
library(tidyverse)
library(here)
library(sf)
library(janitor)
library(rjags)
library(jagsUI)

# Load Data ---------------------------------------------------------------
covars <- readRDS(here("data/processed/occurrence/nw_grid.rds")) 
dets <- readRDS(here("data/processed/detections/nw_nights.rds")) 

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
  arrange(year, sample_unit_id, replicate_id)

# Build Detection Array ---------------------------------------------------

## Pivot to long format: one row per species-cell-year-replicate
dets_long <- dets %>%
  pivot_longer(
    cols = all_of(possible_bats),
    names_to = "spp",
    values_to = "occ"
  ) %>%
  select(spp, sample_unit_id, year, replicate_id, occ) %>%
  arrange(spp, year, sample_unit_id, replicate_id) %>%
  mutate(across(1:4, as.factor))

## Create detection array: dim(spp, sample_unit_id, replicate, year)
y <- tapply(
  dets_long$occ,
  select(dets_long, spp, sample_unit_id, replicate_id, year),
  identity
)

# Build Detection Design Matrix -------------------------------------------

## Detection formula - change this to try different covariate structures
## clutter_percent dummy coded with category 0 (open habitat) as reference
## tmin and dayl scaled to mean 0, SD 1
det_formula <- ~ clutter_percent + scale(tmin) + scale(dayl) + water_ind

## Build design matrix from detection data
design_matrix_det <- model.matrix(det_formula, data = dets)
n_pcovs <- ncol(design_matrix_det)

## Dimensions
n_sites_total <- length(unique(dets$sample_unit_id))
n_years <- length(unique(dets$year))
n_visits_max <- max(table(dets$sample_unit_id, dets$year))
site_ids <- sort(unique(dets$sample_unit_id))
year_ids <- sort(unique(dets$year))

## Build 4D pmat array: dim(sample_unit_id, replicate, year, n_pcovs)
## Each slice [i, j, t, ] is the covariate vector for site i, visit j, year t
pmat <- array(0, dim = c(n_sites_total, n_visits_max, n_years, n_pcovs))

for (i in seq_along(site_ids)) {
  for (t in seq_along(year_ids)) {
    rows <- which(dets$sample_unit_id == site_ids[i] & dets$year == year_ids[t])
    if (length(rows) > 0) {
      for (k in seq_along(rows)) {
        pmat[i, k, t, ] <- design_matrix_det[rows[k], ]
      }
    }
  }
}

# Build Occupancy Design Matrices -----------------------------------------

## Separate sampled and unsampled cells for model fitting vs prediction
nw_grida <- covars_join %>% filter(samp_all == 1) %>% arrange(sample_unit_id)
nw_gridb <- covars_join %>% filter(samp_all == 0) %>% arrange(sample_unit_id)

## Recombine with sampled cells first for consistent indexing
nw_grid_all <- rbind(nw_grida, nw_gridb)

## Standard occupancy design matrix: log forest cover, precip, max elevation
## Scaled across all cells (sampled + unsampled) for consistent prediction
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

# Build Visit Count Matrix ------------------------------------------------

## Count replicates per cell per year, fill 0 for unvisited cell-years
n_visits_jags <- dets %>%
  arrange(sample_unit_id, year) %>%
  group_by(sample_unit_id, year) %>%
  summarise(n_visits = n(), .groups = "drop") %>%
  pivot_wider(
    id_cols = sample_unit_id,
    names_from = year,
    values_from = n_visits,
    values_fill = list(n_visits = 0)
  ) %>%
  select(sample_unit_id, all_of(as.character(year_ids))) %>%
  select(-sample_unit_id) %>%
  simplify2array()

# Prior Configuration -----------------------------------------------------
# Set PRIOR_TYPE to "weakly_informative" (default) or "informative".
# If "informative", specify the path to your CSV file of priors.
PRIOR_TYPE <- "weakly_informative" # Options: "weakly_informative", "informative"
PRIORS_FILE <- "data/raw/Posteriors-for-UpdatingPriors-BatModels.csv" 

# Prior Helpers -----------------------------------------------------------
get_vague_priors <- function(n_years, n_xcovs, n_pcovs) {
  list(
    mu_phi = rep(0, n_years - 1),
    tau_phi = rep(0.1, n_years - 1),
    mu_gamma = rep(0, n_years - 1),
    tau_gamma = rep(0.1, n_years - 1),
    mu_alpha01 = 0,
    tau_alpha01 = 0.1,
    mu_alphas = rep(0, n_xcovs),
    tau_alphas = rep(0.1, n_xcovs),
    mu_betas = rep(0, n_pcovs),
    tau_betas = rep(0.1, n_pcovs)
  )
}

load_informative_priors <- function(file_path, species_name, n_years, n_xcovs, n_pcovs) {
  # Default to weakly informative
  priors <- get_vague_priors(n_years, n_xcovs, n_pcovs)
  
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
  
  # Helper to set mean and precision
  # SD is converted to precision: precision = 1 / (SD^2)
  set_prior <- function(mean_val, sd_val, default_mean = 0, default_prec = 0.1) {
    if (is.na(mean_val) || is.na(sd_val) || sd_val <= 0) {
      list(mean = default_mean, prec = default_prec)
    } else {
      list(mean = mean_val, prec = 1 / (sd_val^2))
    }
  }
  
  # 1. Intercept (alpha01)
  p_intercept <- set_prior(get_val("Intercept"), get_val("InterceptSD"))
  priors$mu_alpha01 <- p_intercept$mean
  priors$tau_alpha01 <- p_intercept$prec
  
  # 2. Gamma (colonization)
  p_gamma <- set_prior(get_val("Gamma"), get_val("GammaSD"))
  priors$mu_gamma <- rep(p_gamma$mean, n_years - 1)
  priors$tau_gamma <- rep(p_gamma$prec, n_years - 1)
  
  # 3. Phi (persistence)
  p_phi <- set_prior(get_val("Phi"), get_val("PhiSD"))
  priors$mu_phi <- rep(p_phi$mean, n_years - 1)
  priors$tau_phi <- rep(p_phi$prec, n_years - 1)
  
  # 4. Forest (alphas[1])
  p_forest <- set_prior(get_val("Forest"), get_val("ForestSD"))
  priors$mu_alphas[1] <- p_forest$mean
  priors$tau_alphas[1] <- p_forest$prec
  
  # 5. Precip (alphas[2])
  p_precip <- set_prior(get_val("Precip"), get_val("PrecipSD"))
  priors$mu_alphas[2] <- p_precip$mean
  priors$tau_alphas[2] <- p_precip$prec
  
  # 6. Elevation (alphas[3])
  p_elev <- set_prior(get_val("Elevation"), get_val("ElevSD"))
  priors$mu_alphas[3] <- p_elev$mean
  priors$tau_alphas[3] <- p_elev$prec
  
  # 7. Cliffs (alphas[4]) - only if xmat has at least 4 columns (cliff-associated species)
  if (n_xcovs >= 4) {
    p_cliffs <- set_prior(get_val("Cliffs"), get_val("CliffsSD"))
    priors$mu_alphas[4] <- p_cliffs$mean
    priors$tau_alphas[4] <- p_cliffs$prec
  }
  
  return(priors)
}

# Bundle Data for JAGS ----------------------------------------------------

## Build species-level list of data objects
## Cliff-associated species get the extended design matrix with cliff cover
occ_data <- list()

for (i in possible_bats) {
  xmat_use <- if (i %in% cliff_spp) xmatc else xmata
  n_xcovs <- ncol(xmat_use)
  
  # Load prior variables
  priors <- if (PRIOR_TYPE == "informative") {
    load_informative_priors(PRIORS_FILE, i, n_years, n_xcovs, n_pcovs)
  } else {
    get_vague_priors(n_years, n_xcovs, n_pcovs)
  }

  occ_data[[i]] <- c(
    list(
      # Detection data
      dets = y[i, , , ], # detection array: dim(cell, replicate, year)
      pmat = pmat, # detection design matrix: dim(cell, replicate, year, n_pcovs)
      n_pcovs = n_pcovs, # number of detection covariates

      # Occupancy data
      xmat = xmat_use, # occupancy design matrix: dim(cell, n_xcovs)
      n_xcovs = n_xcovs, # number of occupancy covariates

      # Dimensions
      n_sites = n_sites_total, # number of sampled cells
      n_visits = n_visits_jags, # visit count matrix: dim(cell, year)
      n_years = n_years # number of study years
    ),
    priors # Append prior parameter vectors/values
  )
}

# Save Data Objects -------------------------------------------------------

## Bundled JAGS data for reproducibility
saveRDS(occ_data, here("data/processed/results/jags/occ_data.rds"))

## Design matrices for parameter naming in 01_ParameterSummaries.R
saveRDS(
  list(
    xmata = xmata, # standard species, sampled
    xmatb = xmatb, # standard species, unsampled (prediction)
    xmatc = xmatc, # cliff species, sampled
    xmatd = xmatd, # cliff species, unsampled (prediction)
    design_matrix_det = design_matrix_det # detection design matrix with column names
  ),
  here("data/processed/results/jags/design_matrices.rds")
)

## Site and year index keys for matching model indices to real IDs
saveRDS(
  list(
    site_ids = site_ids,
    year_ids = year_ids
  ),
  here("data/processed/results/jags/index_keys.rds")
)

# Helper Function ---------------------------------------------------------

## Returns max detection across replicates, handling all-NA cases
## Used to initialize latent occupancy state z
max2 <- function(x) {
  if (sum(is.na(x)) == length(x)) {
    return(0)
  }
  return(max(x, na.rm = TRUE))
}

# Fit Model for Each Species ----------------------------------------------

for (i in 1:1) {
  tmp <- occ_data[[i]]
  cat("\nFitting model for:", names(occ_data)[i], "\n")
  cat("Start time:", format(Sys.time()), "\n")

  ## Initialize latent occupancy state z as naive occupancy
  ## (1 if detected at least once across replicates in a given cell-year, 0 otherwise)
  inits <- function() {
    list(z = apply(tmp$dets, c(1, 3), max2))
  }

  ## Fit dynamic occupancy model in JAGS
  occ_jags <- jagsUI::jags(
    data = tmp,
    inits = inits,
    parameters.to.save = c("alpha01", "alphas", "betas", "phi", "gamma", "psi"),
    model.file = here("src/occ_model_royle.jags"),
    n.chains = 4,
    n.iter = 15000,
    n.burnin = 2000,
    n.thin = 1,
    parallel = TRUE
  )

  ## Save species-level model output
  saveRDS(
    occ_jags,
    here(paste0(
      "data/processed/results/jags/full/fits/",
      names(occ_data)[i],
      "_jagsfit_",
      PRIOR_TYPE,
      ".rds"
    ))
  )

  cat("Finished:", names(occ_data)[i], "at", format(Sys.time()), "\n")
}
