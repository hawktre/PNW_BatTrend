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
library(rjags)
library(jagsUI)

# Load Data ---------------------------------------------------------------
covars <- readRDS(here("data/processed/occurrence/nw_grid.rds")) %>%
  rename(cliff_cover = cliff_canyon)
dets <- readRDS(here("data/processed/detections/nw_nights.rds")) %>%
  rename(tmin = tmin_degrees_c, vp = vp_pa, dayl = dayl_s) %>%
  left_join(select(covars, sample_unit_id, state), by = "sample_unit_id") %>%
  rename(state = state.y) %>%
  select(-state.x)

# Species to model --------------------------------------------------------
# 14 bat species detected in the Pacific Northwest monitoring program
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
## log(p_forest + 1) and log(cliff_cover * 100 + 1) to address right skew
## All continuous covariates scaled to mean 0, SD 1
covars_join <- covars %>%
  select(-c(lat, long, riverlake)) %>%
  mutate(
    p_forest = log(p_forest + 1),
    cliff_cover = log(cliff_cover * 100 + 1),
    across(karst:cliff_cover, ~ scale(.x)[, 1])
  ) %>%
  st_drop_geometry()

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

## Create 4D detection array: dim(spp, sample_unit_id, replicate, year)
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
nw_grida <- covars %>% filter(samp_all == 1) %>% arrange(sample_unit_id)
nw_gridb <- covars %>% filter(samp_all == 0) %>% arrange(sample_unit_id)

## Recombine with sampled cells first for consistent indexing
nw_grid_all <- rbind(nw_grida, nw_gridb)

## Standard occupancy design matrix: log forest cover, precip, max elevation
## No intercept - alpha01 is estimated separately in the JAGS model
## Scaled across all cells (sampled + unsampled) for consistent prediction
xmat_all <- nw_grid_all %>%
  st_drop_geometry() %>%
  mutate(log_fc = log(p_forest + 1)) %>%
  select(log_fc, precip, dem_max) %>%
  scale()

## Cliff species occupancy design matrix: adds log cliff cover
xmat_cliff <- nw_grid_all %>%
  st_drop_geometry() %>%
  mutate(
    log_fc = log(p_forest + 1),
    log_cliff = log(cliff_cover * 100 + 1)
  ) %>%
  select(log_fc, precip, dem_max, log_cliff) %>%
  scale()

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

# Bundle Data for JAGS ----------------------------------------------------

## Build species-level list of data objects
## Cliff-associated species get the extended design matrix with cliff cover
occ_data <- list()

for (i in possible_bats) {
  xmat_use <- if (i %in% cliff_spp) xmatc else xmata

  occ_data[[i]] <- list(
    # Detection data
    dets = y[i, , , ], # detection array: dim(cell, replicate, year)
    pmat = pmat, # detection design matrix: dim(cell, replicate, year, n_pcovs)
    n_pcovs = n_pcovs, # number of detection covariates

    # Occupancy data
    xmat = xmat_use, # occupancy design matrix: dim(cell, n_xcovs)
    n_xcovs = ncol(xmat_use), # number of occupancy covariates

    # Dimensions
    n_sites = n_sites_total, # number of sampled cells
    n_visits = n_visits_jags, # visit count matrix: dim(cell, year)
    n_years = n_years # number of study years
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

for (i in seq_along(occ_data)) {
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
      "_jagsfit.rds"
    ))
  )

  cat("Finished:", names(occ_data)[i], "at", format(Sys.time()), "\n")
}
