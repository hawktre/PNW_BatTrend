## ---------------------------
## Purpose of script: Fit single-species, multi-season occupancy model  WITH spatial random effect
##
## Author: Trent VanHawkins
## ---------------------------
# Load Packages -----------------------------------------------------------
library(tidyverse)
library(here)
library(sf)
library(spOccupancy)

# Load Data ---------------------------------------------------------------
covars <- readRDS(here("data/processed/occurrence/nw_grid.rds")) %>%
  rename(cliff_cover = cliff_canyon)

dets <- readRDS(here("data/processed/detections/nw_nights.rds")) %>%
  rename(tmin = tmin_degrees_c, vp = vp_pa, dayl = dayl_s) %>%
  left_join(
    select(st_drop_geometry(covars), sample_unit_id, state),
    by = "sample_unit_id"
  ) %>%
  rename(state = state.y) %>%
  select(-state.x)

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

cliff_spp <- c("anpa", "euma", "myci", "pahe")

# Prepare Covariates ------------------------------------------------------
covars_scaled <- covars %>%
  st_drop_geometry() %>%
  mutate(
    log_fc = log(p_forest + 1),
    log_cliff = log(cliff_cover * 100 + 1)
  ) %>%
  mutate(across(c(log_fc, precip, dem_max, log_cliff), ~ scale(.x)[, 1]))

covars_sampled <- covars_scaled %>%
  filter(samp_all == 1) %>%
  arrange(sample_unit_id)

# Prepare Detection Data --------------------------------------------------
dets_clean <- dets %>%
  filter(!is.na(state)) %>%
  select(-state) %>%
  left_join(covars_scaled, by = "sample_unit_id") %>%
  group_by(sample_unit_id, year) %>%
  mutate(
    replicate_id = as.numeric(factor(
      location_name,
      levels = unique(location_name)
    ))
  ) %>%
  ungroup() %>%
  drop_na(dem_max, log_fc, precip, log_cliff) %>%
  arrange(year, sample_unit_id, replicate_id)

# Dimensions --------------------------------------------------------------
site_ids <- sort(unique(dets_clean$sample_unit_id))
year_ids <- sort(unique(dets_clean$year))
n_sites <- length(site_ids)
n_years <- length(year_ids)
n_visits_max <- max(table(dets_clean$sample_unit_id, dets_clean$year))

# Build Detection Array ---------------------------------------------------
## spOccupancy expects dim(n_sites, n_years, n_visits)

build_det_array <- function(spp, dets, site_ids, year_ids, n_visits_max) {
  n_sites <- length(site_ids)
  n_years <- length(year_ids)
  y <- array(NA, dim = c(n_sites, n_years, n_visits_max))

  for (i in seq_along(site_ids)) {
    for (t in seq_along(year_ids)) {
      rows <- which(
        dets$sample_unit_id == site_ids[i] &
          dets$year == year_ids[t]
      )
      if (length(rows) > 0) {
        for (k in seq_along(rows)) {
          y[i, t, k] <- dets[[spp]][rows[k]]
        }
      }
    }
  }
  y
}

# Build Detection Covariate Arrays ----------------------------------------
## dim(n_sites, n_years, n_visits) for observation-level covariates

build_det_cov_array <- function(
  cov_name,
  dets,
  site_ids,
  year_ids,
  n_visits_max
) {
  n_sites <- length(site_ids)
  n_years <- length(year_ids)
  arr <- array(NA, dim = c(n_sites, n_years, n_visits_max))

  for (i in seq_along(site_ids)) {
    for (t in seq_along(year_ids)) {
      rows <- which(
        dets$sample_unit_id == site_ids[i] &
          dets$year == year_ids[t]
      )
      if (length(rows) > 0) {
        for (k in seq_along(rows)) {
          arr[i, t, k] <- dets[[cov_name]][rows[k]]
        }
      }
    }
  }
  arr
}

det.covs <- list(
  clutter_percent = build_det_cov_array(
    "clutter_percent",
    dets_clean,
    site_ids,
    year_ids,
    n_visits_max
  ),
  tmin = build_det_cov_array(
    "tmin",
    dets_clean,
    site_ids,
    year_ids,
    n_visits_max
  ),
  dayl = build_det_cov_array(
    "dayl",
    dets_clean,
    site_ids,
    year_ids,
    n_visits_max
  ),
  water_ind = build_det_cov_array(
    "water_ind",
    dets_clean,
    site_ids,
    year_ids,
    n_visits_max
  )
)

## Scale tmin and dayl to match original model
det.covs$tmin <- (det.covs$tmin - mean(det.covs$tmin, na.rm = TRUE)) /
  sd(det.covs$tmin, na.rm = TRUE)
det.covs$dayl <- (det.covs$dayl - mean(det.covs$dayl, na.rm = TRUE)) /
  sd(det.covs$dayl, na.rm = TRUE)

# Build Year Index Matrix -------------------------------------------------
## Integer year index for random year intercept grouping variable
years_index_mat <- matrix(
  rep(seq_len(n_years), each = n_sites),
  nrow = n_sites,
  ncol = n_years,
  byrow = FALSE
)

# Build Spatial Coordinates -----------------------------------------------
## Must be in a projected CRS — Albers Equal Area for PNW
coords <- covars %>%
  filter(sample_unit_id %in% site_ids) %>%
  arrange(sample_unit_id) %>%
  st_centroid() %>%
  st_transform(crs = 5070) %>%
  st_coordinates() %>%
  as.matrix()

## Pairwise distances for prior/init on phi
dist_mat <- dist(coords)

# Formulas ----------------------------------------------------------------
## Same as 06_tPGOcc.R — only addition is spatial random effect via stPGOcc
## (1 | years): random year intercept for flexible non-linear trajectories
## Spatial NNGP random effect captures residual spatial autocorrelation
occ.formula.std <- ~ log_fc + precip + dem_max + (1 | years)
occ.formula.cliff <- ~ log_fc + precip + dem_max + log_cliff + (1 | years)
det.formula <- ~ clutter_percent + tmin + dayl + water_ind

# MCMC Settings -----------------------------------------------------------
## stPGOcc needs more iterations than tPGOcc due to spatial parameters
n.batch <- 1000
batch.length <- 25
n.burn <- 10000
n.thin <- 5
n.chains <- 4

# Output Directory --------------------------------------------------------
dir.create(
  here("data/processed/results/stPGOcc/fits/"),
  recursive = TRUE,
  showWarnings = FALSE
)

# Fit Model for Each Species ----------------------------------------------
for (spp in possible_bats) {
  cat("\nFitting stPGOcc model for:", spp, "\n")
  cat("Start time:", format(Sys.time()), "\n")

  ## Build detection array for this species
  y_spp <- build_det_array(spp, dets_clean, site_ids, year_ids, n_visits_max)

  ## Select occupancy covariates and formula based on species type
  if (spp %in% cliff_spp) {
    occ.covs <- list(
      years = years_index_mat,
      log_fc = covars_sampled$log_fc,
      precip = covars_sampled$precip,
      dem_max = covars_sampled$dem_max,
      log_cliff = covars_sampled$log_cliff
    )
    occ.formula <- occ.formula.cliff
  } else {
    occ.covs <- list(
      years = years_index_mat,
      log_fc = covars_sampled$log_fc,
      precip = covars_sampled$precip,
      dem_max = covars_sampled$dem_max
    )
    occ.formula <- occ.formula.std
  }

  ## Bundle data — coords required for spatial model
  data.list <- list(
    y = y_spp,
    occ.covs = occ.covs,
    det.covs = det.covs,
    coords = coords
  )

  ## Initial values
  ## phi initialized to 3 / mean distance = effective range ~ mean inter-site distance
  z.inits <- apply(y_spp, c(1, 2), function(a) {
    as.numeric(sum(a, na.rm = TRUE) > 0)
  })

  inits <- list(
    beta = 0,
    alpha = 0,
    sigma.sq.psi = 1, # variance of random year intercept
    sigma.sq = 1, # spatial variance
    phi = 3 / mean(dist_mat), # spatial decay
    z = z.inits
  )

  ## Priors
  priors <- list(
    beta.normal = list(mean = 0, var = 2.72),
    alpha.normal = list(mean = 0, var = 2.72),
    sigma.sq.psi.ig = list(a = 0.1, b = 0.1),
    sigma.sq.ig = c(2, 1),
    phi.unif = c(3 / max(dist_mat), 3 / min(dist_mat[dist_mat > 0]))
  )

  ## Fit spatial multi-season occupancy model
  fit <- stPGOcc(
    occ.formula = occ.formula,
    det.formula = det.formula,
    data = data.list,
    inits = inits,
    priors = priors,
    cov.model = "exponential",
    NNGP = TRUE,
    n.neighbors = 5,
    n.batch = n.batch,
    batch.length = batch.length,
    n.burn = n.burn,
    n.thin = n.thin,
    n.chains = n.chains,
    ar1 = FALSE,
    verbose = TRUE,
    n.report = 100
  )

  saveRDS(
    fit,
    here(paste0("data/processed/results/stPGOcc/fits/", spp, "_stPGOccfit.rds"))
  )

  cat("Finished:", spp, "at", format(Sys.time()), "\n")
}
