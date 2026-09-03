renv::activate()
## ---------------------------
## Purpose of script: Evaluate residual spatial autocorrelation (empirical
##                    correlograms of Moran's I across distance bands) for
##                    dynamic occupancy models across all species, survey years,
##                    and sensitivity analyses (informative, weakly informative,
##                    and sensitivity without Idaho).
##
## Author: Trent VanHawkins
## ---------------------------

# Load Packages -----------------------------------------------------------
library(tidyverse)
library(here)
library(sf)
library(cmdstanr)
library(posterior)

# Directories Setup -------------------------------------------------------
fits_dir <- here("data/processed/results/stan/fits")
base_diag_dir <- here("data/processed/results/stan/diagnostics")

prior_types <- c("weakly_informative", "informative", "weakly_informative_sensitivity")

# Create base diagnostics directory, subdirectories for each prior type and species
if (!dir.exists(base_diag_dir)) {
  dir.create(base_diag_dir, recursive = TRUE, showWarnings = FALSE)
}

for (prior in prior_types) {
  prior_dir <- file.path(base_diag_dir, prior)
  if (!dir.exists(prior_dir)) {
    dir.create(prior_dir, recursive = TRUE, showWarnings = FALSE)
  }
}

# Species to model
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

for (spp in possible_bats) {
  spp_dir <- file.path(base_diag_dir, spp)
  if (!dir.exists(spp_dir)) {
    dir.create(spp_dir, recursive = TRUE, showWarnings = FALSE)
  }
}

# Read in Spatial Grids & Distance Matrices --------------------------------
cat("Precomputing spatial distance matrices for full and sensitivity grids...\n")
bat_grid <- readRDS(here("data/processed/occurrence/nw_grid.rds"))
site_keys_full <- readRDS(here("data/processed/results/stan/index_keys.rds"))
site_keys_sens_path <- here("data/processed/results/stan/index_keys_sensitivity.rds")
site_keys_sens <- if (file.exists(site_keys_sens_path)) readRDS(site_keys_sens_path) else NULL

# Precompute distance matrix for full grid (in km)
grid_full <- bat_grid %>%
  filter(sample_unit_id %in% site_keys_full$site_ids) %>%
  arrange(sample_unit_id)
cents_full <- suppressWarnings(st_centroid(grid_full))
site_dists_full <- units::drop_units(st_distance(cents_full)) / 1000
diag(site_dists_full) <- 0

# Precompute distance matrix for sensitivity grid (in km)
site_dists_sens <- if (!is.null(site_keys_sens)) {
  grid_sens <- bat_grid %>%
    filter(sample_unit_id %in% site_keys_sens$site_ids) %>%
    arrange(sample_unit_id)
  cents_sens <- suppressWarnings(st_centroid(grid_sens))
  d_s <- units::drop_units(st_distance(cents_sens)) / 1000
  diag(d_s) <- 0
  d_s
} else {
  site_dists_full
}

# Load stan data list
stan_data_list <- list(
  weakly_informative = if (file.exists(here("data/processed/results/stan/stan_data_weakly_informative.rds"))) {
    readRDS(here("data/processed/results/stan/stan_data_weakly_informative.rds"))
  } else if (file.exists(here("data/processed/results/stan/stan_data.rds"))) {
    readRDS(here("data/processed/results/stan/stan_data.rds"))
  } else {
    NULL
  },
  informative = if (file.exists(here("data/processed/results/stan/stan_data_informative.rds"))) {
    readRDS(here("data/processed/results/stan/stan_data_informative.rds"))
  } else {
    NULL
  },
  weakly_informative_sensitivity = if (file.exists(here("data/processed/results/stan/stan_data_weakly_informative_sensitivity.rds"))) {
    readRDS(here("data/processed/results/stan/stan_data_weakly_informative_sensitivity.rds"))
  } else {
    NULL
  }
)

# Distance Bins Configuration ---------------------------------------------
dist_bins <- seq(10, 150, by = 10)
bin_bounds <- lapply(seq_along(dist_bins), function(d) {
  list(min = if (d == 1) 0 else dist_bins[d - 1], max = dist_bins[d])
})

# Vectorized Moran's I Computation Helpers --------------------------------

# Fully vectorized Moran's I across all posterior draws for occupancy residuals
compute_moran_occ_fast <- function(occ_res_matrix, dist_mat, bin_bounds, dist_bins) {
  S <- nrow(occ_res_matrix)
  N <- ncol(occ_res_matrix)
  Z_c <- occ_res_matrix - rowMeans(occ_res_matrix)
  denom <- rowSums(Z_c^2)
  
  moran_mat <- matrix(NA_real_, nrow = S, ncol = length(dist_bins))
  for (d in seq_along(bin_bounds)) {
    W <- (dist_mat > bin_bounds[[d]]$min) & (dist_mat <= bin_bounds[[d]]$max)
    diag(W) <- FALSE
    storage.mode(W) <- "double"
    w_sum <- sum(W)
    if (w_sum > 0) {
      numer <- rowSums((Z_c %*% W) * Z_c)
      moran_mat[, d] <- ifelse(denom == 0, NA_real_, (N / w_sum) * (numer / denom))
    }
  }
  
  data.frame(
    distance = rep(dist_bins, each = S),
    draw     = rep(1:S, times = length(dist_bins)),
    morans_i = as.vector(moran_mat)
  )
}

# Fully vectorized Moran's I across all posterior draws for detection residuals
compute_moran_det_fast <- function(y_mat_yr, p_mat_yr, z_exp_yr, obs_site_yr, dist_mat, bin_bounds, dist_bins) {
  S <- nrow(p_mat_yr)
  K <- ncol(p_mat_yr)
  
  V <- (z_exp_yr == 1) * 1.0
  n_vec <- rowSums(V)
  
  raw_diff <- (y_mat_yr - p_mat_yr) * V
  mean_diff <- ifelse(n_vec > 0, rowSums(raw_diff) / n_vec, 0)
  
  R_c <- (raw_diff - mean_diff * V)
  denom <- rowSums(R_c^2)
  
  D_obs <- dist_mat[obs_site_yr, obs_site_yr, drop = FALSE]
  moran_mat <- matrix(NA_real_, nrow = S, ncol = length(dist_bins))
  
  for (d in seq_along(bin_bounds)) {
    W <- (D_obs > bin_bounds[[d]]$min) & (D_obs <= bin_bounds[[d]]$max)
    diag(W) <- FALSE
    storage.mode(W) <- "double"
    
    w_sums <- rowSums((V %*% W) * V)
    numer <- rowSums((R_c %*% W) * R_c)
    
    valid <- (n_vec >= 2) & (w_sums > 0) & (denom > 0)
    moran_mat[valid, d] <- (n_vec[valid] / w_sums[valid]) * (numer[valid] / denom[valid])
  }
  
  data.frame(
    distance = rep(dist_bins, each = S),
    draw     = rep(1:S, times = length(dist_bins)),
    morans_i = as.vector(moran_mat)
  )
}

# Storage for Consolidated Results ----------------------------------------
all_moran_draws_list <- list()
all_moran_summary_list <- list()

# Loop Over Species and Sensitivity Analyses / Priors ---------------------
for (spp in possible_bats) {
  for (prior in prior_types) {
    possible_filenames <- c(
      paste0(spp, "_stanfit_", prior, ".rds"),
      paste0(spp, "_stanfit_", prior, "_spatial.rds")
    )
    
    fit_path <- NULL
    for (fn in possible_filenames) {
      target_path <- file.path(fits_dir, fn)
      if (file.exists(target_path)) {
        fit_path <- target_path
        break
      }
    }
    
    if (is.null(fit_path)) {
      cat(sprintf("Skipping spatial diagnostics for %s (%s): file not found in %s\n", spp, prior, fits_dir))
      next
    }
    
    cat(sprintf("\n=== Processing Spatial Diagnostics: %s (%s) ===\n", spp, prior))
    
    spp_stan_data <- stan_data_list[[prior]][[spp]]
    if (is.null(spp_stan_data)) {
      warning(sprintf("Stan data not found for %s (%s)", spp, prior))
      next
    }
    
    fit <- readRDS(fit_path)
    
    tryCatch({
      # Extract posterior quantities
      mod_rvars <- as_draws_rvars(fit$draws(c("z_sim", "psi", "p")))
      z_array <- draws_of(mod_rvars$z_sim)
      psi_array <- draws_of(mod_rvars$psi)
      p_array <- draws_of(mod_rvars$p)
      
      occ_res <- z_array - psi_array
      
      n_draws <- dim(z_array)[1]
      n_sites <- dim(z_array)[2]
      n_years <- dim(z_array)[3]
      n_obs   <- length(spp_stan_data$dets)
      
      # Subsample 100 posterior draws
      set.seed(42)
      draw_sample <- sample(1:n_draws, min(100, n_draws))
      S <- length(draw_sample)
      
      z_sub <- z_array[draw_sample, , , drop = FALSE]
      occ_res_sub <- occ_res[draw_sample, , , drop = FALSE]
      p_sub <- p_array[draw_sample, , drop = FALSE]
      
      # Map observations to sites and years
      obs_site <- integer(n_obs)
      obs_year <- integer(n_obs)
      for (s in 1:n_sites) {
        for (yr in 1:n_years) {
          st <- spp_stan_data$start_idx[s, yr]
          et <- spp_stan_data$end_idx[s, yr]
          if (st > 0 && et >= st) {
            obs_site[st:et] <- s
            obs_year[st:et] <- yr
          }
        }
      }
      
      # Vectorized expansion of z to observation level
      z_expanded <- matrix(NA_integer_, nrow = S, ncol = n_obs)
      for (k in 1:n_obs) {
        z_expanded[, k] <- z_sub[, obs_site[k], obs_year[k]]
      }
      y_matrix <- matrix(spp_stan_data$dets, nrow = S, ncol = n_obs, byrow = TRUE)
      
      # Select appropriate distance matrix
      site_dists <- if (prior == "weakly_informative_sensitivity") site_dists_sens else site_dists_full
      
      # Compute Moran's I across years for Occupancy and Detection
      plot_data_list <- vector("list", n_years)
      for (yr in 1:n_years) {
        occ_res_yr <- occ_res_sub[, , yr]
        occ_df <- compute_moran_occ_fast(occ_res_yr, site_dists, bin_bounds, dist_bins)
        occ_df$component <- "Occupancy"
        occ_df$year <- yr + 2015
        
        yr_cols <- which(obs_year == yr)
        if (length(yr_cols) > 0) {
          det_df <- compute_moran_det_fast(
            y_matrix[, yr_cols, drop = FALSE],
            p_sub[, yr_cols, drop = FALSE],
            z_expanded[, yr_cols, drop = FALSE],
            obs_site[yr_cols],
            site_dists,
            bin_bounds,
            dist_bins
          )
          det_df$component <- "Detection"
          det_df$year <- yr + 2015
          plot_data_list[[yr]] <- rbind(occ_df, det_df)
        } else {
          plot_data_list[[yr]] <- occ_df
        }
      }
      
      plot_df <- do.call(rbind, plot_data_list) %>%
        mutate(species = spp, prior_type = prior)
      
      summary_plot_data <- plot_df %>%
        group_by(species, prior_type, component, year, distance) %>%
        summarise(
          mean_morans_i = mean(morans_i, na.rm = TRUE),
          sd_morans_i   = sd(morans_i, na.rm = TRUE),
          q2.5          = quantile(morans_i, 0.025, na.rm = TRUE, names = FALSE),
          q25           = quantile(morans_i, 0.25, na.rm = TRUE, names = FALSE),
          q50           = quantile(morans_i, 0.50, na.rm = TRUE, names = FALSE),
          q75           = quantile(morans_i, 0.75, na.rm = TRUE, names = FALSE),
          q97.5         = quantile(morans_i, 0.975, na.rm = TRUE, names = FALSE),
          .groups       = "drop"
        )
      
      key_id <- paste(spp, prior, sep = "_")
      all_moran_draws_list[[key_id]] <- plot_df
      all_moran_summary_list[[key_id]] <- summary_plot_data
      
      # Prior label for plot titles
      prior_label <- case_when(
        prior == "informative" ~ "Informative",
        prior == "weakly_informative" ~ "Weakly Informative",
        TRUE ~ "Weakly Informative (No Idaho)"
      )
      
      # Generate Empirical Correlogram Plot
      p <- ggplot(plot_df, aes(x = distance, y = morans_i)) +
        geom_line(aes(group = draw), alpha = 0.08, color = "#2c3e50") +
        geom_line(data = summary_plot_data, aes(x = distance, y = mean_morans_i, group = 1), color = "#e74c3c", linewidth = 0.8) +
        geom_point(data = summary_plot_data, aes(x = distance, y = mean_morans_i), color = "#e74c3c", size = 1.2) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
        facet_grid(component ~ year, scales = "free_y") +
        theme_bw() +
        labs(
          x = "Distance (km)",
          y = "Moran's I",
          title = paste0("Empirical Correlogram of Occupancy & Detection Residuals - ", toupper(spp), " (", prior_label, ")"),
          subtitle = "Spaghetti lines show 100 posterior draws; red line shows posterior mean"
        )
      
      # Save plot to prior folder and species folder
      prior_out_dir <- file.path(base_diag_dir, prior)
      spp_out_dir <- file.path(base_diag_dir, spp)
      
      ggsave(file.path(prior_out_dir, paste0("spatial_correlogram_", spp, "_", prior, ".png")), plot = p, width = 14, height = 6, units = "in", dpi = 300)
      ggsave(file.path(spp_out_dir, paste0("spatial_correlogram_", spp, "_", prior, ".png")), plot = p, width = 14, height = 6, units = "in", dpi = 300)
      
      cat(sprintf("  Saved spatial correlogram plots for %s (%s)\n", spp, prior))
    }, error = function(e) {
      warning(sprintf("Spatial diagnostics calculation failed for %s (%s): %s", spp, prior, e$message))
    })
    
    rm(fit)
    gc()
  }
}

# Bind Consolidated Results -----------------------------------------------
all_moran_draws_df <- bind_rows(all_moran_draws_list)
all_moran_summary_df <- bind_rows(all_moran_summary_list)

# Save Summary Tables -----------------------------------------------------
if (nrow(all_moran_summary_df) > 0) {
  saveRDS(all_moran_draws_df, file.path(base_diag_dir, "spatial_correlograms_draws.rds"))
  saveRDS(all_moran_summary_df, file.path(base_diag_dir, "spatial_correlograms_summary.rds"))
  write_csv(all_moran_summary_df, file.path(base_diag_dir, "spatial_correlograms_summary.csv"))
  cat("\nSpatial correlogram summaries saved to:", base_diag_dir, "\n")
}

# Multi-Species Summary Figures -------------------------------------------
if (nrow(all_moran_summary_df) > 0) {
  cat("Generating multi-species spatial correlogram summary plots...\n")
  
  # 1. Multi-species occupancy correlograms across years (averaged over years per species)
  occ_spp_avg <- all_moran_summary_df %>%
    filter(component == "Occupancy") %>%
    group_by(species, prior_type, distance) %>%
    summarise(
      mean_I = mean(mean_morans_i, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      species_upper = toupper(species),
      prior_label = case_when(
        prior_type == "informative" ~ "Informative",
        prior_type == "weakly_informative" ~ "Weakly Informative",
        TRUE ~ "Weakly Informative (No Idaho)"
      )
    )
  
  for (p_type in prior_types) {
    p_occ_sub <- occ_spp_avg %>% filter(prior_type == p_type)
    if (nrow(p_occ_sub) > 0) {
      p_label <- unique(p_occ_sub$prior_label)
      p_multi <- ggplot(p_occ_sub, aes(x = distance, y = mean_I, group = species_upper, color = species_upper)) +
        geom_line(linewidth = 0.9, alpha = 0.85) +
        geom_point(size = 1.8) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
        theme_bw() +
        labs(
          title = paste0("Mean Residual Spatial Autocorrelation Across Species (", p_label, ")"),
          x = "Distance (km)",
          y = "Posterior Mean Moran's I (Occupancy)",
          color = "Species"
        )
      
      prior_out_dir <- file.path(base_diag_dir, p_type)
      ggsave(file.path(prior_out_dir, paste0("moran_correlogram_", p_type, ".png")), plot = p_multi, width = 10, height = 6, units = "in", dpi = 300)
      ggsave(file.path(base_diag_dir, paste0("moran_correlogram_", p_type, ".png")), plot = p_multi, width = 10, height = 6, units = "in", dpi = 300)
    }
  }
}

cat("\nAll spatial diagnostics and figures completed successfully.\n")
