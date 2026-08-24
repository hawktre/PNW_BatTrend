renv::activate()
## ---------------------------
## Purpose of script: Read in Stan autologistic dynamic occupancy model fits
##                    (dynocc_autologistic.stan) across species and prior types,
##                    compute the site-specific autologistic recursive expectation
##                    across the full Northwest grid (nw_grid.rds), and generate
##                    spatial maps of posterior mean occupancy and uncertainty.
##
## Author: Trent VanHawkins
## ---------------------------

# Load Packages -----------------------------------------------------------
library(tidyverse)
library(here)
library(sf)
library(cmdstanr)
library(posterior)
library(viridis)

# Directories Setup -------------------------------------------------------
fits_dir <- here("data/processed/results/stan/fits")
base_maps_dir <- here("data/processed/results/stan/maps")

prior_types <- c("weakly_informative", "informative", "weakly_informative_sensitivity")

# Create base maps directory and subdirectories for each prior type
if (!dir.exists(base_maps_dir)) {
  dir.create(base_maps_dir, recursive = TRUE, showWarnings = FALSE)
}

for (prior in prior_types) {
  prior_dir <- file.path(base_maps_dir, prior)
  if (!dir.exists(prior_dir)) {
    dir.create(prior_dir, recursive = TRUE, showWarnings = FALSE)
  }
}

# Read in Design Matrices, Data, & Grid -----------------------------------
des_matrix_path <- here("data/processed/results/stan/design_matrices.rds")
des_matrix <- if (file.exists(des_matrix_path)) readRDS(des_matrix_path) else NULL

des_matrix_sens_path <- here("data/processed/results/stan/design_matrices_sensitivity.rds")
des_matrix_sens <- if (file.exists(des_matrix_sens_path)) readRDS(des_matrix_sens_path) else NULL

nw_grid <- readRDS(here("data/processed/occurrence/nw_grid.rds"))

# Species & Prior Setup ---------------------------------------------------
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

# Prepare Trimmed Spatial Geometry & Prediction Grid ----------------------
nw_grid_pred <- nw_grid %>% filter(precip >= 0)
nw_grid_trim <- nw_grid_pred %>% select(sample_unit_id, state, geometry)

# Accumulators for all species --------------------------------------------
all_means_list  <- list()
all_widths_list <- list()

# Loop Over Species and Prior Types ---------------------------------------
for (spp in possible_bats) {
  for (prior in prior_types) {
    # Check for fit file in the single fits directory
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
      cat(sprintf("Skipping Stan predictions for %s (%s): fit file not found in %s\n", spp, prior, fits_dir))
      next
    }
    
    cat(sprintf("\n=== Predicting & Mapping Stan Autologistic: %s (%s) ===\n", spp, prior))
    fit <- readRDS(fit_path)
    
    prior_out_dir <- file.path(base_maps_dir, prior)
    
    des_mat_use <- if (prior == "weakly_informative_sensitivity" && !is.null(des_matrix_sens)) {
      des_matrix_sens
    } else {
      des_matrix
    }
    
    ## Extract appropriate occupancy design matrix (without intercept column)
    xmat <- if (spp %in% cliff_spp) {
      if (!is.null(des_mat_use$xmat_cliff)) des_mat_use$xmat_cliff else des_matrix$xmat_cliff
    } else {
      if (!is.null(des_mat_use$xmat_all)) des_mat_use$xmat_all else des_matrix$xmat_all
    }
    
    site_unit_ids <- rownames(xmat)
    J_pred <- nrow(xmat)
    
    ## Extract parameters for autologistic recurrence
    # Parameters from dynocc_autologistic.stan:
    # alpha0: scalar occupancy intercept
    # alphas: vector of length n_occ_covs (spatial environmental slopes)
    # a: vector of length n_years - 1 (time-varying colonization intercept)
    # b: vector of length n_years - 1 (time-varying autologistic persistence bonus)
    alpha0_draws <- fit$draws(variables = "alpha0", format = "matrix") # dim: n_draws x 1
    alphas_draws <- fit$draws(variables = "alphas", format = "matrix") # dim: n_draws x n_occ_covs
    a_draws      <- fit$draws(variables = "a", format = "matrix")      # dim: n_draws x (n_years - 1)
    b_draws      <- fit$draws(variables = "b", format = "matrix")      # dim: n_draws x (n_years - 1)
    
    n_draws <- nrow(alpha0_draws)
    n_years_fit <- ncol(a_draws) + 1 # e.g. 10 years
    
    cat("  Computing autologistic recursive expectation across grid (J =", J_pred, "cells, draws =", n_draws, ")...\n")
    
    # 1. Precompute spatial environmental effects
    # trans_spatial: J_pred x n_draws (xmat * alphas)
    trans_spatial <- xmat %*% t(alphas_draws)
    
    # psi_spatial: J_pred x n_draws (alpha0 + xmat * alphas)
    psi_spatial <- trans_spatial + matrix(as.vector(alpha0_draws), nrow = J_pred, ncol = n_draws, byrow = TRUE)
    
    # 2. Allocate predictions array (n_draws x J_pred x n_years_fit)
    preds <- array(NA_real_, dim = c(n_draws, J_pred, n_years_fit))
    
    # Year 1: psi_{i, 1} = inv_logit(psi_spatial_i)
    preds[, , 1] <- t(plogis(psi_spatial))
    
    # Years 2..T: Auto-logistic recursive expectation with site-specific transition logits:
    # logit_gamma_{i, t-1} = a_{t-1} + trans_spatial_i
    # logit_phi_{i, t-1}   = a_{t-1} + b_{t-1} + trans_spatial_i
    # psi_{i, t}           = psi_{i, t-1} * phi_{i, t-1} + (1 - psi_{i, t-1}) * gamma_{i, t-1}
    for (t in 2:n_years_fit) {
      a_t <- as.vector(a_draws[, t - 1])
      b_t <- as.vector(b_draws[, t - 1])
      
      logit_phi_mat   <- trans_spatial + matrix(a_t + b_t, nrow = J_pred, ncol = n_draws, byrow = TRUE)
      logit_gamma_mat <- trans_spatial + matrix(a_t, nrow = J_pred, ncol = n_draws, byrow = TRUE)
      
      phi_it   <- t(plogis(logit_phi_mat))   # dim: n_draws x J_pred
      gamma_it <- t(plogis(logit_gamma_mat)) # dim: n_draws x J_pred
      
      preds[, , t] <- preds[, , t - 1] * phi_it + (1 - preds[, , t - 1]) * gamma_it
    }
    
    # 3. Compute summary statistics: mean and 95% CI width
    cat("  Summarizing posterior predictions...\n")
    psi_pred_mean  <- apply(preds, c(2, 3), mean)
    psi_pred_width <- apply(preds, c(2, 3), function(x) {
      quantile2(x, 0.975) - quantile2(x, 0.025)
    })
    
    actual_years <- 2016:(2016 + n_years_fit - 1)
    colnames(psi_pred_mean)  <- paste0("year_", actual_years)
    colnames(psi_pred_width) <- paste0("year_", actual_years)
    
    # Format mean predictions
    mean_df <- as_tibble(psi_pred_mean) %>%
      mutate(sample_unit_id = as.numeric(site_unit_ids)) %>%
      pivot_longer(
        cols = starts_with("year_"),
        names_to = "year",
        names_prefix = "year_",
        values_to = "occ_prob"
      ) %>%
      mutate(
        year = as.numeric(year),
        species = spp,
        prior_type = prior
      ) %>%
      left_join(nw_grid_trim, by = "sample_unit_id") %>%
      st_as_sf()
    
    # Format width (uncertainty) predictions
    width_df <- as_tibble(psi_pred_width) %>%
      mutate(sample_unit_id = as.numeric(site_unit_ids)) %>%
      pivot_longer(
        cols = starts_with("year_"),
        names_to = "year",
        names_prefix = "year_",
        values_to = "ci_width"
      ) %>%
      mutate(
        year = as.numeric(year),
        species = spp,
        prior_type = prior
      ) %>%
      left_join(nw_grid_trim, by = "sample_unit_id") %>%
      st_as_sf()
    
    all_means_list[[paste(spp, prior, sep = "_")]]  <- st_drop_geometry(mean_df)
    all_widths_list[[paste(spp, prior, sep = "_")]] <- st_drop_geometry(width_df)
    
    prior_label <- case_when(
      prior == "informative" ~ "Informative",
      prior == "weakly_informative" ~ "Weakly Informative",
      TRUE ~ "Weakly Informative (Sensitivity)"
    )
    
    # 4. Mean Occupancy Map (2016 vs 2025) -----------------------------------
    cat("  Generating mean occupancy map...\n")
    p_mean <- mean_df %>%
      filter(year %in% c(2016, max(actual_years))) %>%
      ggplot() +
      geom_sf(aes(fill = occ_prob), linewidth = 0) +
      facet_wrap(~year) +
      scale_fill_viridis_c(
        option = "D",
        limits = c(0, 1),
        name = expression("Predicted Mean" ~ hat(bar(psi)))
      ) +
      labs(
        title = paste0("Stan Autologistic: Predicted Posterior Mean Occupancy (", toupper(spp), ")"),
        subtitle = paste0("Priors: ", prior_label, " | Pacific Northwest Grid")
      ) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 13),
        strip.text = element_text(face = "bold", size = 11)
      )
    
    # Save into prior-specific directory
    ggsave(
      filename = file.path(prior_out_dir, paste0("map_mean_psi_", spp, "_", prior, ".png")),
      plot = p_mean,
      width = 9, height = 5, units = "in", dpi = 300
    )
    
    # 5. Posterior Uncertainty / 95% CI Width Map (2016 vs 2025) ------------
    cat("  Generating uncertainty map...\n")
    p_width <- width_df %>%
      filter(year %in% c(2016, max(actual_years))) %>%
      ggplot() +
      geom_sf(aes(fill = ci_width), linewidth = 0) +
      facet_wrap(~year) +
      scale_fill_viridis_c(
        option = "magma",
        limits = c(0, 1),
        name = "95% CI Width"
      ) +
      labs(
        title = paste0("Stan Autologistic: Posterior Uncertainty / 95% CI Width (", toupper(spp), ")"),
        subtitle = paste0("Priors: ", prior_label, " | Pacific Northwest Grid"),
        caption = "Width of 95% Credible Interval (q97.5 - q02.5)"
      ) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 13),
        strip.text = element_text(face = "bold", size = 11)
      )
    
    # Save into prior-specific directory
    ggsave(
      filename = file.path(prior_out_dir, paste0("map_width_psi_", spp, "_", prior, ".png")),
      plot = p_width,
      width = 9, height = 5, units = "in", dpi = 300
    )
    
    rm(fit, preds)
    gc()
  }
}

# Save Combined Prediction Summaries --------------------------------------
if (length(all_means_list) > 0) {
  all_means_df  <- bind_rows(all_means_list)
  all_widths_df <- bind_rows(all_widths_list)
  
  saveRDS(all_means_df, file.path(base_maps_dir, "spatial_mean_predictions.rds"))
  saveRDS(all_widths_df, file.path(base_maps_dir, "spatial_width_predictions.rds"))
  
  cat("\nConsolidated spatial predictions saved to:", base_maps_dir, "\n")
} else {
  cat("\nNo spatial predictions were generated (no model fits found).\n")
}

cat("Done generating spatial maps.\n")
