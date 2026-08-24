renv::activate()
## ---------------------------
## Purpose of script: Read in stPGOcc model fits across species and prior types,
##                    predict occupancy probabilities across the entire regional grid
##                    (nw_grid.rds) for selected years (e.g., 2016 vs 2025), and
##                    generate spatial maps of posterior mean occupancy and uncertainty.
##
## Author: Trent VanHawkins
## ---------------------------

# Load Packages -----------------------------------------------------------
library(tidyverse)
library(here)
library(sf)
library(spOccupancy)
library(janitor)
library(viridis)

# Directories Setup -------------------------------------------------------
fits_dir <- here("data/processed/results/stPGOcc/full/fits")
maps_dir <- here("data/processed/results/stPGOcc/full/maps")

if (!dir.exists(maps_dir)) {
  dir.create(maps_dir, recursive = TRUE, showWarnings = FALSE)
}

# Load Spatial Grid & Metadata --------------------------------------------
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
prior_types <- c("weakly_informative", "informative")

# Year indexing
index_keys_path <- here("data/processed/results/stPGOcc/index_keys.rds")
if (file.exists(index_keys_path)) {
  index_keys <- readRDS(index_keys_path)
  year_ids <- index_keys$year_ids
} else {
  year_ids <- 2016:2025
}

# Prepare Prediction Grid Covariates --------------------------------------
## Filter to valid study grid cells
nw_grid_pred <- nw_grid %>% filter(precip >= 0)

## Transform/scale occupancy covariates matching model fitting
covars_pred <- nw_grid_pred %>%
  select(-c(lat, long, riverlake)) %>%
  rename(cliff_cover = cliff_canyon) %>%
  mutate(
    p_forest = log(p_forest + 1),
    cliff_cover = log(cliff_cover * 100 + 1),
    across(karst:cliff_cover, ~ scale(.x)[, 1])
  ) %>%
  st_drop_geometry() %>%
  clean_names()

## Prediction coordinates in EPSG:5070
coords_pred <- nw_grid_pred %>%
  st_centroid() %>%
  st_transform(crs = 5070) %>%
  st_coordinates() %>%
  as.matrix()

J_pred <- nrow(covars_pred)

# Years to predict (start and end year, or full sequence)
t_cols <- c(1, length(year_ids)) # First (2016) and Last (2025)
pred_years <- year_ids[t_cols]
n_years_pred <- length(t_cols)

year_scaled_pred <- as.vector(scale(year_ids)[t_cols])
year_mat_pred <- matrix(year_scaled_pred, nrow = J_pred, ncol = n_years_pred, byrow = TRUE)

# Trimmed spatial geometry for plotting
nw_grid_trim <- nw_grid_pred %>% select(sample_unit_id, state, geometry)

# Accumulators for all species --------------------------------------------
all_means_list  <- list()
all_widths_list <- list()

# Loop Over Species and Prior Types ---------------------------------------
for (spp in possible_bats) {
  for (prior in prior_types) {
    # Check for fit file
    filename_spatial <- paste0(spp, "_stPGOccfit_", prior, "_spatial.rds")
    filename_plain   <- paste0(spp, "_stPGOccfit_", prior, ".rds")
    
    fit_path <- if (file.exists(file.path(fits_dir, filename_spatial))) {
      file.path(fits_dir, filename_spatial)
    } else if (file.exists(file.path(fits_dir, filename_plain))) {
      file.path(fits_dir, filename_plain)
    } else {
      NULL
    }
    
    if (is.null(fit_path)) {
      cat(sprintf("Skipping %s (%s): fit file not found.\n", spp, prior))
      next
    }
    
    cat(sprintf("\n=== Predicting & Mapping stPGOcc: %s (%s) ===\n", spp, prior))
    fit <- readRDS(fit_path)
    
    is_cliff <- spp %in% cliff_spp
    p_occ <- if (is_cliff) 6 else 5
    
    # Construct 3D prediction design matrix X.0 (J_pred x n_years_pred x p_occ)
    X_0 <- array(1, dim = c(J_pred, n_years_pred, p_occ))
    X_0[, , 2] <- covars_pred$p_forest
    X_0[, , 3] <- covars_pred$precip
    X_0[, , 4] <- covars_pred$dem_max
    
    if (is_cliff) {
      X_0[, , 5] <- covars_pred$cliff_cover
      X_0[, , 6] <- year_mat_pred
    } else {
      X_0[, , 5] <- year_mat_pred
    }
    
    cat("  Running spatial prediction over full grid (J =", J_pred, "cells)...\n")
    out_pred <- tryCatch({
      predict(
        object    = fit,
        X.0       = X_0,
        t.cols    = t_cols,
        ignore.RE = TRUE,
        type      = "occupancy",
        coords.0  = coords_pred
      )
    }, error = function(e) {
      warning(sprintf("Prediction failed for %s (%s): %s", spp, prior, e$message))
      NULL
    })
    
    if (is.null(out_pred)) {
      rm(fit)
      gc()
      next
    }
    
    # Compute summary statistics: mean and 95% CI width
    cat("  Summarizing posterior predictions...\n")
    # psi.0.samples dim: (n_samples, J_pred, n_years_pred)
    psi_samples <- out_pred$psi.0.samples
    
    mean_psi_mat  <- apply(psi_samples, c(2, 3), mean)
    width_psi_mat <- apply(psi_samples, c(2, 3), function(x) {
      quantile(x, probs = 0.975, names = FALSE) - quantile(x, probs = 0.025, names = FALSE)
    })
    sd_psi_mat <- apply(psi_samples, c(2, 3), sd)
    
    colnames(mean_psi_mat)  <- paste0("year_", pred_years)
    colnames(width_psi_mat) <- paste0("year_", pred_years)
    colnames(sd_psi_mat)    <- paste0("year_", pred_years)
    
    # Format mean predictions
    mean_df <- as_tibble(mean_psi_mat) %>%
      mutate(sample_unit_id = nw_grid_pred$sample_unit_id) %>%
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
    width_df <- as_tibble(width_psi_mat) %>%
      mutate(sample_unit_id = nw_grid_pred$sample_unit_id) %>%
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
    
    prior_label <- if_else(prior == "informative", "Informative", "Weakly Informative")
    
    # 1. Mean Occupancy Map -------------------------------------------------
    cat("  Generating mean occupancy map...\n")
    p_mean <- ggplot(mean_df) +
      geom_sf(aes(fill = occ_prob), linewidth = 0) +
      facet_wrap(~year) +
      scale_fill_viridis_c(
        option = "D",
        limits = c(0, 1),
        name = expression("Predicted Mean" ~ hat(bar(psi)))
      ) +
      labs(
        title = paste0("stPGOcc: Predicted Posterior Mean Occupancy (", toupper(spp), ")"),
        subtitle = paste0("Priors: ", prior_label, " | Pacific Northwest Grid"),
        caption = "Predictions with spatial processes and fixed effects"
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
    
    ggsave(
      filename = file.path(maps_dir, paste0("map_mean_psi_", spp, "_", prior, ".png")),
      plot = p_mean,
      width = 9, height = 5, units = "in", dpi = 300
    )
    
    # 2. Posterior Uncertainty / 95% CI Width Map ---------------------------
    cat("  Generating uncertainty map...\n")
    p_width <- ggplot(width_df) +
      geom_sf(aes(fill = ci_width), linewidth = 0) +
      facet_wrap(~year) +
      scale_fill_viridis_c(
        option = "magma",
        limits = c(0, 1),
        name = "95% CI Width"
      ) +
      labs(
        title = paste0("stPGOcc: Posterior Uncertainty / 95% CI Width (", toupper(spp), ")"),
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
    
    ggsave(
      filename = file.path(maps_dir, paste0("map_width_psi_", spp, "_", prior, ".png")),
      plot = p_width,
      width = 9, height = 5, units = "in", dpi = 300
    )
    
    rm(fit, out_pred, psi_samples)
    gc()
  }
}

# Save Combined Prediction Summaries --------------------------------------
if (length(all_means_list) > 0) {
  all_means_df  <- bind_rows(all_means_list)
  all_widths_df <- bind_rows(all_widths_list)
  
  saveRDS(all_means_df, file.path(maps_dir, "spatial_mean_predictions.rds"))
  saveRDS(all_widths_df, file.path(maps_dir, "spatial_width_predictions.rds"))
  
  cat("\nConsolidated spatial predictions saved to:", maps_dir, "\n")
} else {
  cat("\nNo spatial predictions were generated (no model fits found).\n")
}

cat("Done generating spatial maps.\n")
