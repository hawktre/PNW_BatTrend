renv::activate()
## ---------------------------
## Purpose of script: Read in Stan model fits across species and prior types
##                    from the fits directory, run posterior predictive checks (PPC),
##                    extract parameter summaries (occ, det, dyn), compute regional
##                    occupancy trajectories and trends, and generate diagnostic plots
##                    organized into prior-specific subfolders.
##
## Author: Trent VanHawkins
## ---------------------------

# Load Packages -----------------------------------------------------------
library(tidyverse)
library(here)
library(cmdstanr)
library(bayesplot)
library(posterior)
library(tidybayes)
library(sf)

# Directories Setup -------------------------------------------------------
fits_dir <- here("data/processed/results/stan/fits")
base_diag_dir <- here("data/processed/results/stan/diagnostics")

prior_types <- c("weakly_informative", "informative", "weakly_informative_sensitivity")

# Create base diagnostics directory and subdirectories for each prior type
if (!dir.exists(base_diag_dir)) {
  dir.create(base_diag_dir, recursive = TRUE, showWarnings = FALSE)
}

for (prior in prior_types) {
  prior_dir <- file.path(base_diag_dir, prior)
  if (!dir.exists(prior_dir)) {
    dir.create(prior_dir, recursive = TRUE, showWarnings = FALSE)
  }
}

# Read in Metadata & Design Matrices --------------------------------------
des_matrix_path <- here("data/processed/results/stan/design_matrices.rds")
des_matrix <- if (file.exists(des_matrix_path)) readRDS(des_matrix_path) else NULL

des_matrix_sens_path <- here("data/processed/results/stan/design_matrices_sensitivity.rds")
des_matrix_sens <- if (file.exists(des_matrix_sens_path)) readRDS(des_matrix_sens_path) else NULL

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

# Cliff-associated species
cliff_spp <- c("anpa", "euma", "myci", "pahe")

# Load bundled data objects for PPC checks
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

# Spatial Setup for Moran's I --------------------------------------------
nw_grid <- readRDS(here("data/processed/occurrence/nw_grid.rds"))
index_keys <- readRDS(here("data/processed/results/stan/index_keys.rds"))
index_keys_sens_path <- here("data/processed/results/stan/index_keys_sensitivity.rds")
index_keys_sens <- if (file.exists(index_keys_sens_path)) readRDS(index_keys_sens_path) else NULL

# Full grid distance matrix
grid_full <- nw_grid %>% filter(sample_unit_id %in% index_keys$site_ids) %>% arrange(match(sample_unit_id, index_keys$site_ids))
cents_full <- suppressWarnings(sf::st_centroid(grid_full))
d_full <- units::drop_units(sf::st_distance(cents_full)) / 1000 # in km
diag(d_full) <- 0

# Sensitivity grid distance matrix
d_sens <- if (!is.null(index_keys_sens)) {
  grid_sens <- nw_grid %>% filter(sample_unit_id %in% index_keys_sens$site_ids) %>% arrange(match(sample_unit_id, index_keys_sens$site_ids))
  cents_sens <- suppressWarnings(sf::st_centroid(grid_sens))
  d_s <- units::drop_units(sf::st_distance(cents_sens)) / 1000 # in km
  diag(d_s) <- 0
  d_s
} else {
  NULL
}

# Distance bands for spatial correlogram (10 km to 50 km)
distance_bands <- list(
  list(min = 0, max = 10, label = "0-10 km", mid = 5),
  list(min = 10, max = 20, label = "10-20 km", mid = 15),
  list(min = 20, max = 30, label = "20-30 km", mid = 25),
  list(min = 30, max = 40, label = "30-40 km", mid = 35),
  list(min = 40, max = 50, label = "40-50 km", mid = 45)
)

# Vectorized helper function to compute Moran's I across MCMC draws
calc_moran_draws <- function(Z_mat, dist_mat, d_min, d_max, row_standardize = TRUE) {
  n <- ncol(Z_mat)
  W <- matrix(0, nrow = n, ncol = n)
  if (d_min == 0) {
    W[dist_mat <= d_max & dist_mat > 0] <- 1
  } else {
    W[dist_mat > d_min & dist_mat <= d_max] <- 1
  }
  
  if (sum(W) == 0) return(rep(NA_real_, nrow(Z_mat)))
  
  if (row_standardize) {
    rs <- rowSums(W)
    rs[rs == 0] <- 1
    W <- W / rs
  }
  s0 <- sum(W)
  
  Z_centered <- Z_mat - rowMeans(Z_mat)
  numer <- rowSums((Z_centered %*% t(W)) * Z_centered)
  denom <- rowSums(Z_centered^2)
  
  return((n / s0) * (numer / denom))
}

# Storage for results -----------------------------------------------------
all_res_list    <- list()
psi_full_list   <- list()
trend_full_list <- list()
p_overall_list  <- list()
p_annual_list   <- list()
moran_i_list    <- list()

# Loop Over Species and Priors --------------------------------------------
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
      cat(sprintf("Skipping Stan fit for %s (%s): file not found in %s\n", spp, prior, fits_dir))
      next
    }
    
    cat(sprintf("\n=== Processing Stan Diagnostics: %s (%s) ===\n", spp, prior))
    fit <- readRDS(fit_path)
    
    prior_out_dir <- file.path(base_diag_dir, prior)
    
    # 1. Posterior Predictive Check (PPC) -----------------------------------
    cat("  Computing PPC...\n")
    tryCatch({
      y_rep <- fit$draws("y_rep", format = "matrix")
      spp_stan_data <- stan_data_list[[prior]]
      
      y <- if (!is.null(spp_stan_data) && spp %in% names(spp_stan_data)) {
        spp_stan_data[[spp]]$dets
      } else {
        NULL
      }
      
      if (!is.null(y)) {
        ppc <- bayesplot::ppc_bars(y, y_rep) +
          labs(
            title = paste0("Posterior Predictive Check (", toupper(spp), ", ", gsub("_", " ", prior), ")"),
            x = "Detections",
            y = "Count"
          ) +
          theme_bw()
        
        ggsave(
          filename = file.path(prior_out_dir, paste0("ppc_", spp, "_", prior, ".png")),
          plot = ppc, width = 5, height = 4, units = "in", dpi = 300
        )
      }
    }, error = function(e) {
      warning(sprintf("PPC calculation failed for %s (%s): %s", spp, prior, e$message))
    })
    
    # 2. Extract Parameter Summaries ----------------------------------------
    cat("  Extracting parameter summaries...\n")
    occ_params <- c("alpha0", "alphas")
    det_params <- c("beta0", "betas")
    dyn_params <- c("a", "b", "gamma", "phi")
    
    occ_summary <- tryCatch({
      fit$summary(
        variables = occ_params,
        posterior::default_summary_measures()[1:4],
        quantiles = ~ quantile2(., probs = c(0.025, 0.25, 0.5, 0.75, 0.975)),
        posterior::default_convergence_measures()
      )
    }, error = function(e) NULL)
    
    det_summary <- tryCatch({
      fit$summary(
        variables = det_params,
        posterior::default_summary_measures()[1:4],
        quantiles = ~ quantile2(., probs = c(0.025, 0.25, 0.5, 0.75, 0.975)),
        posterior::default_convergence_measures()
      )
    }, error = function(e) NULL)
    
    dyn_summary <- tryCatch({
      fit$summary(
        variables = dyn_params,
        posterior::default_summary_measures()[1:4],
        quantiles = ~ quantile2(., probs = c(0.025, 0.25, 0.5, 0.75, 0.975)),
        posterior::default_convergence_measures()
      )
    }, error = function(e) NULL)
    
    des_mat_use <- if (prior == "weakly_informative_sensitivity" && !is.null(des_matrix_sens)) {
      des_matrix_sens
    } else {
      des_matrix
    }
    
    # Assign variable names if design matrix is available
    if (!is.null(occ_summary) && !is.null(des_mat_use)) {
      if (spp %in% cliff_spp) {
        occ_summary$variable <- c("Intercept", colnames(des_mat_use$xmatc))
      } else {
        occ_summary$variable <- c("Intercept", colnames(des_mat_use$xmata))
      }
      occ_summary$var_type <- "occurrence"
    }
    
    if (!is.null(det_summary) && !is.null(des_mat_use)) {
      det_summary$variable <- c("Intercept", colnames(des_mat_use$vmat_all))
      det_summary$var_type <- "detection"
    }
    
    if (!is.null(dyn_summary)) {
      dyn_summary$var_type <- "dynamics"
    }
    
    tmp_res <- bind_rows(occ_summary, det_summary, dyn_summary) %>%
      mutate(prior_type = prior, species = spp)
    
    all_res_list[[paste(spp, prior, sep = "_")]] <- tmp_res
    
    # 3. Regional Occupancy Trajectories & Trend ----------------------------
    cat("  Computing regional psi and trend...\n")
    tryCatch({
      psi_draws <- tidybayes::spread_draws(fit, psi[site, year])
      
      psi_summary <- psi_draws %>%
        group_by(.draw, year) %>%
        summarise(
          regional_psi = mean(psi),
          .groups = "drop"
        )
      
      psi_regional_summary <- psi_summary %>%
        group_by(year) %>%
        summarise(
          mean_regional = mean(regional_psi),
          sd_regional   = sd(regional_psi),
          q2.5          = quantile2(regional_psi, probs = 0.025),
          q25           = quantile2(regional_psi, probs = 0.25),
          q50           = quantile2(regional_psi, probs = 0.50),
          q75           = quantile2(regional_psi, probs = 0.75),
          q97.5         = quantile2(regional_psi, probs = 0.975),
          .groups       = "drop"
        ) %>%
        mutate(prior_type = prior, species = spp)
      
      psi_full_list[[paste(spp, prior, sep = "_")]] <- psi_regional_summary
      
      # Trend calculation per draw
      trend_draws <- psi_summary %>%
        group_by(.draw) %>%
        summarise(
          delta_trend = coef(lm(regional_psi ~ year))[2],
          .groups = "drop"
        )
      
      delta_trend_summary <- trend_draws %>%
        summarise(
          mean_trend = mean(delta_trend),
          sd_trend   = sd(delta_trend),
          q2.5       = quantile(delta_trend, probs = 0.025, names = FALSE),
          q25        = quantile(delta_trend, probs = 0.25, names = FALSE),
          q50        = quantile(delta_trend, probs = 0.50, names = FALSE),
          q75        = quantile(delta_trend, probs = 0.75, names = FALSE),
          q97.5      = quantile(delta_trend, probs = 0.975, names = FALSE),
          .groups    = "drop"
        ) %>%
        mutate(
          prior_type = prior,
          species = spp,
          trend_status = case_when(
            q2.5 > 0 ~ "Increasing",
            q97.5 < 0 ~ "Decreasing",
            TRUE ~ "Uncertain"
          )
        )
      
      trend_full_list[[paste(spp, prior, sep = "_")]] <- delta_trend_summary
    }, error = function(e) {
      warning(sprintf("Psi/Trend calculation failed for %s (%s): %s", spp, prior, e$message))
    })
    
    # 4. Detection Probability Summaries (Overall & Annual) -----------------
    cat("  Computing detection probability summaries...\n")
    tryCatch({
      p_mat <- fit$draws("p", format = "matrix")
      
      # Extract observation-to-year mapping
      spp_stan_data <- stan_data_list[[prior]]
      if (!is.null(spp_stan_data) && spp %in% names(spp_stan_data)) {
        st <- spp_stan_data[[spp]]$start_idx
        et <- spp_stan_data[[spp]]$end_idx
        n_obs <- spp_stan_data[[spp]]$n_obs
        obs_year <- integer(n_obs)
        for (i in 1:nrow(st)) {
          for (t in 1:ncol(st)) {
            if (st[i, t] > 0) {
              obs_year[st[i, t]:et[i, t]] <- t + 2015
            }
          }
        }
      } else {
        obs_year <- NULL
      }
      
      # (a) Overall detection probability per draw across all survey observations
      overall_p_draws <- rowMeans(p_mat)
      beta0_draws <- fit$draws("beta0", format = "matrix")[, 1]
      p0_draws <- plogis(beta0_draws)
      p_star4_draws <- 1 - (1 - overall_p_draws)^4
      
      p_overall_summary <- tibble(
        species       = spp,
        prior_type    = prior,
        mean_p        = mean(overall_p_draws),
        sd_p          = sd(overall_p_draws),
        q2.5_p        = quantile(overall_p_draws, probs = 0.025, names = FALSE),
        q25_p         = quantile(overall_p_draws, probs = 0.25, names = FALSE),
        q50_p         = quantile(overall_p_draws, probs = 0.50, names = FALSE),
        q75_p         = quantile(overall_p_draws, probs = 0.75, names = FALSE),
        q97.5_p       = quantile(overall_p_draws, probs = 0.975, names = FALSE),
        mean_p0       = mean(p0_draws),
        sd_p0         = sd(p0_draws),
        q2.5_p0       = quantile(p0_draws, probs = 0.025, names = FALSE),
        q97.5_p0      = quantile(p0_draws, probs = 0.975, names = FALSE),
        mean_p_star4  = mean(p_star4_draws),
        sd_p_star4    = sd(p_star4_draws),
        q2.5_p_star4  = quantile(p_star4_draws, probs = 0.025, names = FALSE),
        q97.5_p_star4 = quantile(p_star4_draws, probs = 0.975, names = FALSE)
      )
      p_overall_list[[paste(spp, prior, sep = "_")]] <- p_overall_summary
      
      # (b) Annual detection probability per draw
      if (!is.null(obs_year)) {
        unique_years <- sort(unique(obs_year))
        annual_p_summaries <- map_dfr(unique_years, function(yr) {
          yr_idx <- which(obs_year == yr)
          yr_p_draws <- rowMeans(p_mat[, yr_idx, drop = FALSE])
          tibble(
            species    = spp,
            prior_type = prior,
            year       = yr,
            n_obs      = length(yr_idx),
            mean_p     = mean(yr_p_draws),
            sd_p       = sd(yr_p_draws),
            q2.5       = quantile(yr_p_draws, probs = 0.025, names = FALSE),
            q25        = quantile(yr_p_draws, probs = 0.25, names = FALSE),
            q50        = quantile(yr_p_draws, probs = 0.50, names = FALSE),
            q75        = quantile(yr_p_draws, probs = 0.75, names = FALSE),
            q97.5      = quantile(yr_p_draws, probs = 0.975, names = FALSE)
          )
        })
        p_annual_list[[paste(spp, prior, sep = "_")]] <- annual_p_summaries
      }
      
      rm(p_mat)
    }, error = function(e) {
      warning(sprintf("Detection probability calculation failed for %s (%s): %s", spp, prior, e$message))
    })
    
    # 5. Moran's I Residual Spatial Autocorrelation --------------------------
    cat("  Computing Moran's I correlogram (10-50 km)...\n")
    tryCatch({
      if ("occ_res" %in% fit$metadata()$stan_variables) {
        occ_res_mat <- fit$draws("occ_res", format = "matrix")
        n_draws <- nrow(occ_res_mat)
        
        is_sens <- prior == "weakly_informative_sensitivity"
        d_mat_use <- if (is_sens && !is.null(d_sens)) d_sens else d_full
        n_sites_use <- nrow(d_mat_use)
        n_years_use <- ncol(occ_res_mat) / n_sites_use
        
        occ_res_arr <- array(as.vector(occ_res_mat), dim = c(n_draws, n_sites_use, n_years_use))
        site_res_per_draw <- apply(occ_res_arr, c(1, 2), mean)
        
        moran_spp_list <- map_dfr(distance_bands, function(b) {
          i_draws <- calc_moran_draws(site_res_per_draw, d_mat_use, b$min, b$max)
          tibble(
            species      = spp,
            prior_type   = prior,
            band_label   = b$label,
            dist_min     = b$min,
            dist_max     = b$max,
            dist_mid     = b$mid,
            mean_I       = mean(i_draws, na.rm = TRUE),
            sd_I         = sd(i_draws, na.rm = TRUE),
            q2.5         = quantile(i_draws, 0.025, na.rm = TRUE, names = FALSE),
            q25          = quantile(i_draws, 0.25, na.rm = TRUE, names = FALSE),
            q50          = quantile(i_draws, 0.50, na.rm = TRUE, names = FALSE),
            q75          = quantile(i_draws, 0.75, na.rm = TRUE, names = FALSE),
            q97.5        = quantile(i_draws, 0.975, na.rm = TRUE, names = FALSE)
          )
        })
        moran_i_list[[paste(spp, prior, sep = "_")]] <- moran_spp_list
        rm(occ_res_mat, occ_res_arr, site_res_per_draw)
      }
    }, error = function(e) {
      warning(sprintf("Moran's I calculation failed for %s (%s): %s", spp, prior, e$message))
    })
    
    rm(fit)
    gc()
  }
}

# Bind Consolidated Results -----------------------------------------------
all_res_df     <- bind_rows(all_res_list)
psi_full_df    <- bind_rows(psi_full_list)
trend_full_df  <- bind_rows(trend_full_list)
p_overall_df   <- bind_rows(p_overall_list)
p_annual_df    <- bind_rows(p_annual_list)
moran_i_df     <- bind_rows(moran_i_list)

if (nrow(all_res_df) == 0) {
  stop("No Stan fits were found to summarize. Please check that fits exist in data/processed/results/stan/fits/.")
}

# Save Summary Tables -----------------------------------------------------
saveRDS(all_res_df, file.path(base_diag_dir, "all_res_stan.rds"))
saveRDS(psi_full_df, file.path(base_diag_dir, "psi_full_stan.rds"))
saveRDS(trend_full_df, file.path(base_diag_dir, "trend_full_stan.rds"))
if (nrow(p_overall_df) > 0) saveRDS(p_overall_df, file.path(base_diag_dir, "p_overall_stan.rds"))
if (nrow(p_annual_df) > 0) saveRDS(p_annual_df, file.path(base_diag_dir, "p_annual_stan.rds"))
if (nrow(moran_i_df) > 0)  saveRDS(moran_i_df, file.path(base_diag_dir, "moran_i_stan.rds"))

write_csv(all_res_df, file.path(base_diag_dir, "stan_diagnostics_summary.csv"))
write_csv(trend_full_df, file.path(base_diag_dir, "stan_trend_summary.csv"))
if (nrow(p_overall_df) > 0) write_csv(p_overall_df, file.path(base_diag_dir, "stan_p_overall_summary.csv"))
if (nrow(p_annual_df) > 0) write_csv(p_annual_df, file.path(base_diag_dir, "stan_p_annual_summary.csv"))
if (nrow(moran_i_df) > 0)  write_csv(moran_i_df, file.path(base_diag_dir, "stan_moran_i_summary.csv"))

cat("\nSummary tables saved to:", base_diag_dir, "\n")

# Multi-Species Visualization Figures -------------------------------------

## 1. Occurrence Parameters Summary ---------------------------------------
cat("Generating occurrence parameter summary plots...\n")
occ_plot_data <- all_res_df %>%
  filter(var_type == "occurrence") %>%
  mutate(
    prior_label = case_when(
      prior_type == "informative" ~ "Informative",
      prior_type == "weakly_informative" ~ "Weakly Informative",
      TRUE ~ "Weakly Informative (Sensitivity)"
    ),
    variable_label = factor(
      variable,
      levels = c("Intercept", "log_fc", "precip", "dem_max", "log_cliff"),
      labels = c("Intercept", "Log(Forest Cover %)", "Precip (mm)", "Elevation (m)", "Log(Cliff Cover %)")
    )
  )

# All priors (3 priors)
ggplot(occ_plot_data, aes(x = variable_label, y = mean, group = prior_label, color = prior_label)) +
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0, position = position_dodge(width = 0.35)) +
  geom_errorbar(aes(ymin = q25, ymax = q75), width = 0, position = position_dodge(width = 0.35), linewidth = 1.1) +
  geom_point(position = position_dodge(width = 0.35), size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~species, ncol = 2) +
  labs(
    title = "Occurrence Model Covariate Comparison",
    x = "Covariate",
    y = "Mean (50% and 95% Credible Interval)",
    color = "Prior"
  ) +
  theme_bw() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(base_diag_dir, "occ_pars_summary.png"), width = 14, height = 12, units = "in", dpi = 300)

# Main priors only (without sensitivity)
occ_plot_data_main <- occ_plot_data %>% filter(prior_type %in% c("informative", "weakly_informative"))
ggplot(occ_plot_data_main, aes(x = variable_label, y = mean, group = prior_label, color = prior_label)) +
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0, position = position_dodge(width = 0.35)) +
  geom_errorbar(aes(ymin = q25, ymax = q75), width = 0, position = position_dodge(width = 0.35), linewidth = 1.1) +
  geom_point(position = position_dodge(width = 0.35), size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~species, ncol = 2) +
  labs(
    title = "Occurrence Model Covariate Comparison (Informative vs. Weakly Informative)",
    x = "Covariate",
    y = "Mean (50% and 95% Credible Interval)",
    color = "Prior"
  ) +
  theme_bw() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(base_diag_dir, "occ_pars_summary_main.png"), width = 14, height = 12, units = "in", dpi = 300)
ggsave(file.path(base_diag_dir, "occ_pars_summary_nosens.png"), width = 14, height = 12, units = "in", dpi = 300)

## 2. Detection Parameters Summary ----------------------------------------
cat("Generating detection parameter summary plots (informative prior only)...\n")
det_plot_data <- all_res_df %>%
  filter(var_type == "detection", prior_type == "informative") %>%
  mutate(
    variable_label = factor(
      variable,
      levels = c("Intercept", "clutter_percent1", "clutter_percent2", "clutter_percent3", "clutter_percent4", "scale(tmin)", "scale(dayl)", "water_ind"),
      labels = c("Intercept", "Clutter 0-25%", "Clutter 26-50%", "Clutter 51-75%", "Clutter 76-100%", "Nightly Min Temp", "Day Length", "Waterbody")
    )
  )

# Plot detection parameters for informative prior
ggplot(det_plot_data, aes(x = variable_label, y = mean)) +
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0, color = "#2c3e50") +
  geom_errorbar(aes(ymin = q25, ymax = q75), width = 0, color = "#2980b9", linewidth = 1.1) +
  geom_point(color = "#2c3e50", size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~species, ncol = 2) +
  labs(
    title = "Detection Model Covariates (Informative Prior)",
    x = "Covariate",
    y = "Mean (50% and 95% Credible Interval)"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(base_diag_dir, "det_pars_summary.png"), width = 14, height = 12, units = "in", dpi = 300)
ggsave(file.path(base_diag_dir, "informative", "det_pars_informative.png"), width = 14, height = 12, units = "in", dpi = 300)

## 3. Dynamics Parameters (Phi & Gamma) -----------------------------------
cat("Generating dynamics parameter summary plots...\n")
dyn_phi_data <- all_res_df %>%
  filter(var_type == "dynamics", str_detect(variable, "phi")) %>%
  mutate(
    year = parse_number(variable) + 2016,
    prior_label = case_when(
      prior_type == "informative" ~ "Informative",
      prior_type == "weakly_informative" ~ "Weakly Informative",
      TRUE ~ "Weakly Informative (Sensitivity)"
    )
  )

# Phi (all priors) with dodge and 0-1 limits
ggplot(dyn_phi_data, aes(x = factor(year), y = mean, group = prior_label, color = prior_label)) +
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0.25, position = position_dodge(width = 0.35)) +
  geom_point(position = position_dodge(width = 0.35), size = 1.5) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  facet_wrap(~species, ncol = 2) +
  labs(
    title = "Regional Dynamics (Persistence / Phi)",
    x = "Year",
    y = "Mean Persistence Probability (95% CI)",
    color = "Prior"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(file.path(base_diag_dir, "dyn_phi_summary.png"), width = 14, height = 12, units = "in", dpi = 300)

# Phi (main priors only)
dyn_phi_data_main <- dyn_phi_data %>% filter(prior_type %in% c("informative", "weakly_informative"))
ggplot(dyn_phi_data_main, aes(x = factor(year), y = mean, group = prior_label, color = prior_label)) +
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0.25, position = position_dodge(width = 0.35)) +
  geom_point(position = position_dodge(width = 0.35), size = 1.5) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  facet_wrap(~species, ncol = 2) +
  labs(
    title = "Regional Dynamics (Persistence / Phi - Informative vs. Weakly Informative)",
    x = "Year",
    y = "Mean Persistence Probability (95% CI)",
    color = "Prior"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(file.path(base_diag_dir, "dyn_phi_summary_main.png"), width = 14, height = 12, units = "in", dpi = 300)
ggsave(file.path(base_diag_dir, "dyn_phi_summary_nosens.png"), width = 14, height = 12, units = "in", dpi = 300)

dyn_gamma_data <- all_res_df %>%
  filter(var_type == "dynamics", str_detect(variable, "gamma")) %>%
  mutate(
    year = parse_number(variable) + 2016,
    prior_label = case_when(
      prior_type == "informative" ~ "Informative",
      prior_type == "weakly_informative" ~ "Weakly Informative",
      TRUE ~ "Weakly Informative (Sensitivity)"
    )
  )

# Gamma (all priors) with dodge and 0-1 limits
ggplot(dyn_gamma_data, aes(x = factor(year), y = mean, group = prior_label, color = prior_label)) +
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0.25, position = position_dodge(width = 0.35)) +
  geom_point(position = position_dodge(width = 0.35), size = 1.5) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  facet_wrap(~species, ncol = 2) +
  labs(
    title = "Regional Dynamics (Colonization / Gamma)",
    x = "Year",
    y = "Mean Colonization Probability (95% CI)",
    color = "Prior"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(file.path(base_diag_dir, "dyn_gamma_summary.png"), width = 14, height = 12, units = "in", dpi = 300)

# Gamma (main priors only)
dyn_gamma_data_main <- dyn_gamma_data %>% filter(prior_type %in% c("informative", "weakly_informative"))
ggplot(dyn_gamma_data_main, aes(x = factor(year), y = mean, group = prior_label, color = prior_label)) +
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0.25, position = position_dodge(width = 0.35)) +
  geom_point(position = position_dodge(width = 0.35), size = 1.5) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  facet_wrap(~species, ncol = 2) +
  labs(
    title = "Regional Dynamics (Colonization / Gamma - Informative vs. Weakly Informative)",
    x = "Year",
    y = "Mean Colonization Probability (95% CI)",
    color = "Prior"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(file.path(base_diag_dir, "dyn_gamma_summary_main.png"), width = 14, height = 12, units = "in", dpi = 300)
ggsave(file.path(base_diag_dir, "dyn_gamma_summary_nosens.png"), width = 14, height = 12, units = "in", dpi = 300)

## 4. Occupancy Probability Trajectories over Time -----------------------
cat("Generating regional occupancy trajectory plots...\n")
psi_trend_merged <- psi_full_df %>%
  mutate(cal_year = if_else(year <= 10, year + 2015, year)) %>%
  left_join(
    trend_full_df %>% select(species, prior_type, trend_status, mean_trend),
    by = c("species", "prior_type")
  ) %>%
  group_by(species, prior_type) %>%
  mutate(
    trend_fitted = mean(mean_regional) + mean_trend * (cal_year - mean(cal_year))
  ) %>%
  ungroup()

for (p_type in prior_types) {
  p_label <- case_when(
    p_type == "informative" ~ "Informative",
    p_type == "weakly_informative" ~ "Weakly Informative",
    TRUE ~ "Weakly Informative (Sensitivity)"
  )
  
  psi_p_data <- psi_trend_merged %>% filter(prior_type == p_type)
  
  if (nrow(psi_p_data) > 0) {
    plt <- ggplot(psi_p_data, aes(x = factor(cal_year), y = mean_regional, group = species)) +
      geom_ribbon(aes(ymin = q2.5, ymax = q97.5, fill = trend_status), alpha = 0.25) +
      geom_line(aes(color = trend_status), linewidth = 0.85) +
      geom_line(aes(y = trend_fitted), color = "black", alpha = 0.85, linewidth = 0.85) +
      geom_point(aes(color = trend_status), size = 1.5) +
      facet_wrap(~species, ncol = 2) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
      scale_color_manual(values = c("Increasing" = "#27ae60", "Decreasing" = "#c0392b", "Uncertain" = "#7f8c8d")) +
      scale_fill_manual(values = c("Increasing" = "#27ae60", "Decreasing" = "#c0392b", "Uncertain" = "#7f8c8d")) +
      labs(
        title = paste0("Regional Occupancy Trajectories (", p_label, " Priors)"),
        subtitle = "Colored lines & ribbon: annual posterior mean occupancy (95% CI) | Straight black line: estimated linear trend",
        x = "Year",
        y = "Regional Occupancy Probability (95% CI)",
        color = "Trend Status",
        fill = "Trend Status"
      ) +
      theme_bw() +
      theme(legend.position = "bottom")
    
    prior_out_dir <- file.path(base_diag_dir, p_type)
    ggsave(
      file.path(prior_out_dir, paste0("psi_trend_", p_type, ".png")),
      plot = plt, width = 14, height = 12, units = "in", dpi = 300
    )
    ggsave(
      file.path(base_diag_dir, paste0("psi_trend_", p_type, ".png")),
      plot = plt, width = 14, height = 12, units = "in", dpi = 300
    )
  }
}

## 5. Trend Estimate Comparison Across Prior Types ------------------------
cat("Generating trend estimate comparison plots...\n")
trend_compare_plot <- trend_full_df %>%
  mutate(
    prior_label = case_when(
      prior_type == "informative" ~ "Informative",
      prior_type == "weakly_informative" ~ "Weakly Informative",
      TRUE ~ "Weakly Informative (Sensitivity)"
    )
  ) %>%
  ggplot(aes(x = species, y = mean_trend, color = prior_label)) +
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0.35, position = position_dodge(0.5), linewidth = 0.9) +
  geom_point(position = position_dodge(0.5), size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "Trend Estimate Comparison Across Prior Types",
    x = "Species",
    y = expression(hat(delta)[trend] ~ "(Annual Occupancy Slope)"),
    color = "Prior"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(file.path(base_diag_dir, "delta_trend_compare.png"), plot = trend_compare_plot, width = 12, height = 7, units = "in", dpi = 300)

# Main priors only (without sensitivity)
trend_compare_plot_main <- trend_full_df %>%
  filter(prior_type %in% c("informative", "weakly_informative")) %>%
  mutate(
    prior_label = case_when(
      prior_type == "informative" ~ "Informative",
      TRUE ~ "Weakly Informative"
    )
  ) %>%
  ggplot(aes(x = species, y = mean_trend, color = prior_label)) +
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0.35, position = position_dodge(0.5), linewidth = 0.9) +
  geom_point(position = position_dodge(0.5), size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "Trend Estimate Comparison (Informative vs. Weakly Informative)",
    x = "Species",
    y = expression(hat(delta)[trend] ~ "(Annual Occupancy Slope)"),
    color = "Prior"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(file.path(base_diag_dir, "delta_trend_compare_main.png"), plot = trend_compare_plot_main, width = 12, height = 7, units = "in", dpi = 300)
ggsave(file.path(base_diag_dir, "delta_trend_compare_nosens.png"), plot = trend_compare_plot_main, width = 12, height = 7, units = "in", dpi = 300)

## 6. Overall Species Detection Probabilities Summary --------------------
if (nrow(p_overall_df) > 0) {
  cat("Generating overall detection probability summary plots (informative prior only)...\n")
  p_overall_plot_data <- p_overall_df %>%
    filter(prior_type == "informative") %>%
    mutate(species_upper = toupper(species))
  
  # Informative prior detection probability summary across species
  p_overall_plot <- ggplot(p_overall_plot_data, aes(x = reorder(species_upper, mean_p), y = mean_p)) +
    geom_errorbar(aes(ymin = q2.5_p, ymax = q97.5_p), width = 0.35, color = "#2c3e50", linewidth = 0.8) +
    geom_errorbar(aes(ymin = q25_p, ymax = q75_p), width = 0, color = "#2980b9", linewidth = 1.3) +
    geom_point(color = "#2c3e50", size = 2) +
    coord_flip() +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      title = "Mean Single-Night Detection Probability (p) Across Species (Informative Prior)",
      x = "Species",
      y = "Mean Detection Probability (50% and 95% CI)"
    ) +
    theme_bw()
  
  ggsave(file.path(base_diag_dir, "p_overall_summary.png"), plot = p_overall_plot, width = 10, height = 7, units = "in", dpi = 300)
  ggsave(file.path(base_diag_dir, "informative", "p_overall_informative.png"), plot = p_overall_plot, width = 10, height = 7, units = "in", dpi = 300)
}

## 7. Annual Detection Probability Trajectories over Time -----------------
if (nrow(p_annual_df) > 0) {
  cat("Generating annual detection probability trajectory plots (informative prior only)...\n")
  p_annual_plot_data <- p_annual_df %>%
    filter(prior_type == "informative") %>%
    mutate(species_upper = toupper(species))
  
  # Informative prior annual realized detection probabilities
  p_annual_plot <- ggplot(p_annual_plot_data, aes(x = factor(year), y = mean_p, group = species_upper)) +
    geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.2, fill = "#2980b9") +
    geom_line(color = "#2980b9", linewidth = 0.8) +
    geom_point(color = "#2980b9", size = 1.5) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    facet_wrap(~species_upper, ncol = 2) +
    labs(
      title = "Annual Realized Detection Probabilities (2016-2025, Informative Prior)",
      x = "Year",
      y = "Mean Nightly Detection Probability (95% CI)"
    ) +
    theme_bw()
  
  ggsave(file.path(base_diag_dir, "p_annual_summary.png"), plot = p_annual_plot, width = 14, height = 12, units = "in", dpi = 300)
  ggsave(file.path(base_diag_dir, "informative", "p_annual_informative.png"), plot = p_annual_plot, width = 14, height = 12, units = "in", dpi = 300)
}

## 8. Moran's I Residual Spatial Correlograms (10 to 50 km) ----------------
if (nrow(moran_i_df) > 0) {
  cat("Generating Moran's I residual spatial correlogram plots...\n")
  moran_plot_data <- moran_i_df %>%
    mutate(
      prior_label = case_when(
        prior_type == "informative" ~ "Informative",
        prior_type == "weakly_informative" ~ "Weakly Informative",
        TRUE ~ "Weakly Informative (Sensitivity)"
      ),
      species_upper = toupper(species)
    )
  
  # All priors (3 priors)
  p_moran_all <- ggplot(moran_plot_data, aes(x = band_label, y = mean_I, group = prior_label, color = prior_label)) +
    geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0.25, position = position_dodge(width = 0.35)) +
    geom_point(position = position_dodge(width = 0.35), size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    facet_wrap(~species_upper, ncol = 2) +
    labs(
      title = "Residual Spatial Autocorrelation Correlogram (Moran's I)",
      x = "Neighborhood Distance Band",
      y = "Posterior Mean Moran's I (95% CI)",
      color = "Prior"
    ) +
    theme_bw() +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1))
  
  ggsave(file.path(base_diag_dir, "moran_correlogram_summary.png"), plot = p_moran_all, width = 14, height = 12, units = "in", dpi = 300)
  
  # Main priors only (without sensitivity)
  moran_plot_data_main <- moran_plot_data %>% filter(prior_type %in% c("informative", "weakly_informative"))
  p_moran_main <- ggplot(moran_plot_data_main, aes(x = band_label, y = mean_I, group = prior_label, color = prior_label)) +
    geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0.25, position = position_dodge(width = 0.35)) +
    geom_point(position = position_dodge(width = 0.35), size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    facet_wrap(~species_upper, ncol = 2) +
    labs(
      title = "Residual Spatial Autocorrelation Correlogram (Informative vs. Weakly Informative)",
      x = "Neighborhood Distance Band",
      y = "Posterior Mean Moran's I (95% CI)",
      color = "Prior"
    ) +
    theme_bw() +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1))
  
  ggsave(file.path(base_diag_dir, "moran_correlogram_summary_main.png"), plot = p_moran_main, width = 14, height = 12, units = "in", dpi = 300)
  ggsave(file.path(base_diag_dir, "moran_correlogram_summary_nosens.png"), plot = p_moran_main, width = 14, height = 12, units = "in", dpi = 300)
  
  # Prior-specific correlogram plots
  for (p_type in prior_types) {
    p_sub <- moran_plot_data %>% filter(prior_type == p_type)
    if (nrow(p_sub) > 0) {
      p_label <- case_when(
        p_type == "informative" ~ "Informative",
        p_type == "weakly_informative" ~ "Weakly Informative",
        TRUE ~ "Weakly Informative (Sensitivity)"
      )
      p_moran_sub <- ggplot(p_sub, aes(x = band_label, y = mean_I, group = species_upper)) +
        geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0.25, color = "#2980b9") +
        geom_point(color = "#2980b9", size = 1.5) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
        facet_wrap(~species_upper, ncol = 2) +
        labs(
          title = paste0("Residual Spatial Correlogram - Moran's I (", p_label, ")"),
          x = "Neighborhood Distance Band",
          y = "Posterior Mean Moran's I (95% CI)"
        ) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 30, hjust = 1))
      
      prior_out_dir <- file.path(base_diag_dir, p_type)
      ggsave(file.path(prior_out_dir, paste0("moran_correlogram_", p_type, ".png")), plot = p_moran_sub, width = 14, height = 12, units = "in", dpi = 300)
      ggsave(file.path(base_diag_dir, paste0("moran_correlogram_", p_type, ".png")), plot = p_moran_sub, width = 14, height = 12, units = "in", dpi = 300)
    }
  }
}

cat("\nAll Stan diagnostics and figures completed successfully.\n")
