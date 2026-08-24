renv::activate()
## ---------------------------
## Purpose of script: Read in stPGOcc model fits across species and prior types
##                    (full and sensitivity analyses), run posterior predictive checks (PPC),
##                    extract parameter summaries (occupancy, detection, spatial), compute regional
##                    occupancy trajectories and trends, and generate diagnostic plots.
##
## Author: Trent VanHawkins
## ---------------------------

# Load Packages -----------------------------------------------------------
library(tidyverse)
library(here)
library(spOccupancy)
library(coda)
library(bayesplot)
library(posterior)

# Directories Setup -------------------------------------------------------
diag_dir <- here("data/processed/results/stPGOcc/full/diagnostics")
if (!dir.exists(diag_dir)) {
  dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)
}

# Species & Setup ---------------------------------------------------------
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

# Prior & analysis types to check
prior_types <- c("weakly_informative", "informative", "weakly_informative_sensitivity")

# Year indexing
index_keys_path <- here("data/processed/results/stPGOcc/index_keys.rds")
if (file.exists(index_keys_path)) {
  index_keys <- readRDS(index_keys_path)
  year_ids <- index_keys$year_ids
} else {
  year_ids <- 2016:2025
}

# Storage for results -----------------------------------------------------
all_res_list    <- list()
psi_full_list   <- list()
trend_full_list <- list()
ppc_results     <- list()

# Loop Over Species and Prior Types ---------------------------------------
for (spp in possible_bats) {
  for (prior in prior_types) {
    if (prior == "weakly_informative_sensitivity") {
      filepath <- here("data/processed/results/stPGOcc/sensitivity/fits")
      possible_filenames <- c(
        paste0(spp, "_stPGOccfit_weakly_informative_sensitivity.rds"),
        paste0(spp, "_stPGOccfit_weakly_informative.rds"),
        paste0(spp, "_stPGOccfit_weakly_informative_spatial.rds")
      )
    } else {
      filepath <- here("data/processed/results/stPGOcc/full/fits")
      possible_filenames <- c(
        paste0(spp, "_stPGOccfit_", prior, "_spatial.rds"),
        paste0(spp, "_stPGOccfit_", prior, ".rds")
      )
    }
    
    fit_path <- NULL
    for (fn in possible_filenames) {
      if (file.exists(file.path(filepath, fn))) {
        fit_path <- file.path(filepath, fn)
        break
      }
    }
    
    if (is.null(fit_path)) {
      cat(sprintf("Skipping %s (%s): file not found in %s\n", spp, prior, filepath))
      next
    }
    
    cat(sprintf("\n=== Processing stPGOcc Diagnostics: %s (%s) ===\n", spp, prior))
    fit <- readRDS(fit_path)
    
    # 1. Posterior Predictive Check (PPC) -----------------------------------
    cat("  Computing PPC (Freeman-Tukey)...\n")
    ppc <- tryCatch({
      ppcOcc(fit, fit.stat = "freeman-tukey", group = 1)
    }, error = function(e) {
      warning(sprintf("PPC failed for %s (%s): %s", spp, prior, e$message))
      NULL
    })
    
    if (!is.null(ppc)) {
      bayes_pval <- mean(ppc$fit.stat.y.rep > ppc$fit.stat.y)
      cat(sprintf("  Bayesian p-value: %.3f\n", bayes_pval))
      
      ppc_df <- tibble(
        fit_y     = as.vector(ppc$fit.stat.y),
        fit_y_rep = as.vector(ppc$fit.stat.y.rep)
      )
      
      ppc_plot <- ggplot(ppc_df, aes(x = fit_y, y = fit_y_rep)) +
        geom_point(alpha = 0.4, color = "#2c3e50", size = 1.5) +
        geom_abline(intercept = 0, slope = 1, color = "#e74c3c", linetype = "dashed", linewidth = 1) +
        labs(
          title = paste0("Posterior Predictive Check (", toupper(spp), ", ", gsub("_", " ", prior), ")"),
          subtitle = sprintf("Freeman-Tukey Fit Statistic | Bayesian p-value = %.3f", bayes_pval),
          x = "Discrepancy (Observed Data)",
          y = "Discrepancy (Replicated Data)"
        ) +
        theme_bw() +
        theme(plot.title = element_text(face = "bold", size = 12))
      
      ppc_filename <- paste0("ppc_", spp, "_", prior, ".png")
      ggsave(
        filename = file.path(diag_dir, ppc_filename),
        plot = ppc_plot,
        width = 6, height = 5, units = "in", dpi = 300
      )
      
      ppc_results[[paste(spp, prior, sep = "_")]] <- tibble(
        species = spp,
        prior_type = prior,
        bayesian_pval = bayes_pval
      )
    }
    
    # 2. Extract Parameter Summaries ----------------------------------------
    cat("  Summarizing parameters...\n")
    
    summarize_matrix <- function(samples_mat, var_type_label, rhat_vec = NULL, ess_vec = NULL) {
      varnames <- colnames(samples_mat)
      df_list <- list()
      
      for (k in seq_along(varnames)) {
        v <- varnames[k]
        vals <- samples_mat[, k]
        
        df_list[[k]] <- tibble(
          species    = spp,
          prior_type = prior,
          var_type   = var_type_label,
          variable   = v,
          mean       = mean(vals),
          sd         = sd(vals),
          q2.5       = quantile(vals, probs = 0.025, names = FALSE),
          q25        = quantile(vals, probs = 0.25, names = FALSE),
          q50        = quantile(vals, probs = 0.50, names = FALSE),
          q75        = quantile(vals, probs = 0.75, names = FALSE),
          q97.5      = quantile(vals, probs = 0.975, names = FALSE),
          rhat       = if (!is.null(rhat_vec) && v %in% names(rhat_vec)) rhat_vec[[v]] else NA_real_,
          ess_bulk   = if (!is.null(ess_vec) && v %in% names(ess_vec)) ess_vec[[v]] else NA_real_
        )
      }
      bind_rows(df_list)
    }
    
    occ_sum <- summarize_matrix(
      fit$beta.samples, "occurrence",
      rhat_vec = fit$rhat$beta, ess_vec = fit$ESS$beta
    )
    
    det_sum <- summarize_matrix(
      fit$alpha.samples, "detection",
      rhat_vec = fit$rhat$alpha, ess_vec = fit$ESS$alpha
    )
    
    spatial_sum <- summarize_matrix(
      fit$theta.samples, "spatial",
      rhat_vec = fit$rhat$theta, ess_vec = fit$ESS$theta
    )
    
    all_res_list[[paste(spp, prior, sep = "_")]] <- bind_rows(occ_sum, det_sum, spatial_sum)
    
    # 3. Regional Occupancy Trajectories & Trend ----------------------------
    cat("  Computing regional occupancy trajectory & trend...\n")
    psi_sum_draws <- apply(fit$psi.samples, c(1, 3), mean)
    n_years_fit <- ncol(psi_sum_draws)
    actual_years <- min(year_ids):(min(year_ids) + n_years_fit - 1)
    
    colnames(psi_sum_draws) <- as.character(actual_years)
    
    psi_draws_df <- as.data.frame(psi_sum_draws) %>%
      mutate(draw = row_number()) %>%
      pivot_longer(
        cols = -draw,
        names_to = "year",
        values_to = "regional_psi"
      ) %>%
      mutate(year = as.numeric(year))
    
    psi_yr_summary <- psi_draws_df %>%
      group_by(year) %>%
      summarise(
        mean_regional = mean(regional_psi),
        sd_regional   = sd(regional_psi),
        q2.5          = quantile(regional_psi, probs = 0.025, names = FALSE),
        q25           = quantile(regional_psi, probs = 0.25, names = FALSE),
        q50           = quantile(regional_psi, probs = 0.50, names = FALSE),
        q75           = quantile(regional_psi, probs = 0.75, names = FALSE),
        q97.5         = quantile(regional_psi, probs = 0.975, names = FALSE),
        .groups       = "drop"
      ) %>%
      mutate(species = spp, prior_type = prior)
    
    psi_full_list[[paste(spp, prior, sep = "_")]] <- psi_yr_summary
    
    # Trend calculation per draw
    trend_draws <- psi_draws_df %>%
      group_by(draw) %>%
      summarise(
        beta_trend = coef(lm(regional_psi ~ year))[2],
        .groups = "drop"
      )
    
    trend_summary <- trend_draws %>%
      summarise(
        mean_trend = mean(beta_trend),
        sd_trend   = sd(beta_trend),
        q2.5       = quantile(beta_trend, probs = 0.025, names = FALSE),
        q25        = quantile(beta_trend, probs = 0.25, names = FALSE),
        q50        = quantile(beta_trend, probs = 0.50, names = FALSE),
        q75        = quantile(beta_trend, probs = 0.75, names = FALSE),
        q97.5      = quantile(beta_trend, probs = 0.975, names = FALSE),
        .groups    = "drop"
      ) %>%
      mutate(
        species = spp,
        prior_type = prior,
        trend_status = case_when(
          q2.5 > 0 ~ "Increasing",
          q97.5 < 0 ~ "Decreasing",
          TRUE ~ "Uncertain"
        )
      )
    
    trend_full_list[[paste(spp, prior, sep = "_")]] <- trend_summary
    
    rm(fit)
    gc()
  }
}

# Bind Consolidated Results -----------------------------------------------
all_res_df     <- bind_rows(all_res_list)
psi_full_df    <- bind_rows(psi_full_list)
trend_full_df  <- bind_rows(trend_full_list)
ppc_summary_df <- bind_rows(ppc_results)

if (nrow(all_res_df) == 0) {
  stop("No stPGOcc model fits were found to summarize. Please run 07_dynoccMod_stPGOcc.R first.")
}

# Save Summary Tables -----------------------------------------------------
saveRDS(all_res_df, file.path(diag_dir, "all_res_stPGOcc.rds"))
saveRDS(psi_full_df, file.path(diag_dir, "psi_full_stPGOcc.rds"))
saveRDS(trend_full_df, file.path(diag_dir, "trend_full_stPGOcc.rds"))
write_csv(all_res_df, file.path(diag_dir, "stPGOcc_diagnostics_summary.csv"))
write_csv(trend_full_df, file.path(diag_dir, "stPGOcc_trend_summary.csv"))

if (nrow(ppc_summary_df) > 0) {
  write_csv(ppc_summary_df, file.path(diag_dir, "stPGOcc_ppc_summary.csv"))
}

cat("\nSummary tables saved to:", diag_dir, "\n")

# Multi-Species Visualization Figures -------------------------------------

## 1. Occurrence Parameters Summary ---------------------------------------
cat("Generating occurrence parameter summary plot...\n")
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
      levels = c("(Intercept)", "log_fc", "precip", "dem_max", "log_cliff", "year_eff"),
      labels = c("Intercept", "Log(Forest Cover %)", "Precip (mm)", "Elevation (m)", "Log(Cliff Cover %)", "Year Effect")
    )
  )

ggplot(occ_plot_data, aes(x = variable_label, y = mean, group = prior_label, color = prior_label)) +
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0, position = position_dodge(width = 0.25)) +
  geom_errorbar(aes(ymin = q25, ymax = q75), width = 0, position = position_dodge(width = 0.25), linewidth = 1.1) +
  geom_point(position = position_dodge(width = 0.25), size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~species, scales = "free_y", ncol = 2) +
  labs(
    title = "stPGOcc: Occurrence Model Covariate Comparison",
    x = "Covariate",
    y = "Mean (50% and 95% Credible Interval)",
    color = "Prior"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave(file.path(diag_dir, "occ_pars_summary.png"), width = 14, height = 12, units = "in", dpi = 300)

## 2. Detection Parameters Summary ----------------------------------------
cat("Generating detection parameter summary plot...\n")
det_plot_data <- all_res_df %>%
  filter(var_type == "detection") %>%
  mutate(
    prior_label = case_when(
      prior_type == "informative" ~ "Informative",
      prior_type == "weakly_informative" ~ "Weakly Informative",
      TRUE ~ "Weakly Informative (Sensitivity)"
    ),
    variable_label = factor(
      variable,
      levels = c("(Intercept)", "clutter_percent", "tmin", "dayl", "water_ind"),
      labels = c("Intercept", "Clutter %", "Nightly Min Temp", "Day Length", "Waterbody")
    )
  )

ggplot(det_plot_data, aes(x = variable_label, y = mean, group = prior_label, color = prior_label)) +
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0, position = position_dodge(width = 0.25)) +
  geom_errorbar(aes(ymin = q25, ymax = q75), width = 0, position = position_dodge(width = 0.25), linewidth = 1.1) +
  geom_point(position = position_dodge(width = 0.25), size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~species, scales = "free_y", ncol = 2) +
  labs(
    title = "stPGOcc: Detection Model Covariate Comparison",
    x = "Covariate",
    y = "Mean (50% and 95% Credible Interval)",
    color = "Prior"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave(file.path(diag_dir, "det_pars_summary.png"), width = 14, height = 12, units = "in", dpi = 300)

## 3. Spatial & Temporal Hyperparameter Summary ---------------------------
cat("Generating spatial hyperparameter summary plot...\n")
spatial_plot_data <- all_res_df %>%
  filter(var_type == "spatial") %>%
  mutate(
    prior_label = case_when(
      prior_type == "informative" ~ "Informative",
      prior_type == "weakly_informative" ~ "Weakly Informative",
      TRUE ~ "Weakly Informative (Sensitivity)"
    ),
    variable_label = factor(
      variable,
      levels = c("sigma.sq", "phi", "sigma.sq.t", "rho"),
      labels = c("Spatial Var (sigma^2)", "Spatial Decay (phi)", "Temporal Var (sigma_t^2)", "Temporal AR1 (rho)")
    )
  )

ggplot(spatial_plot_data, aes(x = variable_label, y = mean, group = prior_label, color = prior_label)) +
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0, position = position_dodge(width = 0.25)) +
  geom_errorbar(aes(ymin = q25, ymax = q75), width = 0, position = position_dodge(width = 0.25), linewidth = 1.1) +
  geom_point(position = position_dodge(width = 0.25), size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~species, scales = "free_y", ncol = 2) +
  labs(
    title = "stPGOcc: Spatial and Temporal Hyperparameter Comparison",
    x = "Hyperparameter",
    y = "Mean (50% and 95% Credible Interval)",
    color = "Prior"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave(file.path(diag_dir, "spatial_pars_summary.png"), width = 14, height = 12, units = "in", dpi = 300)

## 4. Occupancy Probability Trajectories over Time -----------------------
cat("Generating regional occupancy trajectory plots...\n")
psi_trend_merged <- psi_full_df %>%
  filter(prior_type != "weakly_informative_sensitivity") %>%
  left_join(
    trend_full_df %>% filter(prior_type != "weakly_informative_sensitivity") %>% select(species, prior_type, trend_status, mean_trend),
    by = c("species", "prior_type")
  )

for (p_type in c("weakly_informative", "informative")) {
  p_label <- if_else(p_type == "informative", "Informative", "Weakly Informative")
  
  psi_p_data <- psi_trend_merged %>% filter(prior_type == p_type)
  
  if (nrow(psi_p_data) > 0) {
    plt <- ggplot(psi_p_data, aes(x = factor(year), y = mean_regional, group = species)) +
      geom_ribbon(aes(ymin = q2.5, ymax = q97.5, fill = trend_status), alpha = 0.25) +
      geom_line(aes(color = trend_status), linewidth = 0.8) +
      geom_point(aes(color = trend_status), size = 1.5) +
      facet_wrap(~species, ncol = 2) +
      scale_y_continuous(limits = c(0, 1)) +
      scale_color_manual(values = c("Increasing" = "#27ae60", "Decreasing" = "#c0392b", "Uncertain" = "#7f8c8d")) +
      scale_fill_manual(values = c("Increasing" = "#27ae60", "Decreasing" = "#c0392b", "Uncertain" = "#7f8c8d")) +
      labs(
        title = paste0("stPGOcc: Regional Occupancy Trajectories (", p_label, " Priors)"),
        x = "Year",
        y = "Regional Occupancy Probability (95% CI)",
        color = "Trend Status",
        fill = "Trend Status"
      ) +
      theme_bw() +
      theme(legend.position = "bottom")
    
    ggsave(
      file.path(diag_dir, paste0("psi_trend_", p_type, ".png")),
      plot = plt, width = 14, height = 12, units = "in", dpi = 300
    )
  }
}

## 5. Trend Estimate Comparison Across Prior Types ------------------------
cat("Generating trend estimate comparison plot...\n")
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
    title = "stPGOcc: Trend Estimate Comparison Across Prior Types",
    x = "Species",
    y = expression(hat(beta)[trend] ~ "(Annual Occupancy Slope)"),
    color = "Prior"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(file.path(diag_dir, "beta_trend_compare.png"), plot = trend_compare_plot, width = 12, height = 7, units = "in", dpi = 300)

cat("\nAll stPGOcc diagnostics and figures completed successfully.\n")
