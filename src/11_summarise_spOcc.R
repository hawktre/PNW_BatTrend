## ---------------------------
## Purpose of script: Summarise posterior estimates from tPGOcc and stPGOcc
##                    model fits for comparison with JAGS results. Outputs
##                    psi_bar, lambda, and parameter summaries in the same
##                    long-format structure as 01_ParameterSummaries.R.
##                    Run after 05_stPGOcc.R and 06_tPGOcc.R.
##
## Author: Trent VanHawkins
##
## Date Created: 2026-04-05
##
## ---------------------------

options(scipen = 6, digits = 4)

# Load Packages -----------------------------------------------------------
library(tidyverse)
library(here)
library(sf)
library(spOccupancy)

# Load Data ---------------------------------------------------------------
nw_grid <- readRDS(here("data/processed/occurrence/nw_grid.rds"))
index_keys <- readRDS(here("data/processed/results/jags/index_keys.rds"))
site_ids <- index_keys$site_ids
year_ids <- index_keys$year_ids

# Create Site/State Key ---------------------------------------------------
site_state_key <- nw_grid %>%
  select(sample_unit_id, state) %>%
  filter(sample_unit_id %in% site_ids) %>%
  st_drop_geometry() %>%
  arrange(sample_unit_id) %>%
  mutate(site_index = row_number())

# Functions ---------------------------------------------------------------

## Compute psi_bar overall and by state, plus lambda trend
## Works for both tPGOcc and stPGOcc — psi.samples has same structure
compute_psi_bar <- function(fit, site_info) {
  psi_array <- fit$psi.samples # dim(n_samples, n_sites, n_years)

  n_iter <- dim(psi_array)[1]
  n_sites <- dim(psi_array)[2]
  n_years <- dim(psi_array)[3]

  psi_df <- as.data.frame.table(psi_array, responseName = "psi") %>%
    rename(iter = Var1, site = Var2, year = Var3) %>%
    mutate(
      iter = as.integer(iter),
      site = as.integer(site),
      year = as.integer(year)
    ) %>%
    left_join(site_info, by = c("site" = "site_index"))

  ## Overall psi_bar per year
  psi_bar_all_iter <- psi_df %>%
    group_by(year, iter) %>%
    summarize(psi_bar = mean(psi), .groups = "drop")

  psi_bar_all <- psi_bar_all_iter %>%
    group_by(year) %>%
    summarize(
      mean = mean(psi_bar),
      lci = quantile(psi_bar, 0.025),
      q25 = quantile(psi_bar, 0.25),
      q75 = quantile(psi_bar, 0.75),
      uci = quantile(psi_bar, 0.975),
      .groups = "drop"
    )

  ## psi_bar by state per year
  psi_bar_state_iter <- psi_df %>%
    group_by(state, year, iter) %>%
    summarize(psi_bar = mean(psi), .groups = "drop")

  psi_bar_by_state <- psi_bar_state_iter %>%
    group_by(state, year) %>%
    summarize(
      mean = mean(psi_bar),
      lci = quantile(psi_bar, 0.025),
      q25 = quantile(psi_bar, 0.25),
      q75 = quantile(psi_bar, 0.75),
      uci = quantile(psi_bar, 0.975),
      .groups = "drop"
    )

  ## Lambda: ratio of psi_bar in last year vs first year
  lambda_tot <- psi_bar_all_iter %>%
    filter(year %in% c(1, n_years)) %>%
    pivot_wider(
      names_from = year,
      values_from = psi_bar,
      names_prefix = "year_"
    ) %>%
    mutate(lambda = .data[[paste0("year_", n_years)]] / year_1) %>%
    summarize(
      mean = mean(lambda),
      lci = quantile(lambda, 0.025),
      q25 = quantile(lambda, 0.25),
      q75 = quantile(lambda, 0.75),
      uci = quantile(lambda, 0.975)
    )

  ## Lambda by state
  lambda_tot_state <- psi_bar_state_iter %>%
    filter(year %in% c(1, n_years)) %>%
    pivot_wider(
      names_from = year,
      values_from = psi_bar,
      names_prefix = "year_"
    ) %>%
    mutate(lambda = .data[[paste0("year_", n_years)]] / year_1) %>%
    group_by(state) %>%
    summarize(
      mean = mean(lambda),
      lci = quantile(lambda, 0.025),
      q25 = quantile(lambda, 0.25),
      q75 = quantile(lambda, 0.75),
      uci = quantile(lambda, 0.975),
      .groups = "drop"
    )

  list(
    overall = list(psi_bar = psi_bar_all, trend = lambda_tot),
    by_state = list(psi_bar = psi_bar_by_state, trend = lambda_tot_state)
  )
}

## Summarize occupancy and detection coefficients
## Works for both tPGOcc and stPGOcc — pulls directly from fit object
summarize_params <- function(fit) {
  occ_names <- names(fit$ESS$beta)
  det_names <- names(fit$ESS$alpha)

  occ_df <- as_tibble(fit$beta.samples) %>%
    pivot_longer(everything(), names_to = "param", values_to = "value")

  det_df <- as_tibble(fit$alpha.samples) %>%
    pivot_longer(everything(), names_to = "param", values_to = "value")

  summary_df <- bind_rows(occ_df, det_df) %>%
    group_by(param) %>%
    summarize(
      mean = mean(value),
      lci = quantile(value, 0.025),
      q25 = quantile(value, 0.25),
      q75 = quantile(value, 0.75),
      uci = quantile(value, 0.975),
      .groups = "drop"
    )

  diag_df <- tibble(
    param = c(occ_names, det_names),
    ess = c(fit$ESS$beta, fit$ESS$alpha),
    rhat = c(fit$rhat$beta, fit$rhat$alpha)
  )

  left_join(summary_df, diag_df, by = "param")
}

## Summarize spatial parameters from theta.samples
## Only relevant for stPGOcc — returns NULL for tPGOcc
summarize_spatial_params <- function(fit) {
  if (is.null(fit$theta.samples)) {
    return(NULL)
  }

  as_tibble(fit$theta.samples) %>%
    pivot_longer(everything(), names_to = "param", values_to = "value") %>%
    group_by(param) %>%
    summarize(
      mean = mean(value),
      lci = quantile(value, 0.025),
      q25 = quantile(value, 0.25),
      q75 = quantile(value, 0.75),
      uci = quantile(value, 0.975),
      ess = coda::effectiveSize(value),
      .groups = "drop"
    )
}

## Helper to run all summary functions for one species
summarize_species <- function(path, site_info) {
  fit <- readRDS(path)
  spp <- str_split_i(tools::file_path_sans_ext(basename(path)), "_", 1)
  cat("Summarizing:", spp, "\n")

  list(
    psi = compute_psi_bar(fit, site_info = site_info),
    params = summarize_params(fit),
    spatial_params = summarize_spatial_params(fit)
  )
}

## Combine species-level results into long format
combine_summaries <- function(res) {
  psi_bar_overall_list <- list()
  trend_overall_list <- list()
  psi_bar_bystate_list <- list()
  trend_bystate_list <- list()
  params_list <- list()
  spatial_params_list <- list()

  for (species in names(res)) {
    spp_data <- res[[species]]

    if (!is.null(spp_data$psi$overall$psi_bar)) {
      df <- spp_data$psi$overall$psi_bar
      df$species <- species
      psi_bar_overall_list[[species]] <- df
    }

    if (!is.null(spp_data$psi$overall$trend)) {
      df <- spp_data$psi$overall$trend
      df$species <- species
      trend_overall_list[[species]] <- df
    }

    if (!is.null(spp_data$psi$by_state$psi_bar)) {
      df <- spp_data$psi$by_state$psi_bar
      df$species <- species
      psi_bar_bystate_list[[species]] <- df
    }

    if (!is.null(spp_data$psi$by_state$trend)) {
      df <- spp_data$psi$by_state$trend
      df$species <- species
      trend_bystate_list[[species]] <- df
    }

    if (!is.null(spp_data$params)) {
      df <- spp_data$params
      df$species <- species
      params_list[[species]] <- df
    }

    if (!is.null(spp_data$spatial_params)) {
      df <- spp_data$spatial_params
      df$species <- species
      spatial_params_list[[species]] <- df
    }
  }

  list(
    psi_bar_overall = bind_rows(psi_bar_overall_list),
    trend_overall = bind_rows(trend_overall_list),
    psi_bar_bystate = bind_rows(psi_bar_bystate_list),
    trend_bystate = bind_rows(trend_bystate_list),
    params = bind_rows(params_list),
    spatial_params = bind_rows(spatial_params_list)
  )
}

## Format combined summary to match all_res structure from 02_TrendFigs.R
format_all_res <- function(summary, analysis_label) {
  list(
    psi = summary$psi_bar_overall %>%
      mutate(
        year = as.integer(year) + min(year_ids) - 1,
        analysis = analysis_label
      ),
    psi_bystate = summary$psi_bar_bystate %>%
      mutate(
        year = as.integer(year) + min(year_ids) - 1,
        analysis = analysis_label
      ),
    trend = summary$trend_overall %>%
      mutate(analysis = analysis_label),
    trend_bystate = summary$trend_bystate %>%
      mutate(analysis = analysis_label),
    params = summary$params %>%
      mutate(analysis = analysis_label),
    spatial_params = summary$spatial_params %>%
      mutate(analysis = analysis_label)
  )
}

# Summarize tPGOcc Results ------------------------------------------------
cat("\n=== Summarizing tPGOcc fits ===\n")

tpg_paths <- list.files(
  here("data/processed/results/tPGOcc/fits/"),
  full.names = TRUE
)
tpg_spps <- str_split_i(tools::file_path_sans_ext(basename(tpg_paths)), "_", 1)

tpg_res <- map(tpg_paths, summarize_species, site_info = site_state_key)
names(tpg_res) <- tpg_spps

saveRDS(tpg_res, here("data/processed/results/tPGOcc/res_summary.rds"))

tpg_summary <- combine_summaries(tpg_res)
all_res_tpg <- format_all_res(tpg_summary, "tPGOcc")

saveRDS(all_res_tpg, here("data/processed/results/tPGOcc/all_res_tpg.rds"))
cat("tPGOcc results saved.\n")

# Summarize stPGOcc Results -----------------------------------------------
cat("\n=== Summarizing stPGOcc fits ===\n")

stpg_paths <- list.files(
  here("data/processed/results/stPGOcc/fits/"),
  full.names = TRUE
)
stpg_spps <- str_split_i(
  tools::file_path_sans_ext(basename(stpg_paths)),
  "_",
  1
)

stpg_res <- map(stpg_paths, summarize_species, site_info = site_state_key)
names(stpg_res) <- stpg_spps

saveRDS(stpg_res, here("data/processed/results/stPGOcc/res_summary.rds"))

stpg_summary <- combine_summaries(stpg_res)
all_res_stpg <- format_all_res(stpg_summary, "stPGOcc")

saveRDS(all_res_stpg, here("data/processed/results/stPGOcc/all_res_stpg.rds"))
cat("stPGOcc results saved.\n")

# Print Spatial Parameter Summary -----------------------------------------
cat("\n--- stPGOcc Spatial Parameter Summary ---\n")
stpg_summary$spatial_params %>%
  mutate(eff_range = 3 / mean) %>% # effective spatial range in coordinate units
  select(species, param, mean, lci, uci, ess, eff_range) %>%
  arrange(species, param) %>%
  print(n = Inf)

cat("\nDone.\n")
cat("tPGOcc: data/processed/results/tPGOcc/all_res_tpg.rds\n")
cat("stPGOcc: data/processed/results/stPGOcc/all_res_stpg.rds\n")
