## ---------------------------
## Purpose of script: Read JAGS model fits, compute convergence diagnostics
##                    (Rhat, ESS) for all and focal parameters, and generate
##                    diagnostic plots (traceplots, ACF, posteriors) saved as PNGs.
##
## Author: Antigravity AI Pair Programmer
## ---------------------------

# Load Packages -----------------------------------------------------------
library(tidyverse)
library(here)
library(cmdstanr)
library(bayesplot)
library(posterior)

# Paths and Setup ---------------------------------------------------------
fits_dir <- here("data/processed/results/stan/full/fits/")
diag_dir <- here("data/processed/results/stan/full/diagnostics")

if (!dir.exists(diag_dir)) {
  dir.create(diag_dir, recursive = TRUE)
}

# Clean up any existing PDF files in the diagnostics directory
old_pdfs <- list.files(diag_dir, pattern = "\\.pdf$", full.names = TRUE)
if (length(old_pdfs) > 0) {
  file.remove(old_pdfs)
}

# Load Metadata for Labels -------------------------------------------------
# Load design matrices to map coefficient indices to covariate names
design_matrices_path <- here("data/processed/results/stan/design_matrices.rds")
if (file.exists(design_matrices_path)) {
  mats <- readRDS(design_matrices_path)
  xcov_names_standard <- colnames(mats$xmata)
  xcov_names_cliff <- colnames(mats$xmatc)
  pcov_names <- colnames(mats$design_matrix_det)
} else {
  warning("design_matrices.rds not found. Using default index labels.")
  xcov_names_standard <- NULL
  xcov_names_cliff <- NULL
  pcov_names <- NULL
}

# Load index keys to map time steps to years
index_keys_path <- here("data/processed/results/stan/index_keys.rds")
if (file.exists(index_keys_path)) {
  keys <- readRDS(index_keys_path)
  year_ids <- keys$year_ids
} else {
  warning("index_keys.rds not found. Using default index labels.")
  year_ids <- NULL
}

# Cliff-associated species list to determine occupancy covariates
cliff_spp <- c("anpa", "euma", "myci", "pahe")

# Helper: Map parameter names to descriptive labels ----------------------
get_clean_label <- function(var_name, species, year_ids, xcov_names_standard, xcov_names_cliff, pcov_names) {
  if (var_name == "alpha01") {
    return("alpha01 (Occupancy Intercept Yr 1)")
  }
  if (var_name == "deviance") {
    return("deviance")
  }
  
  # Determine occupancy covariate list
  is_cliff <- species %in% cliff_spp
  xcov_names <- if (is_cliff) xcov_names_cliff else xcov_names_standard
  
  # Check for alphas[index] using string matching to avoid regex escaping issues
  if (startsWith(var_name, "alphas[")) {
    idx_str <- gsub("[^0-9]", "", var_name)
    idx <- as.integer(idx_str)
    cov_name <- if (!is.null(xcov_names) && idx <= length(xcov_names)) xcov_names[idx] else paste0("Cov", idx)
    return(sprintf("alphas[%d] (%s)", idx, cov_name))
  }
  
  # Check for betas[index]
  if (startsWith(var_name, "betas[")) {
    idx_str <- gsub("[^0-9]", "", var_name)
    idx <- as.integer(idx_str)
    cov_name <- if (!is.null(pcov_names) && idx <= length(pcov_names)) pcov_names[idx] else paste0("Cov", idx)
    return(sprintf("betas[%d] (%s)", idx, cov_name))
  }
  
  # Check for phi[index]
  if (startsWith(var_name, "phi[")) {
    idx_str <- gsub("[^0-9]", "", var_name)
    idx <- as.integer(idx_str)
    if (!is.null(year_ids) && idx < length(year_ids)) {
      yr1 <- year_ids[idx]
      yr2 <- year_ids[idx + 1]
      return(sprintf("phi[%d] (%d -> %d)", idx, yr1, yr2))
    } else {
      return(sprintf("phi[%d] (Yr%d -> Yr%d)", idx, idx, idx + 1))
    }
  }
  
  # Check for gamma[index]
  if (startsWith(var_name, "gamma[")) {
    idx_str <- gsub("[^0-9]", "", var_name)
    idx <- as.integer(idx_str)
    if (!is.null(year_ids) && idx < length(year_ids)) {
      yr1 <- year_ids[idx]
      yr2 <- year_ids[idx + 1]
      return(sprintf("gamma[%d] (%d -> %d)", idx, yr1, yr2))
    } else {
      return(sprintf("gamma[%d] (Yr%d -> Yr%d)", idx, idx, idx + 1))
    }
  }
  
  return(var_name)
}

# Helper: Create safe filenames from parameter names ---------------------
get_safe_filename <- function(var_name) {
  # Replace brackets with underscores and remove closing brackets
  # e.g., "alphas[1]" -> "alphas_1"
  safe_name <- gsub("\\[", "_", var_name)
  safe_name <- gsub("\\]", "", safe_name)
  return(safe_name)
}

# Find All Fit Files ------------------------------------------------------
fit_files <- list.files(fits_dir, pattern = ".*_stanfit_.*\\.rds$", full.names = TRUE)

if (length(fit_files) == 0) {
  stop(sprintf("No JAGS fit files found in: %s", fits_dir))
}

cat("Found", length(fit_files), "fit files to process.\n\n")

# Summary data collector
diagnostics_list <- list()

# Process Each Fit File ---------------------------------------------------

