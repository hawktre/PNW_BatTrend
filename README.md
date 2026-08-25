# PNW_BatTrend

A repository for analyzing multi-species bat population trends in the Pacific Northwest (PNW) across the 2016–2025 monitoring period using Bayesian dynamic occupancy modeling implemented in Stan (`cmdstanr`).

> [!NOTE]
> This project codebase and analytical workflow are developed and maintained with the assistance of Google Gemini.

---

## Environment Setup

This project uses both **R** and **Python**. Environment management and reproducibility are handled by `renv` for R packages and `venv` for Python dependencies.

### R Environment (`renv`)

The R environment contains all necessary packages for data wrangling, spatial processing, and Bayesian dynamic modeling (e.g., `tidyverse`, `sf`, `terra`, `cmdstanr`, `bayesplot`, `posterior`, `tidybayes`).

#### 1. First-time Setup / Restoration
When cloning or opening the project in RStudio or an R session, `renv` automatically activates via [`.Rprofile`](file:///.Rprofile). To restore the exact package versions recorded in [`renv.lock`](file:///renv.lock), run:
```r
renv::restore()
```

#### 2. CmdStan Toolchain Setup
To run the Stan models via `cmdstanr`, ensure the CmdStan C++ toolchain and CmdStan backend are installed:
```r
# Verify or install CmdStan backend
cmdstanr::check_cmdstan_toolchain()
cmdstanr::install_cmdstan()
```

#### 3. Package Management Best Practices
- **Check Environment Health:** Run `renv::status()` to verify that installed packages match `renv.lock`.
- **Install New Packages:** Run `renv::install("package_name")`.
- **Update Lockfile:** After modifying package dependencies, capture updates by running `renv::snapshot()`.

---

### Python Environment (`venv`)

The Python environment is used for downloading daily climate covariates (Daymet) via `pydaymet`.

1. **Create Virtual Environment:**
   ```bash
   python3 -m venv .venv
   ```
2. **Activate Environment:**
   - **macOS / Linux:** `source .venv/bin/activate`
   - **Windows:** `.venv\Scripts\activate`
3. **Install Dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

---

## Pipeline Execution

### Automated Pipeline (`run_pipeline.sh`)
To execute the end-to-end data processing and diagnostic pipeline sequentially, run:
```bash
./run_pipeline.sh
```

### Manual Execution Order
Individual scripts in [`src/`](file:///src/README.md) can be run sequentially as follows:

1. [`src/00_ExportTables.sh`](file:///src/00_ExportTables.sh): Export Access database tables (`data/raw/database/`) to CSV (`data/raw/tables/`).
2. [`src/01_DeploymentsCleaning.R`](file:///src/01_DeploymentsCleaning.R): Clean detector deployment metadata and calculate site-level characteristics.
3. [`src/02_CompileSpatialCovariates.R`](file:///src/02_CompileSpatialCovariates.R): Compile NABat 10km grid and LandFire spatial layers.
4. [`src/03_DetectionsCleaning.R`](file:///src/03_DetectionsCleaning.R): Filter and clean acoustic call records across monitoring years.
5. [`src/04_get_daymet.py`](file:///src/04_get_daymet.py): Download Daymet surface weather and climate data for surveyed locations.
6. [`src/05_detections_modprep.R`](file:///src/05_detections_modprep.R): Merge detections with climate covariates into nightly detection histories.
7. [`src/06_occurrence_modprep.R`](file:///src/06_occurrence_modprep.R): Assemble regional grid-level occurrence covariates.
8. [`src/07_dynocc_autolog.R`](file:///src/07_dynocc_autolog.R): Fit multi-species dynamic occupancy autologistic models in Stan (`dynocc_autologistic.stan`).
9. [`src/08_dynoccSensitivity.R`](file:///src/08_dynoccSensitivity.R): Fit sensitivity models excluding Idaho data.
10. [`src/09_stanDiagnostics.R`](file:///src/09_stanDiagnostics.R): Compute MCMC diagnostics, posterior predictive checks, Moran's I correlograms, parameter summaries, and generate individual species diagnostic figures.
11. [`src/10_stanMaps.R`](file:///src/10_stanMaps.R): Generate spatial prediction maps across the Pacific Northwest grid.

---

## Directory Overview

| Directory | Description |
| :--- | :--- |
| [`src/`](src/README.md) | Shell, R, Python, and Stan scripts for data processing, modeling, and diagnostics. |
| [`data/raw/`](data/raw/ReadMe.md) | Original, unmodified data files (MS Access databases, spatial layers, trait matrices). |
| [`data/processed/`](data/processed/README.md) | Cleaned, formatted datasets, design matrices, and Stan model results. |
| [`output/`](output/) | Manuscript reports, Quarto documents (`model_results.qmd`), and publication-ready figures. |
