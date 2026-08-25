# Source Scripts (`src/`)

This directory contains all data extraction, cleaning, modeling, diagnostic, and visualization scripts for the Pacific Northwest Bat Trend Analysis.

## Active Pipeline Scripts

| File | Function | Input Files | Output Files |
| :--- | :--- | :--- | :--- |
| [`00_ExportTables.sh`](file:///src/00_ExportTables.sh) | Extracts database tables from MS Access (`.accdb`) to CSV using `mdbtools`. | `data/raw/database/*.accdb` | `data/raw/tables/*.csv` |
| [`01_DeploymentsCleaning.R`](file:///src/01_DeploymentsCleaning.R) | Cleans acoustic detector deployment records, extracts metadata, and exports surveyed locations. | `data/raw/tables/` (`tblDeployment.csv`, `tblPointLocation.csv`, `tblSite.csv`, `tluClutterType.csv`, `tluWaterBodyType.csv`) | `data/processed/detections/deployments_to2025.rds`, `data/raw/covariates/daymet/daymet_sites.csv` |
| [`02_CompileSpatialCovariates.R`](file:///src/02_CompileSpatialCovariates.R) | Extracts and aggregates spatial covariates (forest cover, karst, cliff/canyon, elevation, precipitation) for the NABat 10km grid. | `data/raw/covariates/`, `data/raw/batgrid/` | `data/processed/occurrence/batgrid_covars.shp` |
| [`03_DetectionsCleaning.R`](file:///src/03_DetectionsCleaning.R) | Cleans and harmonizes species call identifications across all acoustic output CSVs (2016–2025). | `data/processed/detections/deployments_to2025.rds`, `data/raw/tables/calls_*.csv` | `data/processed/detections/detection_histories.rds` |
| [`04_get_daymet.py`](file:///src/04_get_daymet.py) | Downloads daily Daymet weather data (tmin, dayl, prcp) for surveyed sites via `pydaymet`. | `data/raw/covariates/daymet/daymet_sites.csv` | `data/raw/covariates/daymet/daymet_output.csv` |
| [`05_detections_modprep.R`](file:///src/05_detections_modprep.R) | Merges cleaned species detections with climate covariates and formats nightly survey replicates. | `data/processed/detections/detection_histories.rds`, `data/raw/covariates/daymet/daymet_output.csv` | `data/processed/detections/nw_nights.rds` |
| [`06_occurrence_modprep.R`](file:///src/06_occurrence_modprep.R) | Finalizes grid-level occurrence covariate tables and joins historical survey status. | `data/processed/occurrence/batgrid_covars.shp`, `data/processed/detections/nw_nights.rds` | `data/processed/occurrence/nw_grid.rds` |
| [`07_dynocc_autolog.R`](file:///src/07_dynocc_autolog.R) | Formats Stan data matrices and fits multi-species autologistic dynamic occupancy models across priors. | `data/processed/occurrence/nw_grid.rds`, `data/processed/detections/nw_nights.rds` | `data/processed/results/stan/stan_data_*.rds`, `design_matrices.rds`, `index_keys.rds`, fits in `data/processed/results/stan/fits/` |
| [`08_dynoccSensitivity.R`](file:///src/08_dynoccSensitivity.R) | Runs sensitivity analysis by fitting models excluding Idaho survey data. | `data/processed/occurrence/nw_grid.rds`, `data/processed/detections/nw_nights.rds` | Sensitivity design matrices and fits in `data/processed/results/stan/fits/` |
| [`09_stanDiagnostics.R`](file:///src/09_stanDiagnostics.R) | Evaluates MCMC convergence, computes posterior predictive checks, calculates Moran's I correlograms (10–50 km), estimates regional occupancy trajectories/slopes, and outputs multi-species and individual species figures. | Stan fits in `data/processed/results/stan/fits/` | Diagnostic tables (`all_res_stan.rds`, `trend_full_stan.rds`, etc.) and plots in `data/processed/results/stan/diagnostics/` |
| [`10_stanMaps.R`](file:///src/10_stanMaps.R) | Projects posterior parameter distributions across the Northwest grid to create regional occupancy prediction maps and uncertainty plots. | Stan fits in `data/processed/results/stan/fits/`, `nw_grid.rds` | Spatial prediction maps in `data/processed/results/stan/maps/` |
| [`dynocc_autologistic.stan`](file:///src/dynocc_autologistic.stan) | Stan model file implementing a first-order autologistic dynamic occupancy model with time-varying survival ($\phi_t$), colonization ($\gamma_t$), and logit-linear detection and initial occurrence covariates. | N/A (Compiled by `cmdstanr`) | N/A |

---

## Subdirectories

- **`archive/`**: Contains legacy exploratory scripts from earlier iterations (e.g., `spOccupancy`, `jagsUI`, `stPGOcc`, `TV_dynocc_forward.stan`).

