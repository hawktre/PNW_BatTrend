# Source Scripts (src/)

This directory contains all scripts required to clean raw data, compile covariates, and run occupancy models.

## Scripts Overview

| File | Function | Input Files | Output Files |
| :--- | :--- | :--- | :--- |
| `00_ExportTables.sh` | Exports tables from MS Access databases to CSV. | `.accdb` files in `data/raw/database/` | `.csv` files in `data/raw/tables/` |
| `01_DeploymentsCleaning.R` | Cleans deployment records and calculates site-level metrics. | `data/raw/tables/` (various .csv) | `data/processed/detections/deployments_to2024.rds`, `data/raw/covariates/daymet/daymet_sites.csv` |
| `01a_get_daymet.py` | Downloads Daymet climate data for surveyed sites. | `data/raw/covariates/daymet/daymet_sites.csv` | `data/raw/covariates/daymet/daymet_output.csv` |
| `02_CompileSpatialCovariates.R` | Compiles and aggregates NABat grid and LandFire spatial covariates. | `data/raw/covariates/`, `data/raw/batgrid/` | `data/processed/occurrence/batgrid_covars.shp` |
| `03_DetectionsCleaning.R` | Formats species detection data for occupancy modeling. | `data/processed/detections/deployments_to2024.rds`, `data/raw/tables/calls_*.csv` | `data/processed/detections/detection_histories.rds` |
| `04_detections_modprep.R` | Joins detections with climate data and prepares temporal replicates. | `data/processed/detections/detection_histories.rds`, `data/raw/covariates/daymet/daymet_output.csv` | `data/processed/detections/nw_nights.rds` |
| `05_occurrence_modprep.R` | Finalizes grid-level covariate preparation for modeling. | `data/processed/occurrence/batgrid_covars.shp`, `data/processed/detections/nw_nights.rds` | `data/processed/occurrence/nw_grid.rds` |
| `07_jagsMod.R` | Runs occupancy models using JAGS. | `data/processed/` (.rds files) | Model results in `data/processed/results/` |
| `08_tPGOcc.R` | Runs multi-season occupancy models using `spOccupancy`. | `data/processed/` (.rds files) | Model results in `data/processed/results/` |
| `09_stPGocc.R` | Runs spatial-temporal occupancy models using `spOccupancy`. | `data/processed/` (.rds files) | Model results in `data/processed/results/` |
| `10_MethodComparisons.R` | Compares results from different modeling approaches. | Model outputs | Comparison plots and tables |
| `11_summarise_spOcc.R` | Summarizes and visualizes `spOccupancy` model results. | Model outputs | Summaries and maps |
| `occ_model_royle.jags` | JAGS model specification file. | N/A | N/A |
