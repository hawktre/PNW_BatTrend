# Detection Data (data/processed/detections/)

Cleaned and formatted detection data for occupancy modeling.

## Files Overview

| File | Description | Input Files | Output Files |
| :--- | :--- | :--- | :--- |
| `deployments_to2024.rds` | Cleaned deployment records with spatial and habitat attributes. | `data/raw/tables/` | N/A (Produced by `01_DeploymentsCleaning.R`) |
| `detection_histories.rds` | Wide-format detection data for all species. | `data/raw/tables/calls_*.csv` | N/A (Produced by `03_DetectionsCleaning.R`) |
| `nw_nights.rds` | Final detection dataset joined with Daymet climate covariates. | `detection_histories.rds`, `daymet_output.csv` | N/A (Produced by `04_detections_modprep.R`) |
| `sites_missing_covars.csv` | Sites dropped due to missing habitat data. | N/A | Produced by `01_DeploymentsCleaning.R` |
| `missing_clutterpercent.png` | Diagnostic map of sites with missing covariates. | N/A | Produced by `01_DeploymentsCleaning.R` |
