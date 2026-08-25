# Detection Data (`data/processed/detections/`)

Cleaned and formatted detection data for multi-species dynamic occupancy modeling.

## Files Overview

| File | Description | Input Files | Output Files |
| :--- | :--- | :--- | :--- |
| `deployments_to2025.rds` | Cleaned acoustic deployment records with point coordinates, site IDs, clutter percent, and water body types across 2016–2025. | `data/raw/tables/*.csv` | N/A (Produced by [`src/01_DeploymentsCleaning.R`](file:///src/01_DeploymentsCleaning.R)) |
| `detection_histories.rds` | Wide-format detection data for all species across surveyed deployments. | `data/raw/tables/calls_*.csv`, `deployments_to2025.rds` | N/A (Produced by [`src/03_DetectionsCleaning.R`](file:///src/03_DetectionsCleaning.R)) |
| `nw_nights.rds` | Final observation-level detection dataset joined with Daymet climate covariates (tmin, dayl) and grid attributes. | `detection_histories.rds`, `daymet_output.csv` | N/A (Produced by [`src/05_detections_modprep.R`](file:///src/05_detections_modprep.R)) |
| `sites_missing_covars.csv` | List of sites dropped due to missing habitat / clutter data. | `data/raw/tables/` | Produced by [`src/01_DeploymentsCleaning.R`](file:///src/01_DeploymentsCleaning.R) |
| `missing_clutterpercent.png` | Diagnostic map of sites dropped due to missing clutter covariate values. | `data/raw/tables/` | Produced by [`src/01_DeploymentsCleaning.R`](file:///src/01_DeploymentsCleaning.R) |

