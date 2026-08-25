# Covariate Data (`data/raw/covariates/`)

This directory contains spatial and climate covariate layers from various sources used in detection and occurrence model preparation.

## Subdirectories

| Directory | Description |
| :--- | :--- |
| `daymet/` | Daily surface weather data downloaded via `src/04_get_daymet.py` (`daymet_sites.csv`, `daymet_output.csv` containing minimum temperature, day length, and precipitation). |
| `LandFire/` | Existing Vegetation Type (EVT) raster/spatial data for ID, OR, and WA used to identify cliff and canyon habitats. |
| `NABat_grid_covariates/` | Pre-aggregated NABat CONUS 10km grid-level spatial covariates (forest cover, karst, mean elevation, mean precipitation). |

## Files Overview

| File | Description |
| :--- | :--- |
| `README.md` | Overview of raw covariate datasets and directory structure. |

