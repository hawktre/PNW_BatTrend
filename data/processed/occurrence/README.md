# Occurrence Data (data/processed/occurrence/)

Grid-level spatial covariates and sampling history.

## Files Overview

| File | Description | Input Files | Output Files |
| :--- | :--- | :--- | :--- |
| `batgrid_covars.shp` | NABat 10km grid with aggregated spatial covariates (forest, karst, cliff/canyon, etc.). | `data/raw/covariates/`, `data/raw/batgrid/` | N/A (Produced by `02_CompileSpatialCovariates.R`) |
| `nw_grid.rds` | Final grid dataset with sampling history, ready for occupancy modeling. | `batgrid_covars.shp`, `nw_nights.rds` | N/A (Produced by `05_occurrence_modprep.R`) |
