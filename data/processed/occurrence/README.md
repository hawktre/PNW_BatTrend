# Occurrence Data (`data/processed/occurrence/`)

Grid-level spatial covariates, NABat sampling history, and spatial frameworks.

## Files Overview

| File | Description | Input Files | Output Files |
| :--- | :--- | :--- | :--- |
| `batgrid_covars.shp` | NABat 10km grid shapefile for the Pacific Northwest with aggregated spatial covariates (forest cover %, karst %, cliff/canyon %, mean elevation, mean precipitation). | `data/raw/covariates/`, `data/raw/batgrid/` | N/A (Produced by [`src/02_CompileSpatialCovariates.R`](file:///src/02_CompileSpatialCovariates.R)) |
| `nw_grid.rds` | Final grid dataset joined with annual sampling indicators across 2016–2025, ready for dynamic occupancy modeling. | `batgrid_covars.shp`, `nw_nights.rds` | N/A (Produced by [`src/06_occurrence_modprep.R`](file:///src/06_occurrence_modprep.R)) |

