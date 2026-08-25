# Processed Data (`data/processed/`)

This directory contains intermediate and finalized data files that have been cleaned, formatted, and prepared for dynamic occupancy modeling.

## Subdirectories

| Directory | Description |
| :--- | :--- |
| [`detections/`](detections/README.md) | Cleaned deployment records, detection histories, and merged nightly observation datasets. |
| [`occurrence/`](occurrence/README.md) | Regional NABat 10km grid shapefiles and joined occurrence covariates. |
| `results/` | Stan model input data, design matrices, MCMC fit objects (`results/stan/fits/`), diagnostic summaries and figures (`results/stan/diagnostics/`), and prediction maps (`results/stan/maps/`). |

