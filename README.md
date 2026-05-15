# PNW_BatTrend

A repository for analyzing bat population trends in the Pacific Northwest using occupancy modeling.

## Environment Setup

This project uses both R and Python. Environment management is handled by `renv` for R and `venv` for Python.

### R Environment (renv)
The R environment contains all necessary packages for spatial analysis, data cleaning, and Bayesian modeling (e.g., `tidyverse`, `sf`, `terra`, `rjags`, `spOccupancy`).

1.  Open the project in RStudio or your preferred R IDE.
2.  The `renv` environment should automatically initialize via `.Rprofile`.
3.  Run the following command to restore the project library:
    ```r
    renv::restore()
    ```

### Python Environment (venv)
The Python environment is primarily used for downloading climate data (Daymet) via `pydaymet`.

1.  Create the virtual environment (if not already present):
    ```bash
    python3 -m venv .venv
    ```
2.  Activate the environment:
    - **macOS/Linux:** `source .venv/bin/activate`
    - **Windows:** `.venv\Scripts\activate`
3.  Install dependencies:
    ```bash
    pip install -r requirements.txt
    ```

### Automated Pipeline
To execute the entire analysis pipeline from data extraction to figure generation, run:
```bash
./run_pipeline.sh
```
This script will run each step in sequence, providing progress updates and stopping if an error occurs.

---

## Directory Overview

| Directory | Description |
| :--- | :--- |
| [src/](src/README.md) | Shell, R, and Python scripts for data processing and modeling. |
| [data/raw/](data/raw/ReadMe.md) | Original, unmodified data files including databases and spatial layers. |
| [data/processed/](data/processed/README.md) | Cleaned and formatted data files ready for analysis. |
