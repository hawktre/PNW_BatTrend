#!/bin/bash

# Master script to run the PNW Bat Trend Analysis Pipeline sequentially.
# It ensures all data extraction, cleaning, and modeling steps are executed in order.

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Function to run a step and check for success
run_step() {
    local script_path=$1
    local description=$2
    local command=$3

    echo -e "${GREEN}>>> Running: $script_path${NC}"
    echo -e "Function: $description"
    
    # Execute the command
    if eval "$command"; then
        echo -e "${GREEN}✔ Success: $script_path${NC}\n"
    else
        echo -e "${RED}✘ Error: $script_path failed.${NC}"
        echo "Please investigate the error in the output above before continuing."
        exit 1
    fi
}

echo "===================================================="
echo "      PNW Bat Trend Analysis Pipeline"
echo "===================================================="

# Step 0: Export Data
run_step "src/00_ExportTables.sh" \
    "Exporting tables from MS Access databases to CSV" \
    "bash src/00_ExportTables.sh data/raw/database data/raw/tables"

# Step 1: Clean Deployments
run_step "src/01_DeploymentsCleaning.R" \
    "Cleaning deployment records and calculating site metrics" \
    "Rscript src/01_DeploymentsCleaning.R"


# Step 2: Compile Spatial Covariates
run_step "src/02_CompileSpatialCovariates.R" \
    "Compiling and aggregating NABat grid and LandFire spatial covariates" \
    "Rscript src/02_CompileSpatialCovariates.R"

# Step 3: Clean Detections
run_step "src/03_DetectionsCleaning.R" \
    "Formatting species detection data for occupancy modeling" \
    "Rscript src/03_DetectionsCleaning.R"

# Step 4: Get Daymet (Python)
run_step "src/04_get_daymet.py" \
    "Downloading Daymet climate data for surveyed sites via pydaymet" \
    ".venv/bin/python src/04_get_daymet.py"

# Step 5: Detection Model Prep
run_step "src/05_detections_modprep.R" \
    "Joining detections with climate data and preparing temporal replicates" \
    "Rscript src/05_detections_modprep.R"

# Step 6: Occurrence Model Prep
run_step "src/06_occurrence_modprep.R" \
    "Finalizing grid-level covariate preparation for modeling" \
    "Rscript src/06_occurrence_modprep.R"

# --- Modeling & Diagnostics Phase ---
# Stan model fitting is computationally intensive. If fits already exist,
# fitting steps can be run or skipped depending on workflow needs.

# Step 7: Stan Dynamic Occupancy Modeling (Autologistic)
# run_step "src/07_dynocc_autolog.R" \
#     "Fitting dynamic occupancy autologistic models in Stan across species" \
#     "Rscript src/07_dynocc_autolog.R"

# Step 8: Stan Sensitivity Analysis (Excluding Idaho)
# run_step "src/08_dynoccSensitivity.R" \
#     "Fitting sensitivity analysis models in Stan without Idaho data" \
#     "Rscript src/08_dynoccSensitivity.R"

# Step 9: Stan Diagnostics, Summaries, and Diagnostic Plots
if [ -d "data/processed/results/stan/fits" ] && [ "$(find data/processed/results/stan/fits -maxdepth 1 -name '*_stanfit_*.rds' | wc -l)" -gt 0 ]; then
    run_step "src/09_stanDiagnostics.R" \
        "Computing MCMC diagnostics, posterior predictive checks, Moran's I, and trend summaries" \
        "Rscript src/09_stanDiagnostics.R"
else
    echo -e "${GREEN}>>> Skipping: src/09_stanDiagnostics.R (No Stan fit files found in data/processed/results/stan/fits)${NC}\n"
fi

# Step 10: Stan Spatial Prediction Maps
if [ -d "data/processed/results/stan/fits" ] && [ "$(find data/processed/results/stan/fits -maxdepth 1 -name '*_stanfit_*.rds' | wc -l)" -gt 0 ]; then
    run_step "src/10_stanMaps.R" \
        "Generating regional spatial prediction maps of occupancy and uncertainty" \
        "Rscript src/10_stanMaps.R"
else
    echo -e "${GREEN}>>> Skipping: src/10_stanMaps.R (No Stan fit files found in data/processed/results/stan/fits)${NC}\n"
fi

echo "===================================================="
echo -e "${GREEN}Pipeline completed successfully!${NC}"
echo "===================================================="

