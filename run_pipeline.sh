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

# Step 1a: Get Daymet (Python)
# Use the virtual environment python if it exists, otherwise fall back to system python3
PYTHON_CMD="python3"
if [ -f ".venv/bin/python" ]; then
    PYTHON_CMD=".venv/bin/python"
elif [ -f ".venv/bin/python3" ]; then
    PYTHON_CMD=".venv/bin/python3"
fi

run_step "src/01a_get_daymet.py" \
    "Downloading Daymet climate data for surveyed sites" \
    "$PYTHON_CMD src/01a_get_daymet.py"

# Step 2: Compile Spatial Covariates
run_step "src/02_CompileSpatialCovariates.R" \
    "Compiling and aggregating NABat grid and LandFire spatial covariates" \
    "Rscript src/02_CompileSpatialCovariates.R"

# Step 3: Clean Detections
run_step "src/03_DetectionsCleaning.R" \
    "Formatting species detection data for occupancy modeling" \
    "Rscript src/03_DetectionsCleaning.R"

# Step 4: Detection Model Prep
run_step "src/04_detections_modprep.R" \
    "Joining detections with climate data and preparing temporal replicates" \
    "Rscript src/04_detections_modprep.R"

# Step 5: Occurrence Model Prep
run_step "src/05_occurrence_modprep.R" \
    "Finalizing grid-level covariate preparation for modeling" \
    "Rscript src/05_occurrence_modprep.R"

# --- Modeling Phase ---
# These steps are computationally intensive and may take significant time.

# Step 7: JAGS Modeling
run_step "src/07_jagsMod.R" \
    "Running occupancy models using JAGS" \
    "Rscript src/07_jagsMod.R"

# Step 8: tPGOcc Modeling
run_step "src/08_tPGOcc.R" \
    "Running multi-season occupancy models using spOccupancy (non-spatial)" \
    "Rscript src/08_tPGOcc.R"

# Step 9: stPGOcc Modeling
run_step "src/09_stPGocc.R" \
    "Running spatial-temporal occupancy models using spOccupancy" \
    "Rscript src/09_stPGocc.R"

# --- Summarization Phase ---

# Step 11: Summarize spOccupancy
run_step "src/10_summarise_spOcc.R" \
    "Summarizing and visualizing spOccupancy model results" \
    "Rscript src/10_summarise_spOcc.R"

# Step 10: Method Comparisons
run_step "src/11_MethodComparisons.R" \
    "Comparing results from different modeling approaches and generating final figures" \
    "Rscript src/11_MethodComparisons.R"

echo "===================================================="
echo -e "${GREEN}Pipeline completed successfully!${NC}"
echo "===================================================="
