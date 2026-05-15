#!/bin/bash

# Check if both required arguments were provided 
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Usage: $0 <input_db_directory> <output_table_directory>" 
    echo "Example: $0 DataRaw/database DataRaw/tables"
    exit 1
fi

# Assign arguments to descriptive variables 
DB_DIR="$1"
OUT_DIR="$2"

# Define filenames using the input directory 
SRC="$DB_DIR/PNW_BatHub_Database_20251004.accdb"
SRC_CALLS="$DB_DIR/PNW_BatHub_Database_AcousticOutput_20251004.accdb"
SRC_CALLS2="$DB_DIR/PNW_BatHub_Database_AcousticOutput_2024_to_Present_20251004.accdb"

# Safety check: Does the input directory exist? [cite: 5]
if [ ! -d "$DB_DIR" ]; then
    echo "Error: Input directory '$DB_DIR' not found."
    exit 1
fi

# Ensure the output folder exists [cite: 3]
mkdir -p "$OUT_DIR"

echo "Extracting Bat Hub data from $DB_DIR to $OUT_DIR..." 

# Export selected tables [cite: 2]
mdb-export "$SRC" tblDeployment    > "$OUT_DIR/tblDeployment.csv"
mdb-export "$SRC" tblPointLocation > "$OUT_DIR/tblPointLocation.csv"
mdb-export "$SRC" tblSite          > "$OUT_DIR/tblSite.csv"
mdb-export "$SRC" tluClutterType   > "$OUT_DIR/tluClutterType.csv"
mdb-export "$SRC" tluWaterBodyType > "$OUT_DIR/tluWaterBodyType.csv"
mdb-export "$SRC_CALLS" tblDeploymentDetection7 > "$OUT_DIR/calls_to_2024.csv"
mdb-export "$SRC_CALLS2" tblDeploymentDetection8 > "$OUT_DIR/calls_from_2024.csv"

echo "Export complete."