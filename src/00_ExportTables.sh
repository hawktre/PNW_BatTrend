#!/bin/bash

# Check if both required arguments were provided 
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Usage: $0 <input_db_directory> <output_table_directory>" 
    echo "Example: $0 data/raw/database data/raw/tables"
    exit 1
fi

# Assign arguments to descriptive variables 
DB_DIR="$1"
OUT_DIR="$2"

# Safety check: Does the input directory exist?
if [ ! -d "$DB_DIR" ]; then
    echo "Error: Input directory '$DB_DIR' not found."
    exit 1
fi

# Ensure the output folder exists
mkdir -p "$OUT_DIR"

echo "Extracting Bat Hub data from $DB_DIR to $OUT_DIR..." 

# Track if we found any database files
found_db=false

# Loop over all .accdb files in the input directory
for db_file in "$DB_DIR"/*.accdb; do
    # Skip if no .accdb files are found (glob yields the literal string)
    [ -e "$db_file" ] || continue
    found_db=true
    
    filename=$(basename "$db_file")
    
    # Decide database type based on filename
    if [[ "$filename" == *AcousticOutput* ]]; then
        # Acoustic output database
        # Extract suffix (e.g. 2016-2023, 2024, Idaho_2025)
        suffix="${filename##*AcousticOutput_}"
        suffix="${suffix%.accdb}"
        
        # Find table starting with tblDeploymentDetection
        tables=$(mdb-tables -1 "$db_file" | grep '^tblDeploymentDetection')
        
        if [ -n "$tables" ]; then
            for table in $tables; do
                echo "Exporting acoustic output table '$table' from '$filename' to 'calls_${suffix}.csv'..."
                mdb-export "$db_file" "$table" > "$OUT_DIR/calls_${suffix}.csv"
            done
        else
            echo "Warning: No table starting with 'tblDeploymentDetection' found in '$filename'."
        fi
    else
        # Normal database containing site metadata
        echo "Exporting site metadata tables from '$filename'..."
        mdb-export "$db_file" tblDeployment    > "$OUT_DIR/tblDeployment.csv"
        mdb-export "$db_file" tblPointLocation > "$OUT_DIR/tblPointLocation.csv"
        mdb-export "$db_file" tblSite          > "$OUT_DIR/tblSite.csv"
        mdb-export "$db_file" tluClutterType   > "$OUT_DIR/tluClutterType.csv"
        mdb-export "$db_file" tluWaterBodyType > "$OUT_DIR/tluWaterBodyType.csv"
    fi
done

if [ "$found_db" = false ]; then
    echo "Warning: No .accdb database files found in '$DB_DIR'."
fi

echo "Export complete."