# Bat Hub Database (`data/raw/database/`)

This directory contains the MS Access databases (`.accdb`) provided by the PNW Bat Hub covering monitoring data from 2016 through 2025.

## Files Overview

| File | Description | Function in Analysis |
| :--- | :--- | :--- |
| `PNW_BatHub_Database_20260414.accdb` | Primary metadata database containing site, deployment, location, clutter, and waterbody records. | Exported to metadata CSVs via `src/00_ExportTables.sh`. |
| `PNW_BatHub_Database_AcousticOutput_2016-2023.accdb` | Database containing acoustic call detections for 2016–2023. | Exported to `calls_2016-2023.csv` via `src/00_ExportTables.sh`. |
| `PNW_BatHub_Database_AcousticOutput_2024.accdb` | Database containing acoustic call detections for 2024. | Exported to `calls_2024.csv` via `src/00_ExportTables.sh`. |
| `PNW_BatHub_Database_AcousticOutput_2025.accdb` | Database containing acoustic call detections for 2025 (WA/OR). | Exported to `calls_2025.csv` via `src/00_ExportTables.sh`. |
| `PNW_BatHub_Database_AcousticOutput_Idaho_2025.accdb` | Database containing acoustic call detections for Idaho in 2025. | Exported to `calls_Idaho_2025.csv` via `src/00_ExportTables.sh`. |

## Subdirectories

- **`archive/`**: Older legacy versions of database files and metadata.

