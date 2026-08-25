# Database Tables (`data/raw/tables/`)

This directory contains CSV exports extracted from the MS Access databases in `data/raw/database/` via [`src/00_ExportTables.sh`](file:///src/00_ExportTables.sh).

## Files Overview

| File | Description | Function in Analysis |
| :--- | :--- | :--- |
| `calls_2016-2023.csv` | Acoustic detections from 2016 through 2023. | Cleaned in [`src/03_DetectionsCleaning.R`](file:///src/03_DetectionsCleaning.R) to build detection histories. |
| `calls_2024.csv` | Acoustic detections from 2024. | Cleaned in [`src/03_DetectionsCleaning.R`](file:///src/03_DetectionsCleaning.R). |
| `calls_2025.csv` | Acoustic detections from 2025 (WA/OR). | Cleaned in [`src/03_DetectionsCleaning.R`](file:///src/03_DetectionsCleaning.R). |
| `calls_Idaho_2025.csv` | Acoustic detections from Idaho in 2025. | Cleaned in [`src/03_DetectionsCleaning.R`](file:///src/03_DetectionsCleaning.R). |
| `tblDeployment.csv` | Records of acoustic detector deployments. | Used in [`src/01_DeploymentsCleaning.R`](file:///src/01_DeploymentsCleaning.R) to link deployments to dates and locations. |
| `tblPointLocation.csv` | Spatial coordinates and location names for deployment points. | Used in [`src/01_DeploymentsCleaning.R`](file:///src/01_DeploymentsCleaning.R) for spatial mapping and Daymet site exports. |
| `tblSite.csv` | Metadata regarding monitoring grid cells and sites. | Joined in [`src/01_DeploymentsCleaning.R`](file:///src/01_DeploymentsCleaning.R). |
| `tluClutterType.csv` | Lookup table for habitat clutter categories (0–25%, 26–50%, 51–75%, 76–100%). | Used in [`src/01_DeploymentsCleaning.R`](file:///src/01_DeploymentsCleaning.R) to format detection clutter covariates. |
| `tluWaterBodyType.csv` | Lookup table for nearby water body types. | Used in [`src/01_DeploymentsCleaning.R`](file:///src/01_DeploymentsCleaning.R) to format water indicator covariates. |

## Subdirectories

- **`archive/`**: Older legacy exports.

