# Database Tables (data/raw/tables/)

This directory contains CSV exports from the MS Access databases in `data/raw/database/`.

## Files Overview

| File | Description | Function in Analysis |
| :--- | :--- | :--- |
| `calls_from_2024.csv` | Acoustic detections from 2024 to present. | Used in `03_DetectionsCleaning.R` to build detection histories. |
| `calls_to_2024.csv` | Acoustic detections prior to 2024. | Used in `03_DetectionsCleaning.R` to build detection histories. |
| `tblDeployment.csv` | Records of acoustic detector deployments. | Used in `01_DeploymentsCleaning.R` to link sites to dates. |
| `tblPointLocation.csv` | Spatial coordinates for deployment points. | Used in `01_DeploymentsCleaning.R` for mapping and spatial joins. |
| `tblSite.csv` | Information about monitoring sites. | Joined in `01_DeploymentsCleaning.R`. |
| `tluClutterType.csv` | Lookup table for habitat clutter categories. | Used in `01_DeploymentsCleaning.R` for covariate formatting. |
| `tluWaterBodyType.csv` | Lookup table for nearby water body types. | Used in `01_DeploymentsCleaning.R` for covariate formatting. |
