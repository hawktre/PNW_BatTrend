
# Overview:
# Build csv tables from Access database files
# Databases are originally stored in Box folder (HERS_Working/Bats/DatabaseBE)
# Databases used in the analysis are moved to data/raw/database
# This script pulls data out of database files and saves to data/raw/tables

# load packages
library(tidyverse)
library(here)
library(RODBC)

# ----- set input/output settings ----------------------------------------------
## set file path
filepath <- here("data/raw/database")

## set db file locations
### adjust file locations if different than below
db_names <- tribble(
  ~"db",         ~"filename",
  "SRC",        "PNW_BatHub_Database_20260414.accdb",
  "SRC_CALLS1", "PNW_BatHub_Database_AcousticOutput_2016-2023.accdb",
  "SRC_CALLS2", "PNW_BatHub_Database_AcousticOutput_2024.accdb",
  "SRC_CALLS3", "PNW_BatHub_Database_AcousticOutput_2025.accdb",
  "SRC_CALLS4", "PNW_BatHub_Database_AcousticOutput_Idaho_2025.accdb",
  )

## build filepaths
db_names$db_locs <- paste0(rep(filepath, nrow(db_names)),
                           "/", 
                           db_names$filename)

## set output location for raw data tables
OUT <- "data/raw/tables"

## check for input directory
if(!dir.exists(filepath)){
  stop("Input directory '", filepath, "' not found.")
}

## ensure output folder exists
if(!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)

cat("Extracting Bat Hub data from", filepath, "to", OUT, "...\n")

# ----- load databases and convert to csv --------------------------------------
## site metadata db -----------------------------------
db <- db_names %>% filter(db == "SRC") %>% pull(db_locs)

if(!file.exists(db)){
  warning("Database file not found: ", db) 
} else{
  cat("Exporting site metadata tables from '", 
      basename(db), 
      "'...\n", 
      sep = "")
  
  # load connection
  src <- RODBC::odbcConnectAccess2007(db)
  
  ## read in Access Tables
  tblDeployment <- RODBC::sqlFetch(src,"tblDeployment")
  tblPointLocation <- RODBC::sqlFetch(src, "tblPointLocation")
  tblSite <- RODBC::sqlFetch(src, "tblSite")
  tluClutterType <- RODBC::sqlFetch(src, "tluClutterType")
  tluWaterBodyType <-RODBC::sqlFetch(src, "tluWaterBodyType")
  
  ## close db
  RODBC::odbcClose(src)
  
  # write to csv
  ## write out to csv
  write_csv(tblDeployment,    file.path(OUT, "tblDeployment.csv"))
  write_csv(tblPointLocation, file.path(OUT, "tblPointLocation.csv"))
  write_csv(tblSite,          file.path(OUT, "tblSite.csv"))
  write_csv(tluClutterType,   file.path(OUT, "tluClutterType.csv"))
  write_csv(tluWaterBodyType, file.path(OUT, "tluWaterBodyType.csv"))
  
  ## remove object
  rm(src, tblDeployment, tblPointLocation, tblSite, tluClutterType, tluWaterBodyType) 
}

## Acoustic output databases ---------------------------------
calls_db_locs <- db_names %>% 
  filter(str_detect(filename, "AcousticOutput")) %>% 
  pull(db_locs)

if (length(calls_db_locs) == 0) {
  warning("No acoustic output databases found in db_names.")
}

for(db in calls_db_locs) {
  filename <- basename(db)
  
  if(!file.exists(db)) {
    warning("Database file not found: ", db)
    next
  }
  
  # suffix
  suffix <- str_remove(filename, "^.*AcousticOutput_")
  suffix <- str_remove(suffix, "\\.accdb$")
  
  dbTables <- RODBC::odbcConnectAccess2007(db)
  
  # find all tables that start with tblDeploymentDetection
  all_tables <- RODBC::sqlTables(dbTables, tableType = "TABLE")$TABLE_NAME
  tbls <- all_tables[str_starts(all_tables, "tblDeploymentDetection")]
  
  if (length(tbls) == 0) {
    warning("No table starting with 'tblDeploymentDetection' found in '", filename, "'.")
  } else {
    for (tbl in tbls) {
      cat("Exporting acoustic output table '", tbl, "' from '", filename,
          "' to 'calls_", suffix, ".csv'...\n", sep = "")
      
      dat <- RODBC::sqlFetch(dbTables, tbl)
      # export to csv in output 
      write_csv(dat, file.path(OUT, paste0("calls_", suffix, ".csv")))
      rm(dat)
    }
  }

  RODBC::odbcClose(dbTables)
  rm(dbTables, db, filename, suffix, all_tables, tbls)
  }
cat("Export complete.\n")

