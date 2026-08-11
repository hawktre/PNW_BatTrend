## ---------------------------
##
## Purpose of script: Format Detections for Occupancy Modeling
##
## Author: Trent VanHawkins
##
## ---------------------------

## load up the packages we will need:  (uncomment as required)
library(tidyverse)
library(here)
library(janitor)

# Read in deloyment data and join with daymet  --------------------------------------------------
deployment <- readRDS(here("data/processed/detections/deployments_to2025.rds"))

# Read in detection data --------------------------------------------------
## Read in all matching calls files dynamically
calls_files <- list.files(path = here("data/raw/tables"), pattern = "^calls_.*\\.csv$", full.names = TRUE)
all_raw_acoustics <- map_df(calls_files, function(f) {
  df <- data.table::fread(f, select = c("DeploymentID", "Night", "ManualIDSpp1"))
  df$DeploymentID <- as.integer(df$DeploymentID)
  df$Night <- as.character(df$Night)
  df$ManualIDSpp1 <- as.character(df$ManualIDSpp1)
  df
})

## Clean up
acoustics <- all_raw_acoustics %>%
  ## Drop Blanks from Manual SPP ID
  drop_na(ManualIDSpp1) %>%
  ## Fix values and ensure all in same case
  mutate(
    ManualIDSpp1 = case_when(
      ManualIDSpp1 == 'LASCIN' ~ 'LACI',
      ManualIDSpp1 == 'LASNOC' ~ 'LANO',
      ManualIDSpp1 == 'MYOCIL' ~ 'MYCI',
      ManualIDSpp1 == 'MYOEVO' ~ 'MYEV',
      ManualIDSpp1 == 'MYOLUC' ~ 'MYLU',
      ManualIDSpp1 == 'MYOYUM' ~ 'MYYU',
      ManualIDSpp1 == 'MYOCAL' ~ 'MYCA',
      ManualIDSpp1 == 'EPTFUS' ~ 'EPFU',
      ManualIDSpp1 == 'MYOTHY' ~ 'MYTH',
      TRUE ~ ManualIDSpp1
    ),
    ManualIDSpp1 = tolower(ManualIDSpp1),
    Night = mdy_hms(Night),
    Year = year(Night)
  )

## Create a list of possible bat IDs
possible_bats <- c(
  "laci",
  "lano",
  "myev",
  "epfu",
  "myyu",
  "myth",
  "myci",
  "myvo",
  "tabr",
  "anpa",
  "pahe",
  "euma",
  "myca",
  "mylu",
  "coto"
)

# Remove non-bats ---------------------------------------------------------
acoustics <- acoustics %>%
  filter(ManualIDSpp1 %in% possible_bats)


# Join with deployments ---------------------------------------------------

all_detections <- left_join(
  acoustics,
  deployment,
  by = c("DeploymentID" = "id")
)

## Select just the columns we want
detections <- all_detections %>%
  select(names(deployment)[-1], ManualIDSpp1, Night, DeploymentID) %>%
  select(sample_unit_id, location_name, Night, everything()) %>%
  select(-c(deployment_date, recovery_date)) %>%
  clean_names()

##Take the first night in the case of multiple nights
detections <- detections %>%
  group_by(location_name, year) %>%
  slice_min(night) %>%
  ungroup() %>%
  drop_na()

## Drop NA's

## Check that we don't have any more temporal replicates in a year
detections %>%
  select(location_name, year, night) %>%
  distinct() %>%
  group_by(location_name, year) %>%
  summarise(N = n(), .groups = "drop") %>%
  filter(N > 1)

## Save out sites for daymet
daymet_sites <- detections %>%
  select(location_name, latitude, longitude, night) %>%
  distinct() %>%
  mutate(night = date(night))

write.csv(
  daymet_sites,
  here("data/raw/covariates/daymet/daymet_sites.csv"),
  row.names = F
)
# Remove WA TABR ----------------------------------------------------------

##find the record
bad_tabr <- detections %>%
  filter(manual_id_spp1 == "tabr") %>%
  slice_max(latitude, n = 1) %>%
  pull(deployment_id)
##remove bad record
detections <- detections %>% filter(deployment_id != bad_tabr)

# Pivot Wider to get Spp Richness -------------------------------------------------------------
## Pivot Wider
detections_wide <- detections %>%
  drop_na() %>%
  distinct() %>%
  # Create a presence indicator column
  mutate(present = 1) %>%
  pivot_wider(
    names_from = manual_id_spp1,
    values_from = present,
    values_fill = 0
  ) %>%
  select(-deployment_id)


#Write out
saveRDS(
  detections_wide,
  here("data/processed/detections/detection_histories.rds")
)
