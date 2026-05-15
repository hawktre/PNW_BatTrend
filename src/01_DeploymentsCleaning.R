## ---------------------------
## Purpose of script: Read-in and format data for trend analysis
##
## Author: Trent VanHawkins
## ---------------------------

## ---------------------------

## load up the packages we will need:  (uncomment as required)

library(tidyverse)
library(here)
library(janitor)
library(sf)
library(tigris)

# Read in Required Data ---------------------------------------------------
tblDeployment <- read.csv(here("data/raw/tables/tblDeployment.csv"))
tblPointLocation <- read.csv(here("data/raw/tables/tblPointLocation.csv"))
tblSite <- read.csv(here("data/raw/tables/tblSite.csv"))
tluClutter <- read.csv(here("data/raw/tables/tluClutterType.csv"))
tluWaterBodyType <- read.csv(here("data/raw/tables/tluWaterBodyType.csv"))

# Join all tables together ------------------------------------------------
all_join <- left_join(
  tblDeployment,
  tblPointLocation,
  by = join_by(PointLocationID == ID)
) %>%
  left_join(tblSite, by = join_by(SiteID == ID)) %>%
  left_join(tluClutter, by = join_by(ClutterTypeID == ID)) %>%
  left_join(tluWaterBodyType, by = join_by(WaterBodyTypeID == ID))

# Select only the columns we need -----------------------------------------

deployment <- all_join %>%
  select(
    ID,
    SampleUnitID,
    LocationName,
    Latitude,
    Longitude,
    DeploymentDate,
    RecoveryDate,
    Label.x,
    Label.y,
    ClutterPercent
  ) %>%
  rename("ClutterType" = "Label.x", "WaterBodyType" = "Label.y")


# Make Dates ----------------------------------------
deployment$DeploymentDate <- as_date(as_datetime(
  deployment$DeploymentDate,
  format = "%m/%d/%y %H:%M:%S"
))
deployment$RecoveryDate <- as_date(as_datetime(
  deployment$RecoveryDate,
  format = "%m/%d/%y %H:%M:%S"
))
deployment$year <- year(deployment$DeploymentDate)

# Create State column -----------------------------------------------------
all_states <- states(cb = TRUE)
states_map <- all_states %>%
  filter(NAME %in% c("Oregon", "Washington", "Idaho"))

deployment <- deployment %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_join(states_map %>% select(NAME) %>% st_transform(4326)) %>%
  filter(!is.na(NAME)) %>%
  rename(state = NAME) %>%
  mutate(
    Longitude = st_coordinates(.)[, 1],
    Latitude = st_coordinates(.)[, 2]
  ) %>%
  st_drop_geometry()

# Correct  NAs and wrong values -----------------------------
## Clutter Percent
unique(deployment$ClutterPercent)

deployment <- deployment %>%
  mutate(
    ClutterPercent = factor(case_when(
      ClutterPercent ==
        "0% (no structural interference, e.g., open habitat)" ~ "0",
      ClutterPercent == "1 to 25%" ~ "1",
      ClutterPercent == "26 to 50%" | ClutterPercent == "26-50" ~ "2",
      ClutterPercent == "<null>" ~ NA,
      TRUE ~ ClutterPercent
    ))
  )
sum(is.na(deployment$ClutterPercent))

## Clutter Type
unique(deployment$ClutterType)
deployment <- deployment %>%
  mutate(ClutterType = if_else(ClutterType == "<null>", NA, ClutterType))

## Water Bodies
unique(deployment$WaterBodyType)
### Create Water Indicator
deployment <- deployment %>%
  mutate(water_ind = if_else(WaterBodyType == "None", 0, 1))
### Correct Water Indicator NA's
deployment <- deployment %>%
  mutate(
    water_ind = if_else(
      is.na(water_ind) & str_detect(ClutterType, "Water"),
      1,
      water_ind
    )
  )

# How many rows are we dropping?  -----------------------------------------
## Write CSV for missing locations
missing_sites <- deployment %>% filter(is.na(ClutterPercent) | is.na(water_ind))
write.csv(
  missing_sites,
  here("data/processed/detections/sites_missing_covars.csv")
)

## Create plot of missing sites
states_map <- states(cb = TRUE) %>%
  filter(NAME %in% c("Oregon", "Washington", "Idaho"))

missing_sites %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
  ggplot() +
  geom_sf(data = states_map, fill = "gray90", color = "black") +
  geom_sf(aes(color = factor(year)))

ggsave(
  filename = "missing_clutterpercent.png",
  path = here("data/processed/detections/"),
  height = 6,
  width = 6,
  units = "in"
)

## How many sites?
n_drop <- sum(is.na(deployment$ClutterPercent) | is.na(deployment$water_ind))
p_drop <- n_drop / nrow(deployment)

cat("Dropping", n_drop, "rows missing clutter percent or waterbody indicator")

## Drop the missing rows
deployment <- deployment %>%
  drop_na(ClutterPercent, water_ind) %>%
  select(-c(ClutterType, WaterBodyType))

deployment <- clean_names(deployment) %>%
  arrange(sample_unit_id, deployment_date)

deployment %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  ggplot() +
  geom_sf(data = states_map, fill = "gray90", color = "black") +
  geom_sf(aes(color = state)) +
  facet_wrap(~year) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(color = "State")

# Save out for verification in excel --------------------------------------
saveRDS(deployment, here("data/processed/detections/deployments_to2024.rds"))


# Save out sites for daymet -----------------------------------------------
daymet_sites <- deployment %>%
  select(location_name, deployment_date, latitude, longitude) %>%
  distinct() |>
  rename("night" = deployment_date)

## Save out for daymet
write_csv(daymet_sites, here("data/raw/covariates/daymet/daymet_sites.csv"))
