## ---------------------------
## Purpose of script: Compare JAGS, tPGOcc, and stPGOcc results for all
##                    species to assess sensitivity to model parameterization
##                    and spatial random effects. Includes trend classification
##                    comparison across all three models.
##
## Author: Trent VanHawkins
## ---------------------------

# Load Packages -----------------------------------------------------------
library(tidyverse)
library(here)

# Color Palettes ----------------------------------------------------------
analysis_cols <- wesanderson::wes_palette("Darjeeling1")
state_colors <- wesanderson::wes_palette("FantasticFox1")

trend_cols <- c(
  "Decreasing" = '#D7191C',
  "Increasing" = '#1A9641',
  "Uncertain" = '#D3D3D3'
)

# Load Data ---------------------------------------------------------------
jags_res <- readRDS(here("data/processed/results/jags/all_res.rds"))
tpg_res <- readRDS(here("data/processed/results/tPGOcc/all_res_tpg.rds"))
stpg_res <- readRDS(here("data/processed/results/stPGOcc/all_res_stpg.rds"))

# Combine Results ---------------------------------------------------------
## Filter JAGS to full analysis only, then bind all three models

comparison_res <- list(
  psi = bind_rows(
    jags_res$psi %>% filter(analysis == "OR|WA|ID"),
    tpg_res$psi,
    stpg_res$psi
  ),
  psi_bystate = bind_rows(
    jags_res$psi_bystate %>% filter(analysis == "OR|WA|ID"),
    tpg_res$psi_bystate,
    stpg_res$psi_bystate
  ),
  trend = bind_rows(
    jags_res$trend %>% filter(analysis == "OR|WA|ID"),
    tpg_res$trend,
    stpg_res$trend
  ),
  trend_bystate = bind_rows(
    jags_res$trend_bystate %>% filter(analysis == "OR|WA|ID"),
    tpg_res$trend_bystate,
    stpg_res$trend_bystate
  ),
  params = bind_rows(
    jags_res$params %>% filter(analysis == "OR|WA|ID"),
    tpg_res$params,
    stpg_res$params
  )
)

# Dynamic Year Range ------------------------------------------------------
year_range <- comparison_res$psi %>%
  summarise(min_year = min(year), max_year = max(year))

year_label <- paste0(year_range$min_year, "\u2013", year_range$max_year)

# Model label key ---------------------------------------------------------
model_labels <- c(
  "OR|WA|ID" = "JAGS",
  "tPGOcc" = "tPGOcc",
  "stPGOcc" = "stPGOcc"
)

# Trend Key ---------------------------------------------------------------
## Derived from JAGS full model results for species ordering consistency
trend_key <- jags_res$trend %>%
  filter(analysis == "OR|WA|ID") %>%
  mutate(
    trend = case_when(
      lci > 1 ~ "Increasing",
      uci < 1 ~ "Decreasing",
      TRUE ~ "Uncertain"
    )
  ) %>%
  select(species, trend) %>%
  distinct()

# Output Directory --------------------------------------------------------
dir.create(
  here("Background/presentation_figs/comparison/"),
  recursive = TRUE,
  showWarnings = FALSE
)

# Trend Classification Comparison -----------------------------------------

classify_trend <- function(trend_df, analysis_label) {
  trend_df %>%
    mutate(
      trend = case_when(
        uci < 1 ~ "Decreasing",
        lci > 1 ~ "Increasing",
        TRUE ~ "Uncertain"
      ),
      model = analysis_label
    ) %>%
    select(species, model, mean, lci, uci, trend)
}

trend_comparison <- bind_rows(
  classify_trend(jags_res$trend %>% filter(analysis == "OR|WA|ID"), "JAGS"),
  classify_trend(tpg_res$trend, "tPGOcc"),
  classify_trend(stpg_res$trend, "stPGOcc")
) %>%
  select(species, model, trend) %>%
  pivot_wider(names_from = model, values_from = trend) %>%
  mutate(
    all_agree = JAGS == tPGOcc & JAGS == stPGOcc,
    jags_tpg = JAGS == tPGOcc,
    jags_stpg = JAGS == stPGOcc,
    tpg_stpg = tPGOcc == stPGOcc
  ) %>%
  arrange(JAGS, species)

cat("--- Trend Classification: All Three Models ---\n")
print(
  trend_comparison %>% select(species, JAGS, tPGOcc, stPGOcc, all_agree),
  n = Inf
)

cat("\nSpecies where all three models agree:\n")
trend_comparison %>%
  filter(all_agree) %>%
  select(species, JAGS) %>%
  print(n = Inf)

cat("\nSpecies where models disagree:\n")
trend_comparison %>%
  filter(!all_agree) %>%
  select(species, JAGS, tPGOcc, stPGOcc) %>%
  print(n = Inf)

saveRDS(trend_comparison, here("data/processed/results/trend_comparison.rds"))

## Trend classification heatmap
trend_comparison %>%
  select(species, JAGS, tPGOcc, stPGOcc) %>%
  pivot_longer(
    cols = c(JAGS, tPGOcc, stPGOcc),
    names_to = "model",
    values_to = "trend"
  ) %>%
  mutate(model = factor(model, levels = c("JAGS", "tPGOcc", "stPGOcc"))) %>%
  ggplot(aes(x = model, y = species, fill = trend)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(values = trend_cols) +
  labs(
    title = paste("Trend Classification Comparison", year_label),
    x = "Model",
    y = "Species",
    fill = "Trend"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "bottom"
  )

ggsave(
  "comparison_trend_heatmap.png",
  path = here("Background/presentation_figs/comparison/"),
  width = 5,
  height = 7,
  dpi = 300
)

# Parameter Figures -------------------------------------------------------

## Occupancy coefficients
shared_occ_params <- c("log_fc", "precip", "dem_max", "log_cliff")

comparison_res$params %>%
  filter(param %in% shared_occ_params) %>%
  mutate(
    param = case_when(
      param == "log_fc" ~ "Log(Forest Cover)",
      param == "precip" ~ "Precipitation",
      param == "dem_max" ~ "Elevation",
      param == "log_cliff" ~ "Log(Cliff Cover)"
    )
  ) %>%
  ggplot(aes(x = species, y = mean, group = analysis)) +
  geom_errorbar(
    aes(ymin = lci, ymax = uci),
    width = 0,
    linewidth = 0.75,
    position = position_dodge(width = 0.6)
  ) +
  geom_errorbar(
    aes(ymin = q25, ymax = q75, color = analysis),
    width = 0,
    linewidth = 1.5,
    position = position_dodge(width = 0.6)
  ) +
  geom_point(position = position_dodge(width = 0.6)) +
  facet_grid(param ~ ., scales = "free") +
  geom_hline(yintercept = 0, lty = 2) +
  scale_color_manual(values = analysis_cols[1:3], labels = model_labels) +
  labs(
    title = "Occupancy Coefficients: All Models",
    subtitle = year_label,
    x = "Species",
    y = "Posterior Distribution (Log-Odds)",
    color = "Model"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "bottom"
  )

ggsave(
  "comparison_alphas.png",
  path = here("Background/presentation_figs/comparison/"),
  width = 8,
  height = 7,
  dpi = 300
)

## Detection coefficients
shared_det_params <- c(
  "(Intercept)",
  "tmin",
  "dayl",
  "water_ind",
  "scale(tmin)",
  "scale(dayl)",
  "clutter_percent",
  "clutter_percent1",
  "clutter_percent2",
  "clutter_percent3",
  "clutter_percent4"
)

comparison_res$params %>%
  filter(param %in% shared_det_params) %>%
  mutate(
    param = case_when(
      param %in% c("tmin", "scale(tmin)") ~ "Min. Temp",
      param %in% c("dayl", "scale(dayl)") ~ "Day Length",
      param == "water_ind" ~ "Waterbody",
      param %in%
        c(
          "clutter_percent",
          "clutter_percent1",
          "clutter_percent2",
          "clutter_percent3",
          "clutter_percent4"
        ) ~ "Clutter",
      TRUE ~ "(Intercept)"
    )
  ) %>%
  ggplot(aes(x = species, y = mean, group = analysis)) +
  geom_errorbar(
    aes(ymin = lci, ymax = uci),
    width = 0,
    linewidth = 0.75,
    position = position_dodge(width = 0.6)
  ) +
  geom_errorbar(
    aes(ymin = q25, ymax = q75, color = analysis),
    width = 0,
    linewidth = 1.5,
    position = position_dodge(width = 0.6)
  ) +
  geom_point(position = position_dodge(width = 0.6)) +
  facet_grid(param ~ ., scales = "free") +
  geom_hline(yintercept = 0, lty = 2) +
  scale_color_manual(values = analysis_cols[1:3], labels = model_labels) +
  labs(
    title = "Detection Coefficients: All Models",
    subtitle = year_label,
    x = "Species",
    y = "Posterior Distribution (Log-Odds)",
    color = "Model"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "bottom"
  )

ggsave(
  "comparison_betas.png",
  path = here("Background/presentation_figs/comparison/"),
  width = 8,
  height = 8,
  dpi = 300
)

# Lambda Figures ----------------------------------------------------------

## Overall lambda
comparison_res$trend %>%
  ggplot(aes(x = species, y = mean, group = analysis)) +
  geom_errorbar(
    aes(ymin = lci, ymax = uci),
    width = 0,
    linewidth = 0.75,
    position = position_dodge(width = 0.6)
  ) +
  geom_errorbar(
    aes(ymin = q25, ymax = q75, color = analysis),
    width = 0,
    linewidth = 1.5,
    position = position_dodge(width = 0.6)
  ) +
  geom_point(position = position_dodge(width = 0.6)) +
  geom_hline(yintercept = 1, lty = 2) +
  scale_color_manual(values = analysis_cols[1:3], labels = model_labels) +
  labs(
    title = paste("Trend Comparison: All Models", year_label),
    x = "Species",
    y = expression("Posterior Distribution (" * lambda[tot] * ")"),
    color = "Model"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "bottom"
  )

ggsave(
  "comparison_lambda.png",
  path = here("Background/presentation_figs/comparison/"),
  width = 8,
  height = 4,
  dpi = 300
)

## Lambda by state
comparison_res$trend_bystate %>%
  ggplot(aes(x = species, y = mean, group = interaction(analysis, state))) +
  geom_errorbar(
    aes(ymin = lci, ymax = uci),
    width = 0,
    linewidth = 0.75,
    position = position_dodge(width = 0.6)
  ) +
  geom_errorbar(
    aes(ymin = q25, ymax = q75, color = state),
    width = 0,
    linewidth = 1.5,
    position = position_dodge(width = 0.6)
  ) +
  geom_point(aes(shape = analysis), position = position_dodge(width = 0.6)) +
  geom_hline(yintercept = 1, lty = 2) +
  scale_color_manual(values = state_colors[3:5]) +
  scale_shape_manual(
    values = c("OR|WA|ID" = 16, "tPGOcc" = 17, "stPGOcc" = 15),
    labels = model_labels
  ) +
  labs(
    title = paste("Trend by State: All Models", year_label),
    x = "Species",
    y = expression("Posterior Distribution (" * lambda[tot] * ")"),
    color = "State",
    shape = "Model"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "bottom"
  )

ggsave(
  "comparison_lambda_bystate.png",
  path = here("Background/presentation_figs/comparison/"),
  width = 10,
  height = 5,
  dpi = 300
)

# Psi Trend Figures -------------------------------------------------------

## Overall psi by species
comparison_res$psi %>%
  left_join(trend_key, by = "species") %>%
  ggplot(aes(
    x = year,
    y = mean,
    color = analysis,
    fill = analysis,
    group = analysis
  )) +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  lims(y = c(0, 1)) +
  facet_wrap(~species, ncol = 3) +
  scale_color_manual(values = analysis_cols[1:3], labels = model_labels) +
  scale_fill_manual(values = analysis_cols[1:3], labels = model_labels) +
  labs(
    title = paste("Occupancy Trends: All Models", year_label),
    x = "Year",
    y = expression("Posterior Distribution (" * bar(psi)[t] * ")"),
    color = "Model",
    fill = "Model"
  ) +
  theme_bw() +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(
  "comparison_psi_overall.png",
  path = here("Background/presentation_figs/comparison/"),
  width = 10,
  height = 10,
  dpi = 300
)

## Psi by state
comparison_res$psi_bystate %>%
  left_join(trend_key, by = "species") %>%
  ggplot(aes(
    x = year,
    y = mean,
    color = analysis,
    fill = analysis,
    group = analysis
  )) +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  lims(y = c(0, 1)) +
  facet_grid(species ~ state) +
  scale_color_manual(values = analysis_cols[1:3], labels = model_labels) +
  scale_fill_manual(values = analysis_cols[1:3], labels = model_labels) +
  labs(
    title = paste("Occupancy Trends by State: All Models", year_label),
    x = "Year",
    y = expression("Posterior Distribution (" * bar(psi)[t] * ")"),
    color = "Model",
    fill = "Model"
  ) +
  theme_bw() +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(
  "comparison_psi_bystate.png",
  path = here("Background/presentation_figs/comparison/"),
  width = 10,
  height = 14,
  dpi = 300
)
