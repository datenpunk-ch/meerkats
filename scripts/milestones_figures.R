#!/usr/bin/env Rscript

# Generate milestone-age figures from three brms models:
# - SEMI: age when DIG > REP
# - DIG_PEAK: age at peak DIG proportion
# - FULL: age when CC > DIG
#
# Outputs (PNG) in figures/:
# - Milestone_SEMI.png
# - Milestone_DIG_PEAK.png
# - Milestone_FULL.png
# - Milestones_Combined.png

options(brms.backend = "cmdstanr")

suppressPackageStartupMessages({
  library(cmdstanr)
  library(brms)
  library(dplyr)
  library(tidyr)
  library(tidybayes)
  library(ggplot2)
  library(ggdist)
})

paths <- list(
  figures_dir = "figures",
  semi_model = file.path("models", "SEMI_anova.rds"),
  dig_peak_model = file.path("models", "DIG_PEAK_anova.rds"),
  full_model = file.path("models", "FULL_anova.rds"),
  semi_data = file.path("data", "SEMI_data.rds"),
  dig_peak_data = file.path("data", "PEAK_data.rds"),
  full_data = file.path("data", "FULL_data.rds")
)

dir.create(paths$figures_dir, showWarnings = FALSE, recursive = TRUE)

required_cols <- c("TREATMENT", "SEX")

read_fit <- function(p) {
  if (!file.exists(p)) stop("Missing model file: ", p, "\n(Models are ignored by git; ensure it exists locally.)")
  readRDS(p)
}

read_dat <- function(p) {
  if (!file.exists(p)) stop("Missing data file: ", p)
  d <- readRDS(p)
  miss <- setdiff(required_cols, names(d))
  if (length(miss)) stop("Data file ", p, " missing columns: ", paste(miss, collapse = ", "))
  d
}

summarize_milestone <- function(fit, dat, milestone_id, milestone_label) {
  newdata <- tidyr::expand_grid(
    TREATMENT = levels(as.factor(dat$TREATMENT)),
    SEX = levels(as.factor(dat$SEX))
  )

  # Use epred as the posterior expected milestone age, marginalizing random effects.
  draws <- fit %>%
    tidybayes::epred_draws(newdata = newdata, re_formula = NA, robust = TRUE) %>%
    rename(Milestone_Age = .epred) %>%
    mutate(
      Milestone = milestone_id,
      Milestone_label = milestone_label,
      SEX = as.factor(SEX),
      TREATMENT = factor(TREATMENT, levels = c("DC", "SC", "DT"))
    )

  draws
}

plot_single <- function(draws, title) {
  ggplot(draws, aes(x = TREATMENT, y = Milestone_Age, color = TREATMENT)) +
    ggdist::stat_pointinterval(.width = 0.95, point_size = 2, interval_size = 0.8) +
    facet_wrap(~SEX) +
    labs(x = "Maternal treatment", y = "Milestone age (days)", title = title) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
}

plot_combined <- function(draws_all) {
  ggplot(
    draws_all,
    aes(x = Milestone_label, y = Milestone_Age, color = TREATMENT)
  ) +
    ggdist::stat_pointinterval(
      .width = 0.95,
      point_size = 2,
      interval_size = 0.8,
      position = position_dodge(width = 0.55)
    ) +
    facet_wrap(~SEX) +
    labs(x = "Developmental milestone", y = "Milestone age (days)", color = "Maternal treatment") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
}

semi_fit <- read_fit(paths$semi_model)
dig_peak_fit <- read_fit(paths$dig_peak_model)
full_fit <- read_fit(paths$full_model)

semi_dat <- read_dat(paths$semi_data)
dig_peak_dat <- read_dat(paths$dig_peak_data)
full_dat <- read_dat(paths$full_data)

semi_draws <- summarize_milestone(semi_fit, semi_dat, "SEMI", "SEMI (DIG > REP)")
dig_peak_draws <- summarize_milestone(dig_peak_fit, dig_peak_dat, "DIG_PEAK", "PEAK DIG")
full_draws <- summarize_milestone(full_fit, full_dat, "FULL", "FULL (CC > DIG)")

draws_all <- bind_rows(semi_draws, dig_peak_draws, full_draws)

ggsave(
  filename = file.path(paths$figures_dir, "Milestone_SEMI.png"),
  plot = plot_single(semi_draws, "Milestone age — SEMI (DIG > REP)"),
  width = 10,
  height = 6,
  units = "in",
  dpi = 200
)

ggsave(
  filename = file.path(paths$figures_dir, "Milestone_DIG_PEAK.png"),
  plot = plot_single(dig_peak_draws, "Milestone age — PEAK DIG"),
  width = 10,
  height = 6,
  units = "in",
  dpi = 200
)

ggsave(
  filename = file.path(paths$figures_dir, "Milestone_FULL.png"),
  plot = plot_single(full_draws, "Milestone age — FULL (CC > DIG)"),
  width = 10,
  height = 6,
  units = "in",
  dpi = 200
)

ggsave(
  filename = file.path(paths$figures_dir, "Milestones_Combined.png"),
  plot = plot_combined(draws_all),
  width = 12,
  height = 6,
  units = "in",
  dpi = 200
)

message("Wrote milestone figures to: ", normalizePath(paths$figures_dir, winslash = "/"))

