#!/usr/bin/env Rscript

# Single entry point: load fitted brms model + original data, generate an interactive
# "model predictions" figure (call-type proportions over age), and save it as a
# static image (PNG) that can be embedded in `meerkats-against-patriarchy.html`.

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

args <- commandArgs(trailingOnly = TRUE)
open_after <- any(args %in% c("--open", "-o"))

paths <- list(
  data_rds = file.path("data", "PROP_data.rds"),
  model_rds = file.path("models", "B_prop.rds"),
  out_png = file.path("figures", "CallProportions.png")
)

if (!file.exists(paths$data_rds)) stop("Missing: ", paths$data_rds)
if (!file.exists(paths$model_rds)) stop("Missing: ", paths$model_rds)

dat <- readRDS(paths$data_rds)
fit <- readRDS(paths$model_rds)

# Match your original plotting logic:
# - REC_AGE is the human-readable age (days)
# - AGE_z is the scaled/centered version used in the model
# - Predictions are made over TREATMENT × SEX × AGE_z, with other covariates held at their means
# - Multi-category output (.category) corresponds to REP / DIG / CC call types

required_cols <- c(
  "REC_AGE",
  "TREATMENT",
  "SEX",
  "AGE_z",
  "WEIGHT_z",
  "COMP_NORM_z",
  "GS_z",
  "RAIN_z"
)
missing <- setdiff(required_cols, names(dat))
if (length(missing)) stop("PROP_data is missing required columns: ", paste(missing, collapse = ", "))

sd_age <- stats::sd(dat$REC_AGE, na.rm = TRUE)
mean_age <- mean(dat$REC_AGE, na.rm = TRUE)

rec_age_c <- seq(30, 130, by = 1)
age_z_vals <- (rec_age_c - mean_age) / sd_age
age_z_vals <- seq(min(age_z_vals), max(age_z_vals), length.out = 20)

newdata <- tidyr::expand_grid(
  TREATMENT = levels(as.factor(dat$TREATMENT)),
  SEX = levels(as.factor(dat$SEX)),
  AGE_z = age_z_vals,
  WEIGHT_z = mean(dat$WEIGHT_z, na.rm = TRUE),
  COMP_NORM_z = mean(dat$COMP_NORM_z, na.rm = TRUE),
  GS_z = mean(dat$GS_z, na.rm = TRUE),
  RAIN_z = mean(dat$RAIN_z, na.rm = TRUE),
  Total_calls = 1
)

prop_pred <- fit %>%
  tidybayes::epred_draws(newdata = newdata, re_formula = NA, robust = TRUE)

prop_pred <- prop_pred %>%
  mutate(
    REC_AGE = AGE_z * sd_age + mean_age,
    Call_prop = .epred,
    Call_type = as.factor(.category),
    SEX = as.factor(SEX),
    TREATMENT = factor(TREATMENT, levels = c("DC", "SC", "DT"))
  )

# If the model has the categories in a different order/name, keep your labeling stable.
call_linetypes <- c("solid", "dotted", "twodash")
names(call_linetypes) <- levels(prop_pred$Call_type)[seq_len(min(3, nlevels(prop_pred$Call_type)))]

p <- ggplot(prop_pred, aes(x = REC_AGE, y = Call_prop, color = TREATMENT, fill = TREATMENT, linetype = Call_type)) +
  ggdist::stat_lineribbon(.width = 0.95, alpha = 0.3) +
  scale_linetype_manual(values = call_linetypes, name = "Call type", labels = c("REP", "DIG", "CC")) +
  labs(x = "Age (days)", y = "Call proportion") +
  scale_x_continuous(breaks = seq(30, 130, by = 10)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  facet_wrap(~SEX) +
  theme_minimal(base_size = 12)

ggplot2::ggsave(
  filename = paths$out_png,
  plot = p,
  width = 12,
  height = 7,
  units = "in",
  dpi = 200
)
message("Wrote: ", paths$out_png)

if (open_after) {
  utils::browseURL(normalizePath(paths$out_png, winslash = "/"))
}

