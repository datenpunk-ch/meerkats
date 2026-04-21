#!/usr/bin/env Rscript

# Interactive call-proportions plot (custom HTML + Plotly).
# - Posterior predictions on REC_AGE = 30, 35, …, 130 days (5-day steps)
# - Single plot with checkbox filters (including call type)
# - Treatment = color (Okabe–Ito mapping; DC/SC/DT uses order c(2,1,3))
# - Sex = marker symbol + line dash
# - Checkbox filters for Treatment and Sex, plus "Averaged" (auto-managed)
# - Optional 95% credible-interval ribbons toggle (behind lines)
# - Custom tooltip + continuous crosshair line

options(brms.backend = "cmdstanr")

suppressPackageStartupMessages({
  library(cmdstanr)
  library(brms)
  library(dplyr)
  library(tidyr)
  library(tidybayes)
  library(jsonlite)
})

paths <- list(
  figures_dir = "figures",
  data_rds = file.path("data", "PROP_data.rds"),
  model_rds = file.path("models", "B_prop.rds"),
  out_html = file.path("figures", "CallProportions.html")
)

dir.create(paths$figures_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(paths$data_rds)) stop("Missing: ", paths$data_rds)
if (!file.exists(paths$model_rds)) stop("Missing model file: ", paths$model_rds, "\n(Models are ignored by git; ensure it exists locally.)")

dat <- readRDS(paths$data_rds)
fit <- readRDS(paths$model_rds)

# Avoid brms trying to touch rstan backends.
options(brms.backend = "cmdstanr")

required_cols <- c("REC_AGE", "TREATMENT", "SEX", "AGE_z")
missing <- setdiff(required_cols, names(dat))
if (length(missing)) stop("PROP_data is missing required columns: ", paste(missing, collapse = ", "))

# Optional covariates: if present, hold at mean.
# Note: different extracts may name competition normalization slightly differently.
maybe_covariates <- c(
  "WEIGHT_z",
  "COMP_NORM_z",
  "COMP_z",
  "GS_z",
  "RAIN_z",
  "Total_calls",
  "TOTAL_CALLS"
)
covariates_present <- intersect(maybe_covariates, names(dat))
# Prefer COMP_NORM_z if both exist.
if (all(c("COMP_NORM_z", "COMP_z") %in% covariates_present)) {
  covariates_present <- setdiff(covariates_present, "COMP_z")
}
# Prefer Total_calls naming if both exist.
if (all(c("Total_calls", "TOTAL_CALLS") %in% covariates_present)) {
  covariates_present <- setdiff(covariates_present, "TOTAL_CALLS")
}

sex_levels <- levels(as.factor(dat[["SEX"]]))
treat_levels <- levels(as.factor(dat[["TREATMENT"]]))

# Prefer a stable, meaningful treatment order when applicable.
preferred_treat_order <- c("DC", "SC", "DT")
if (all(preferred_treat_order %in% treat_levels)) {
  treat_levels <- c(preferred_treat_order, setdiff(treat_levels, preferred_treat_order))
}

make_palette <- function(levels) {
  okabe_ito <- c(
    orange = "#E69F00",
    skyblue = "#56B4E9",
    bluishgreen = "#009E73",
    yellow = "#F0E442",
    blue = "#0072B2",
    vermillion = "#D55E00",
    reddishpurple = "#CC79A7",
    black = "#000000"
  )

  # Match: scale_color_okabe_ito(order = c(2, 1, 3), labels = c('DC','SC','DT'))
  if (all(c("DC", "SC", "DT") %in% levels)) {
    out <- stats::setNames(rep(NA_character_, length(levels)), levels)
    out[["DC"]] <- okabe_ito[["skyblue"]]
    out[["SC"]] <- okabe_ito[["orange"]]
    out[["DT"]] <- okabe_ito[["bluishgreen"]]

    rest <- setdiff(levels, c("DC", "SC", "DT"))
    pool <- unname(okabe_ito[c("yellow", "blue", "vermillion", "reddishpurple", "black")])
    if (length(rest)) out[rest] <- rep(pool, length.out = length(rest))
    return(out)
  }

  cols <- unname(okabe_ito)
  stats::setNames(rep(cols, length.out = length(levels)), levels)
}

palette <- make_palette(treat_levels)

# Build age grid in REC_AGE space, then map to AGE_z using sample mean/sd.
sd_age <- stats::sd(dat[["REC_AGE"]], na.rm = TRUE)
mean_age <- mean(dat[["REC_AGE"]], na.rm = TRUE)
if (!is.finite(sd_age) || sd_age <= 0) stop("REC_AGE has invalid sd; cannot map to AGE_z.")

# Fixed prediction window: 30, 35, …, 130 days (5-day steps) for interactive display.
pred_age_lo <- 30
pred_age_hi <- 130
pred_age_by <- 5
rec_age_grid <- seq(pred_age_lo, pred_age_hi, by = pred_age_by)
age_z_grid <- (rec_age_grid - mean_age) / sd_age
age_lookup <- data.frame(AGE_z = age_z_grid, REC_AGE = as.numeric(rec_age_grid))

newdata <- tidyr::expand_grid(
  TREATMENT = treat_levels,
  SEX = sex_levels,
  AGE_z = age_z_grid
)

if (length(covariates_present)) {
  for (v in covariates_present) {
    if (v %in% c("Total_calls", "TOTAL_CALLS")) {
      newdata[[v]] <- 1
    } else {
      newdata[[v]] <- mean(dat[[v]], na.rm = TRUE)
    }
  }
}

draws <- fit %>%
  # Limit draws for responsiveness: interactive plots don't need full posterior size.
  tidybayes::epred_draws(newdata = newdata, re_formula = NA, robust = TRUE, ndraws = 600) %>%
  dplyr::left_join(age_lookup, by = "AGE_z") %>%
  mutate(
    Call_prop = .epred,
    Call_type = as.factor(.category),
    SEX = as.factor(SEX),
    TREATMENT = factor(TREATMENT, levels = treat_levels)
  )

# Ensure proportions sum to 1 at each age within each draw/treatment/sex.
draws <- draws %>%
  group_by(.draw, TREATMENT, SEX, AGE_z, REC_AGE) %>%
  mutate(Call_prop = Call_prop / sum(Call_prop, na.rm = TRUE)) %>%
  ungroup()

call_levels <- levels(droplevels(draws$Call_type))
preferred_call_order <- c("REP", "DIG", "CC")
if (all(preferred_call_order %in% call_levels)) {
  call_levels <- c(preferred_call_order, setdiff(call_levels, preferred_call_order))
}

label_for_call_level <- function(level) {
  lvl <- as.character(level)
  u <- toupper(lvl)
  if (u == "REP" || grepl("^REP", u) || grepl("REPEAT", u) || grepl("SUMBEG", u) || grepl("BEG", u)) {
    return("REP")
  }
  # brms category labels may be SumDIG / SumCC (not prefixed with DIG / CC).
  if (u == "DIG" || grepl("^DIG", u) || grepl("SUMDIG", u)) {
    return("DIG")
  }
  if (u == "CC" || grepl("^CC", u) || grepl("CLOSE", u) || grepl("SUMCC", u)) {
    return("CC")
  }
  lvl
}

call_labels <- stats::setNames(vapply(call_levels, label_for_call_level, character(1)), call_levels)

summ <- draws %>%
  group_by(TREATMENT, SEX, Call_type, REC_AGE) %>%
  tidybayes::median_qi(Call_prop, .width = 0.95) %>%
  ungroup() %>%
  mutate(Call_type = factor(Call_type, levels = call_levels))

# After summarising, re-scale point estimates + ribbons so the three categories
# always sum to 1 at each age (compositional constraint), while staying in [0,1].
# Use the *same* scale for median and quantiles (sum of medians); scaling .lower
# and .upper by sum(.lower) / sum(.upper) separately breaks interval coherence
# and can invert bands or shift ribbon edges (mushy / wrong-looking CIs).
rescale_compositional <- function(df) {
  df %>%
    group_by(TREATMENT, SEX, REC_AGE) %>%
    mutate(
      s_m = sum(Call_prop, na.rm = TRUE),
      Call_prop = dplyr::if_else(s_m > 0, Call_prop / s_m, Call_prop),
      .lower = dplyr::if_else(s_m > 0, .lower / s_m, .lower),
      .upper = dplyr::if_else(s_m > 0, .upper / s_m, .upper)
    ) %>%
    ungroup() %>%
    mutate(
      .lower = pmax(0, pmin(1, .lower)),
      .upper = pmax(0, pmin(1, .upper)),
      .upper = pmax(.upper, .lower),
      Call_prop = pmax(0, pmin(1, Call_prop))
    ) %>%
    select(-s_m)
}

summ <- rescale_compositional(summ)

# Sex-marginalised summary ("Averaged") by averaging within each posterior draw.
summ_allsex <- draws %>%
  group_by(.draw, TREATMENT, Call_type, REC_AGE) %>%
  summarise(Call_prop = mean(Call_prop, na.rm = TRUE), .groups = "drop") %>%
  group_by(.draw, TREATMENT, REC_AGE) %>%
  mutate(Call_prop = Call_prop / sum(Call_prop, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(TREATMENT, Call_type, REC_AGE) %>%
  tidybayes::median_qi(Call_prop, .width = 0.95) %>%
  ungroup() %>%
  mutate(
    SEX = factor("ALL", levels = "ALL"),
    Call_type = factor(Call_type, levels = call_levels)
  )

summ_allsex <- rescale_compositional(summ_allsex)

make_traces <- function(summ, summ_allsex, sex_levels, treat_levels, call_levels, palette) {
  traces <- list()

  sex_symbol <- setNames(rep("circle", length(sex_levels)), sex_levels)
  sex_dash <- setNames(rep("solid", length(sex_levels)), sex_levels)
  if ("F" %in% sex_levels) {
    sex_symbol[["F"]] <- "circle"
    sex_dash[["F"]] <- "solid"
  }
  if ("M" %in% sex_levels) {
    sex_symbol[["M"]] <- "diamond"
    sex_dash[["M"]] <- "dot"
  }
  sex_symbol[["ALL"]] <- "square"
  sex_dash[["ALL"]] <- "solid"

  # Prediction ages are every 5 d (30 … 130); show a marker at each predicted age.
  msz_line <- function(n) rep(6, n)

  for (ct in call_levels) {
    for (sx in sex_levels) {
      for (tr in treat_levels) {
        df <- summ %>% filter(Call_type == ct, SEX == sx, TREATMENT == tr) %>% arrange(REC_AGE)
        if (nrow(df) == 0) next

        col <- palette[[as.character(tr)]]
        if (is.null(col)) col <- "#333333"

        x0 <- df$REC_AGE
        # Plot on a 0–100% scale (still compositional in the underlying [0,1] space).
        y0 <- round(df$Call_prop * 100, 4)
        lo0 <- round(df$.lower * 100, 4)
        hi0 <- round(df$.upper * 100, 4)

        msz <- msz_line(length(x0))

        # Ribbon (linear segments between ages — avoids spline overshoot on sparse x)
        traces[[length(traces) + 1]] <- list(
          x = c(x0, rev(x0)),
          y = c(lo0, rev(hi0)),
          type = "scatter",
          mode = "lines",
          line = list(color = "rgba(0,0,0,0)", width = 0, shape = "linear"),
          fill = "toself",
          fillcolor = paste0(col, "18"),
          hoverinfo = "skip",
          name = paste0(tr, " — ", sx, " 95% CI"),
          meta = list(treatment = as.character(tr), sex = as.character(sx), call_type = as.character(ct), kind = "ribbon")
        )

        # Line+markers
        traces[[length(traces) + 1]] <- list(
          x = x0,
          y = y0,
          type = "scatter",
          mode = "lines+markers",
          cliponaxis = "ticks",
          line = list(color = col, width = 2, dash = sex_dash[[sx]], shape = "linear"),
          marker = list(color = col, size = msz, symbol = sex_symbol[[sx]], line = list(width = 0)),
          name = paste0(tr, " — ", sx),
          hoverinfo = "none",
          hovertemplate = "",
          meta = list(treatment = as.character(tr), sex = as.character(sx), call_type = as.character(ct), kind = "line")
        )
      }
    }
  }

  for (ct in call_levels) {
    for (tr in treat_levels) {
      df <- summ_allsex %>% filter(Call_type == ct, TREATMENT == tr) %>% arrange(REC_AGE)
      if (nrow(df) == 0) next

      col <- palette[[as.character(tr)]]
      if (is.null(col)) col <- "#333333"

      x0 <- df$REC_AGE
      y0 <- round(df$Call_prop * 100, 4)
      lo0 <- round(df$.lower * 100, 4)
      hi0 <- round(df$.upper * 100, 4)

      msz <- msz_line(length(x0))

      traces[[length(traces) + 1]] <- list(
        x = c(x0, rev(x0)),
        y = c(lo0, rev(hi0)),
        type = "scatter",
        mode = "lines",
        line = list(color = "rgba(0,0,0,0)", width = 0, shape = "linear"),
        fill = "toself",
        fillcolor = paste0(col, "18"),
        hoverinfo = "skip",
        name = paste0(tr, " — ALL 95% CI"),
        meta = list(treatment = as.character(tr), sex = "ALL", call_type = as.character(ct), kind = "ribbon")
      )

      traces[[length(traces) + 1]] <- list(
        x = x0,
        y = y0,
        type = "scatter",
        mode = "lines+markers",
        cliponaxis = "ticks",
        line = list(color = col, width = 2, dash = "solid", shape = "linear"),
        marker = list(color = col, size = msz, symbol = "square", line = list(width = 0)),
        name = paste0(tr, " — ALL"),
        hoverinfo = "none",
        hovertemplate = "",
        meta = list(treatment = as.character(tr), sex = "ALL", call_type = as.character(ct), kind = "line")
      )
    }
  }

  traces
}

traces <- make_traces(summ, summ_allsex, sex_levels, treat_levels, call_levels, palette)

payload <- list(
  traces = traces,
  palette = as.list(palette),
  sex_levels = sex_levels,
  treat_levels = treat_levels,
  call_levels = call_levels,
  call_labels = as.list(call_labels),
  x_range = c(pred_age_lo, pred_age_hi)
)

payload_json <- jsonlite::toJSON(payload, auto_unbox = TRUE, digits = 6)

html <- paste0(
  "<!doctype html>\n",
  "<html lang=\"en\">\n",
  "<head>\n",
  "  <meta charset=\"utf-8\" />\n",
  "  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\" />\n",
  "  <title>Call proportions — interactive</title>\n",
  "  <script src=\"https://cdn.plot.ly/plotly-2.35.2.min.js\"></script>\n",
  "  <style>\n",
  "    html, body { height: 100%; }\n",
  "    body { margin: 0; font-family: system-ui, -apple-system, Segoe UI, Roboto, Arial, sans-serif; }\n",
  "    .wrap { padding: 12px 12px 8px; box-sizing: border-box; }\n",
  "    .controls { display: flex; flex-wrap: wrap; gap: 10px 18px; align-items: flex-start; margin-bottom: 10px; }\n",
  "    fieldset { border: 1px solid #ddd; padding: 8px 10px; }\n",
  "    legend { padding: 0 6px; font-weight: 600; }\n",
  "    label { display: inline-flex; align-items: center; gap: 6px; margin-right: 10px; white-space: nowrap; }\n",
  "    input:disabled { cursor: not-allowed; }\n",
  "    label:has(input:disabled) { opacity: 0.55; }\n",
  "    label:has(input:disabled:checked) { opacity: 1; }\n",
  "    label:has(input:disabled:checked) input { accent-color: #000; }\n",
  "    label:has(input:disabled:checked) .sexsample { outline: 2px solid #000; outline-offset: 2px; border-radius: 3px; }\n",
  "    label:has(input:disabled:checked) .sexline { border-top-color: #000; }\n",
  "    label:has(input:disabled:checked) .sexmark { background: #000; }\n",
  "    .plotwrap { position: relative; }\n",
  "    #plot { width: 100%; height: 640px; }\n",
  "    #tooltip { position: absolute; z-index: 20; padding: 8px 10px; border: 1px solid #ddd; background: rgba(255,255,255,0.95); font-size: 12px; color: #111; display: none; max-width: none; width: max-content; pointer-events: none; }\n",
  "    #tooltip table.tttbl { border-collapse: collapse; line-height: 1.35; }\n",
  "    #tooltip .tttbl .ttage { text-align: left; font-weight: 700; padding: 0 0 6px; }\n",
  "    #tooltip .tttbl td.ttpct-max { opacity: 1; }\n",
  "    #tooltip .tttbl td.ttpct-max strong { font-weight: 750; letter-spacing: 0.01em; }\n",
  "    #tooltip .tttbl th.ttpcth { text-align: right; font-weight: 650; padding: 2px 0 5px 1.15em; font-variant-numeric: tabular-nums; color: #333; }\n",
  "    #tooltip .tttbl th.tt-ident-h, #tooltip .tttbl td.tt-ident { padding: 3px 10px 3px 0; text-align: left; vertical-align: middle; white-space: nowrap; }\n",
  "    #tooltip .tttbl td.ttpct { text-align: right; padding: 3px 0 3px 1.15em; font-variant-numeric: tabular-nums; min-width: 3.25em; opacity: 0.95; }\n",
  "    #tooltip .tt-ident-inner { display: inline-flex; align-items: baseline; gap: 6px; }\n",
  "    #tooltip .sym { width: 1.1em; display: inline-block; text-align: center; }\n",
  "    .key { display: grid; grid-template-columns: 90px 1fr; column-gap: 12px; row-gap: 6px; align-items: center; }\n",
  "    .keyitems { display: flex; flex-wrap: wrap; gap: 10px 16px; align-items: center; }\n",
  "    .swatch { display: inline-flex; align-items: center; gap: 8px; }\n",
  "    .chip { width: 14px; height: 14px; background: #000; display: inline-block; }\n",
  "    .sexsample { position: relative; width: 34px; height: 14px; display: inline-block; }\n",
  "    .sexline { position: absolute; left: 0; right: 0; top: 50%; transform: translateY(-50%); border-top: 3px solid #111; }\n",
  "    .sexline--dash { border-top-style: dashed; }\n",
  "    .sexline--dot { border-top-style: dotted; }\n",
  "    .sexmark { position: absolute; left: 50%; top: 50%; width: 10px; height: 10px; background: #111; transform: translate(-50%, -50%); }\n",
  "    .sexmark--diamond { transform: translate(-50%, -50%) rotate(45deg); }\n",
  "  </style>\n",
  "</head>\n",
  "<body>\n",
  "  <div class=\"wrap\">\n",
  "    <div class=\"controls\">\n",
  "      <fieldset>\n",
  "        <legend>Call type</legend>\n",
  "        <div class=\"key\">\n",
  "          <strong style=\"font-size: 12px; letter-spacing: 0.02em\">Show</strong>\n",
  "          <span class=\"keyitems\" id=\"calltype-boxes\"></span>\n",
  "        </div>\n",
  "      </fieldset>\n",
  "      <fieldset>\n",
  "        <legend>Uncertainty</legend>\n",
  "        <label><input type=\"checkbox\" id=\"show-ribbons\" /> Show 95% CI ribbons</label>\n",
  "      </fieldset>\n",
  "      <fieldset>\n",
  "        <legend>Legend</legend>\n",
  "        <div class=\"key\" style=\"margin-bottom: 8px\">\n",
  "          <strong style=\"font-size: 12px; letter-spacing: 0.02em\">Treatment</strong>\n",
  "          <span class=\"keyitems\" id=\"treatment-boxes\"></span>\n",
  "        </div>\n",
  "        <div class=\"key\">\n",
  "          <strong style=\"font-size: 12px; letter-spacing: 0.02em\">Sex</strong>\n",
  "          <span class=\"keyitems\" id=\"sex-boxes\"></span>\n",
  "        </div>\n",
  "      </fieldset>\n",
  "    </div>\n",
  "    <div class=\"plotwrap\">\n",
  "      <div id=\"plot\"></div>\n",
  "      <div id=\"tooltip\"></div>\n",
  "    </div>\n",
  "  </div>\n",
  "  <script>\n",
  "    const payload = ", payload_json, ";\n",
  "    const xRange = (Array.isArray(payload.x_range) && payload.x_range.length === 2)\n",
  "      ? payload.x_range.map(Number)\n",
  "      : null;\n",
  "    const gd = document.getElementById('plot');\n",
  "    const tooltip = document.getElementById('tooltip');\n",
  "    const treatments = payload.treat_levels;\n",
  "    const sexes = payload.sex_levels;\n",
  "    const callTypes = payload.call_levels;\n",
  "    const callLabels = payload.call_labels || {};\n",
  "    const palette = payload.palette || {};\n",
  "    const baseTraces = payload.traces;\n",
  "\n",
  "    function safeId(v) { return String(v).replace(/[^a-zA-Z0-9_-]/g, '_'); }\n",
  "    function addCheckboxRow(container, id, value, labelHtml, checked, disabled=false) {\n",
  "      const lab = document.createElement('label');\n",
  "      lab.className = 'swatch';\n",
  "      const input = document.createElement('input');\n",
  "      input.type = 'checkbox';\n",
  "      input.id = id;\n",
  "      input.value = value;\n",
  "      input.checked = checked;\n",
  "      input.disabled = !!disabled;\n",
  "      lab.appendChild(input);\n",
  "      const span = document.createElement('span');\n",
  "      span.innerHTML = labelHtml;\n",
  "      lab.appendChild(span);\n",
  "      container.appendChild(lab);\n",
  "    }\n",
  "\n",
  "    // Build legend checkboxes\n",
  "    const trBox = document.getElementById('treatment-boxes');\n",
  "    treatments.forEach(tr => {\n",
  "      const col = palette[tr] || '#111';\n",
  "      addCheckboxRow(trBox, `tr-${safeId(tr)}`, tr, `<span class=\\\"chip\\\" style=\\\"background:${col}\\\"></span> ${tr}`, true);\n",
  "    });\n",
  "\n",
  "    const sxBox = document.getElementById('sex-boxes');\n",
  "    function sexSample(label) {\n",
  "      if (label === 'M') return `<span class=\\\"sexsample\\\"><span class=\\\"sexline sexline--dot\\\"></span><span class=\\\"sexmark sexmark--diamond\\\"></span></span> M`;\n",
  "      if (label === 'F') return `<span class=\\\"sexsample\\\"><span class=\\\"sexline\\\"></span><span class=\\\"sexmark\\\" style=\\\"border-radius:50%\\\"></span></span> F`;\n",
  "      return `<span class=\\\"sexsample\\\"><span class=\\\"sexline\\\"></span><span class=\\\"sexmark\\\"></span></span> Averaged`;\n",
  "    }\n",
  "    sexes.forEach(sx => addCheckboxRow(sxBox, `sx-${safeId(sx)}`, sx, sexSample(sx), true, false));\n",
  "    addCheckboxRow(sxBox, 'sx-ALL', 'ALL', sexSample('ALL'), false, true);\n",
  "\n",
  "    // Call type checkboxes\n",
  "    const ctBox = document.getElementById('calltype-boxes');\n",
  "    callTypes.forEach(ct => {\n",
  "      const lbl = callLabels[ct] || ct;\n",
  "      addCheckboxRow(ctBox, `ct-${safeId(ct)}`, ct, lbl, true, false);\n",
  "    });\n",
  "\n",
  "    function getSelected(prefix, values) {\n",
  "      const selected = new Set();\n",
  "      values.forEach(v => {\n",
  "        const cb = document.getElementById(`${prefix}-${safeId(v)}`);\n",
  "        if (cb && cb.checked) selected.add(v);\n",
  "      });\n",
  "      return selected;\n",
  "    }\n",
  "\n",
  "    function syncSexSelection() {\n",
  "      const cbA = document.getElementById('sx-ALL');\n",
  "      if (!cbA) return;\n",
  "      cbA.disabled = true;\n",
  "      const anySexChecked = sexes.some(sx => {\n",
  "        const cb = document.getElementById(`sx-${safeId(sx)}`);\n",
  "        return cb && cb.checked;\n",
  "      });\n",
  "      if (anySexChecked) { cbA.checked = false; return; }\n",
  "      cbA.checked = true;\n",
  "    }\n",
  "\n",
  "    function getActiveSexes() {\n",
  "      const cbA = document.getElementById('sx-ALL');\n",
  "      if (cbA && cbA.checked) return new Set(['ALL']);\n",
  "      const sel = getSelected('sx', sexes);\n",
  "      return sel.size ? sel : new Set(['ALL']);\n",
  "    }\n",
  "\n",
  "    function sexOrderForTooltip(useAllSex, selSx) {\n",
  "      if (useAllSex) return ['ALL'];\n",
  "      const out = [];\n",
  "      ['F','M'].forEach(sx => { if (selSx.has(sx)) out.push(sx); });\n",
  "      // Fallback: if something unexpected is selected, preserve insertion order.\n",
  "      sexes.forEach(sx => { if (selSx.has(sx) && !out.includes(sx)) out.push(sx); });\n",
  "      return out;\n",
  "    }\n",
  "\n",
  "    // Match milestones tooltip: marker glyph in treatment colour + \"DC · F\" (no chip / no mini line art).\n",
  "    function sexSymbolChar(sx) {\n",
  "      if (sx === 'ALL') return '■';\n",
  "      if (sx === 'M') return '◆';\n",
  "      return '●';\n",
  "    }\n",
  "\n",
  "    function fmtTooltipIdent(tr, sx) {\n",
  "      const col = palette[tr] || '#111';\n",
  "      const sym = `<span class=\\\"sym\\\" style=\\\"color:${col}; font-weight: 900\\\">${sexSymbolChar(sx)}</span>`;\n",
  "      return `<span class=\\\"tt-ident-inner\\\">${sym}<span>${tr} · ${sx}</span></span>`;\n",
  "    }\n",
  "\n",
  "    function pctCells(map, tr, sx, selCtOrdered) {\n",
  "      const cells = [];\n",
  "      for (const ct of selCtOrdered) {\n",
  "        const v = map.get(`${tr}|${sx}|${ct}`);\n",
  "        cells.push(Number.isFinite(v) ? `${String(Math.round(v))}%` : '—');\n",
  "      }\n",
  "      return cells;\n",
  "    }\n",
  "\n",
  "    function pctCellsHtml(map, tr, sx, selCtOrdered) {\n",
  "      const nums = [];\n",
  "      const labels = [];\n",
  "      for (const ct of selCtOrdered) {\n",
  "        const v = map.get(`${tr}|${sx}|${ct}`);\n",
  "        if (Number.isFinite(v)) {\n",
  "          nums.push(v);\n",
  "          labels.push(`${String(Math.round(v))}%`);\n",
  "        } else {\n",
  "          nums.push(null);\n",
  "          labels.push('—');\n",
  "        }\n",
  "      }\n",
  "      const finite = nums.filter(x => x !== null && Number.isFinite(x));\n",
  "      const mx = finite.length ? Math.max(...finite.map(v => Math.round(v))) : null;\n",
  "      return labels.map((lab, i) => {\n",
  "        if (lab === '—' || mx === null) return `<td class=\\\"ttpct\\\">${lab}</td>`;\n",
  "        const r = Math.round(nums[i]);\n",
  "        const isMax = r === mx;\n",
  "        if (isMax) return `<td class=\\\"ttpct ttpct-max\\\"><strong>${lab}</strong></td>`;\n",
  "        return `<td class=\\\"ttpct\\\">${lab}</td>`;\n",
  "      }).join('');\n",
  "    }\n",
  "\n",
  "    function rowHasPctData(cells) {\n",
  "      return cells.some(c => c !== '—');\n",
  "    }\n",
  "\n",
  "    function buildTraces({showRibbons, selTr, selSx, selCt}) {\n",
  "      const ribbons = [];\n",
  "      const lines = [];\n",
  "      const useAllSex = selSx.has('ALL');\n",
  "      for (const t of baseTraces) {\n",
  "        const meta = t.meta || {};\n",
  "        if (!selCt.has(meta.call_type)) continue;\n",
  "        const ok = selTr.has(meta.treatment) && (useAllSex ? (meta.sex === 'ALL') : selSx.has(meta.sex));\n",
  "        if (!ok) continue;\n",
  "        if ((meta.kind || '') === 'ribbon') {\n",
  "          if (!showRibbons) continue;\n",
  "          ribbons.push(t);\n",
  "          continue;\n",
  "        }\n",
  "        lines.push(t);\n",
  "      }\n",
  "      return ribbons.concat(lines);\n",
  "    }\n",
  "\n",
  "    const layout = {\n",
  "      autosize: true,\n",
  "      margin: { l: 48, r: 18, t: 10, b: 50 },\n",
  "      showlegend: false,\n",
  "      hovermode: 'x',\n",
  "      hoverdistance: 15,\n",
  "      xaxis: Object.assign(\n",
  "        {\n",
  "          title: 'Age (days)',\n",
  "          showgrid: false,\n",
  "          tickmode: 'linear',\n",
  "          tick0: 30,\n",
  "          dtick: 5,\n",
  "          ticks: 'outside',\n",
  "          ticklen: 6,\n",
  "          // Minor ticks every 1 day (no labels, no extra gridlines)\n",
  "          minor: { dtick: 1, ticks: 'outside', ticklen: 3, showgrid: false }\n",
  "        },\n",
  "        xRange && xRange.every(Number.isFinite)\n",
  "          ? { range: xRange, autorange: false }\n",
  "          : {}\n",
  "      ),\n",
  "      yaxis: {\n",
  "        title: 'Call %',\n",
  "        range: [0, 100],\n",
  "        autorange: false,\n",
  "        // Major ticks/labels/grid every 5%\n",
  "        dtick: 5,\n",
  "        showgrid: true,\n",
  "        gridcolor: 'rgba(0,0,0,0.08)',\n",
  "        ticks: 'outside',\n",
  "        ticklen: 6,\n",
  "        // Minor ticks every 1% (no labels, no extra gridlines)\n",
  "        minor: { dtick: 1, ticks: 'outside', ticklen: 3, showgrid: false }\n",
  "      },\n",
  "      uirevision: 'lock-axes'\n",
  "    };\n",
  "\n",
  "    function render() {\n",
  "      const selTr = getSelected('tr', treatments);\n",
  "      syncSexSelection();\n",
  "      const selSx = getActiveSexes();\n",
  "      const selCt = getSelected('ct', callTypes);\n",
  "      const showRibbons = document.getElementById('show-ribbons').checked;\n",
  "      const traces = buildTraces({showRibbons, selTr, selSx, selCt});\n",
  "      Plotly.react(gd, traces, layout, {responsive: true, displayModeBar: false});\n",
  "    }\n",
  "\n",
  "    document.querySelectorAll('input[type=checkbox]').forEach(cb => cb.addEventListener('change', render));\n",
  "    render();\n",
  "\n",
  "    // Crosshair: horizontal follows cursor y; vertical snaps to hovered age (x).\n",
  "    let rafPending = false;\n",
  "    let lastMouseEvent = null;\n",
  "    let lastHoverX = null;\n",
  "    function setCrosshairShapes(yVal, xVal) {\n",
  "      if (!Number.isFinite(yVal)) return;\n",
  "      const line = { color: 'rgba(120,120,120,0.22)', width: 1 };\n",
  "      const shapes = [\n",
  "        { type: 'line', xref: 'paper', x0: 0, x1: 1, yref: 'y', y0: yVal, y1: yVal, line }\n",
  "      ];\n",
  "      if (Number.isFinite(xVal)) {\n",
  "        shapes.push({ type: 'line', xref: 'x', x0: xVal, x1: xVal, yref: 'paper', y0: 0, y1: 1, line });\n",
  "      }\n",
  "      Plotly.relayout(gd, { shapes });\n",
  "    }\n",
  "    function clearCrosshair() {\n",
  "      lastHoverX = null;\n",
  "      Plotly.relayout(gd, { shapes: [] });\n",
  "    }\n",
  "    function refreshCrosshair() {\n",
  "      if (!lastMouseEvent) return;\n",
  "      const full = gd._fullLayout;\n",
  "      if (!full || !full.yaxis || !full._size) return;\n",
  "      const rect = gd.getBoundingClientRect();\n",
  "      const yClient = lastMouseEvent.clientY - rect.top;\n",
  "      const plotY = yClient - full._size.t;\n",
  "      if (plotY < 0 || plotY > full._size.h) return;\n",
  "      const yVal = full.yaxis.p2l(plotY);\n",
  "      setCrosshairShapes(yVal, lastHoverX);\n",
  "    }\n",
  "    function onMouseMove(e) {\n",
  "      lastMouseEvent = e;\n",
  "      if (rafPending) return;\n",
  "      rafPending = true;\n",
  "      window.requestAnimationFrame(() => {\n",
  "        rafPending = false;\n",
  "        refreshCrosshair();\n",
  "      });\n",
  "    }\n",
  "    gd.addEventListener('mousemove', onMouseMove);\n",
  "    gd.addEventListener('mouseleave', () => { tooltip.style.display = 'none'; clearCrosshair(); });\n",
  "\n",
  "    gd.on('plotly_hover', (ev) => {\n",
  "      const selTr = getSelected('tr', treatments);\n",
  "      syncSexSelection();\n",
  "      const selSx = getActiveSexes();\n",
  "      const selCt = getSelected('ct', callTypes);\n",
  "      const useAllSex = selSx.has('ALL');\n",
  "      const pts = (ev && ev.points) ? ev.points : [];\n",
  "      const e = ev && ev.event ? ev.event : null;\n",
  "      const map = new Map(); // key: `${tr}|${sx}|${ct}` -> percent (0-100)\n",
  "      let ageHeader = null;\n",
  "      lastHoverX = null;\n",
  "      for (const p of pts) {\n",
  "        const t = p.data;\n",
  "        const meta = (t && t.meta) ? t.meta : {};\n",
  "        if ((meta.kind || '') !== 'line') continue;\n",
  "        if (!selCt.has(meta.call_type)) continue;\n",
  "        if (!selTr.has(meta.treatment) || !(useAllSex ? (meta.sex === 'ALL') : selSx.has(meta.sex))) continue;\n",
  "        const age = Number(p.x);\n",
  "        const prop = Number(p.y);\n",
  "        if (!Number.isFinite(age) || !Number.isFinite(prop)) continue;\n",
  "        if (ageHeader === null) ageHeader = Math.round(age);\n",
  "        if (lastHoverX === null) lastHoverX = age;\n",
  "        map.set(`${meta.treatment}|${meta.sex}|${meta.call_type}`, prop);\n",
  "      }\n",
  "      if (!map.size) { tooltip.style.display = 'none'; lastHoverX = null; return; }\n",
  "\n",
  "      const selCtOrdered = callTypes.filter(ct => selCt.has(ct));\n",
  "      const sexOrder = sexOrderForTooltip(useAllSex, selSx);\n",
  "      if (!selCtOrdered.length) {\n",
  "        tooltip.style.display = 'none';\n",
  "        lastHoverX = null;\n",
  "        refreshCrosshair();\n",
  "        return;\n",
  "      }\n",
  "      const colspan = 1 + selCtOrdered.length;\n",
  "      const hdrCells = selCtOrdered.map(ct => `<th class=\\\"ttpcth\\\">${callLabels[ct] || ct}</th>`).join('');\n",
  "      const rows = [];\n",
  "      for (const tr of treatments) {\n",
  "        if (!selTr.has(tr)) continue;\n",
  "        for (const sx of sexOrder) {\n",
  "          const cells = pctCells(map, tr, sx, selCtOrdered);\n",
  "          if (!rowHasPctData(cells)) continue;\n",
  "          rows.push(\n",
  "            `<tr><td class=\\\"tt-ident\\\">${fmtTooltipIdent(tr, sx)}</td>` +\n",
  "              pctCellsHtml(map, tr, sx, selCtOrdered) +\n",
  "              `</tr>`\n",
  "          );\n",
  "        }\n",
  "      }\n",
  "      if (!rows.length) {\n",
  "        tooltip.style.display = 'none';\n",
  "        refreshCrosshair();\n",
  "        return;\n",
  "      }\n",
  "      const ageTop = (ageHeader !== null) ? `<tr><th class=\\\"ttage\\\" colspan=\\\"${colspan}\\\">${ageHeader} d</th></tr>` : '';\n",
  "      const ageRow = `<thead>${ageTop}<tr><th class=\\\"tt-ident-h\\\"></th>${hdrCells}</tr></thead>`;\n",
  "      tooltip.innerHTML = `<table class=\\\"tttbl\\\">${ageRow}<tbody>${rows.join('')}</tbody></table>`;\n",
  "      tooltip.style.display = 'block';\n",
  "      const tipW = tooltip.offsetWidth || 0;\n",
  "      const tipH = tooltip.offsetHeight || 0;\n",
  "      if (e) {\n",
  "        const rect = gd.getBoundingClientRect();\n",
  "        const pad = 8;\n",
  "        let x = e.clientX - rect.left + 14;\n",
  "        let y = e.clientY - rect.top + 14;\n",
  "        x = Math.max(pad, Math.min(x, rect.width - pad - tipW));\n",
  "        y = Math.max(pad, Math.min(y, rect.height - pad - tipH));\n",
  "        tooltip.style.left = x + 'px';\n",
  "        tooltip.style.top = y + 'px';\n",
  "      } else {\n",
  "        tooltip.style.left = '12px';\n",
  "        tooltip.style.top = '12px';\n",
  "      }\n",
  "      refreshCrosshair();\n",
  "    });\n",
  "    gd.on('plotly_unhover', () => {\n",
  "      tooltip.style.display = 'none';\n",
  "      lastHoverX = null;\n",
  "      refreshCrosshair();\n",
  "    });\n",
  "  </script>\n",
  "</body>\n",
  "</html>\n"
)

writeLines(html, con = paths$out_html, useBytes = TRUE)
message("Wrote: ", paths$out_html)

