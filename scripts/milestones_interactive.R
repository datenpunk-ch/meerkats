#!/usr/bin/env Rscript

# Interactive combined milestone plot (custom HTML + Plotly).
# - Single plot (no facets)
# - Treatment = color
# - Sex = marker symbol + line dash
# - Lines = trajectories across SEMI → PEAK → FULL
# - Checkbox filters for Treatment and Sex
# - Optional 95% credible-interval ribbons toggle (off by default)

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
  semi_model = file.path("models", "SEMI_anova.rds"),
  dig_peak_model = file.path("models", "DIG_PEAK_anova.rds"),
  full_model = file.path("models", "FULL_anova.rds"),
  semi_data = file.path("data", "SEMI_data.rds"),
  dig_peak_data = file.path("data", "PEAK_data.rds"),
  full_data = file.path("data", "FULL_data.rds"),
  out_html = file.path("figures", "Milestones_Combined.html")
)

dir.create(paths$figures_dir, showWarnings = FALSE, recursive = TRUE)

read_fit <- function(p) {
  if (!file.exists(p)) stop("Missing model file: ", p, "\n(Models are ignored by git; ensure it exists locally.)")
  readRDS(p)
}

read_dat <- function(p) {
  if (!file.exists(p)) stop("Missing data file: ", p)
  readRDS(p)
}

summarize_milestone <- function(fit, dat, milestone_label) {
  newdata <- tidyr::expand_grid(
    TREATMENT = levels(as.factor(dat$TREATMENT)),
    SEX = levels(as.factor(dat$SEX))
  )

  fit %>%
    tidybayes::epred_draws(newdata = newdata, re_formula = NA, robust = TRUE) %>%
    rename(Milestone_Age = .epred) %>%
    mutate(
      Milestone_label = milestone_label,
      SEX = as.factor(SEX),
      TREATMENT = as.factor(TREATMENT)
    )
}

semi_fit <- read_fit(paths$semi_model)
dig_peak_fit <- read_fit(paths$dig_peak_model)
full_fit <- read_fit(paths$full_model)

semi_dat <- read_dat(paths$semi_data)
dig_peak_dat <- read_dat(paths$dig_peak_data)
full_dat <- read_dat(paths$full_data)

milestone_order <- c("SEMI", "PEAK", "FULL")

draws_all <- bind_rows(
  summarize_milestone(semi_fit, semi_dat, "SEMI"),
  summarize_milestone(dig_peak_fit, dig_peak_dat, "PEAK"),
  summarize_milestone(full_fit, full_dat, "FULL")
) %>%
  mutate(Milestone_label = factor(Milestone_label, levels = milestone_order))

summ <- draws_all %>%
  group_by(TREATMENT, SEX, Milestone_label) %>%
  tidybayes::median_qi(Milestone_Age, .width = 0.95) %>%
  ungroup()

# Also compute a sex-marginalised summary ("no sex interaction" view) by averaging
# over SEX within each posterior draw, then summarising uncertainty.
summ_allsex <- draws_all %>%
  group_by(.draw, TREATMENT, Milestone_label) %>%
  summarise(Milestone_Age = mean(Milestone_Age, na.rm = TRUE), .groups = "drop") %>%
  group_by(TREATMENT, Milestone_label) %>%
  tidybayes::median_qi(Milestone_Age, .width = 0.95) %>%
  ungroup() %>%
  mutate(SEX = factor("ALL", levels = "ALL"))

sex_levels <- levels(droplevels(summ$SEX))
treat_levels <- levels(droplevels(summ$TREATMENT))

# Prefer a stable, meaningful treatment order when applicable.
preferred_treat_order <- c("DC", "SC", "DT")
if (all(preferred_treat_order %in% treat_levels)) {
  treat_levels <- c(preferred_treat_order, setdiff(treat_levels, preferred_treat_order))
}

# Lock y-axis range so toggling ribbons doesn't change axes.
y_min <- floor(min(summ$.lower, na.rm = TRUE))
y_max <- ceiling(max(summ$.upper, na.rm = TRUE))

# Okabe-Ito-ish (works well on white), assigned in treatment order.
make_palette <- function(levels) {
  # Standard Okabe–Ito palette (common ordering used in many gg* helpers)
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
  # => DC = skyblue, SC = orange, DT = bluishgreen
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

make_traces <- function(summ, summ_allsex, sex_levels, treat_levels) {
  traces <- list()

  # Symbol/dash maps for sex (fallbacks if levels differ)
  sex_symbol <- setNames(rep("circle", length(sex_levels)), sex_levels)
  sex_dash <- setNames(rep("solid", length(sex_levels)), sex_levels)
  if ("F" %in% sex_levels) {
    sex_symbol[["F"]] <- "circle"
    sex_dash[["F"]] <- "solid"
  }
  if ("M" %in% sex_levels) {
    sex_symbol[["M"]] <- "diamond"
    sex_dash[["M"]] <- "dash"
  }
  # Fallback symbol for the sex-marginalised view.
  sex_symbol[["ALL"]] <- "square"
  sex_dash[["ALL"]] <- "solid"

  for (sx in sex_levels) {
    for (tr in treat_levels) {
      df <- summ %>% filter(SEX == sx, TREATMENT == tr) %>% arrange(Milestone_label)
      if (nrow(df) == 0) next

      col <- palette[[as.character(tr)]]
      if (is.null(col)) col <- "#333333"

      x0 <- as.character(df$Milestone_label)
      y0 <- df$Milestone_Age
      lo0 <- df$.lower
      hi0 <- df$.upper

      legend_name <- paste0(tr, " — ", sx)

      # Ribbon trace (hidden by default; can be toggled via checkbox)
      traces[[length(traces) + 1]] <- list(
        x = c(x0, rev(x0)),
        y = c(lo0, rev(hi0)),
        type = "scatter",
        mode = "lines",
        fill = "toself",
        fillcolor = paste0(col, "22"), # light fill
        line = list(color = "rgba(0,0,0,0)"),
        hoverinfo = "skip",
        name = paste0(legend_name, " 95% CI"),
        meta = list(treatment = as.character(tr), sex = as.character(sx), kind = "ribbon")
      )

      # Median line+points trace
      traces[[length(traces) + 1]] <- list(
        x = x0,
        y = y0,
        type = "scatter",
        mode = "lines+markers",
        line = list(color = col, width = 3, dash = sex_dash[[sx]]),
        marker = list(color = col, size = 8, symbol = sex_symbol[[sx]]),
        name = legend_name,
        # Disable Plotly's native hover box; we render our own tooltip.
        hoverinfo = "none",
        hovertemplate = "",
        meta = list(treatment = as.character(tr), sex = as.character(sx), kind = "line")
      )
    }
  }

  # Add sex-marginalised traces with meta.sex = "ALL"
  for (tr in treat_levels) {
    df <- summ_allsex %>% filter(TREATMENT == tr) %>% arrange(Milestone_label)
    if (nrow(df) == 0) next

    col <- palette[[as.character(tr)]]
    if (is.null(col)) col <- "#333333"

    x0 <- as.character(df$Milestone_label)
    y0 <- df$Milestone_Age
    lo0 <- df$.lower
    hi0 <- df$.upper

    legend_name <- paste0(tr, " — ALL")

    traces[[length(traces) + 1]] <- list(
      x = c(x0, rev(x0)),
      y = c(lo0, rev(hi0)),
      type = "scatter",
      mode = "lines",
      fill = "toself",
      fillcolor = paste0(col, "22"),
      line = list(color = "rgba(0,0,0,0)"),
      hoverinfo = "skip",
      name = paste0(legend_name, " 95% CI"),
      meta = list(treatment = as.character(tr), sex = "ALL", kind = "ribbon")
    )

    traces[[length(traces) + 1]] <- list(
      x = x0,
      y = y0,
      type = "scatter",
      mode = "lines+markers",
      line = list(color = col, width = 3, dash = sex_dash[["ALL"]]),
      marker = list(color = col, size = 8, symbol = sex_symbol[["ALL"]]),
      name = legend_name,
      hoverinfo = "none",
      hovertemplate = "",
      meta = list(treatment = as.character(tr), sex = "ALL", kind = "line")
    )
  }

  traces
}

traces <- make_traces(summ, summ_allsex, sex_levels, treat_levels)

layout <- list(
  autosize = TRUE,
  margin = list(l = 60, r = 25, t = 45, b = 60),
  showlegend = FALSE,
  height = 740,
  yaxis = list(
    title = "Milestone age (days)",
    range = c(y_min, y_max),
    autorange = FALSE,
    # Major ticks/labels/grid every 5 days
    dtick = 5,
    showgrid = TRUE,
    # Show the axis ticks
    ticks = "outside",
    ticklen = 6,
    # Minor ticks every day (no labels)
    minor = list(dtick = 1, ticks = "outside", ticklen = 3, showgrid = FALSE)
  ),
  # Categorical axis; padding/range will be derived in JS (generic output).
  xaxis = list(
    title = "Milestone (ontogeny)",
    type = "category",
    categoryorder = "array",
    categoryarray = milestone_order
  ),
  hovermode = "x",
  hoverdistance = 15,
  uirevision = "lock-axes"
)

payload <- list(
  traces = traces,
  layout = layout,
  # Serialize as an object in JSON (name -> hex color).
  palette = as.list(palette),
  sex_levels = sex_levels,
  treat_levels = treat_levels,
  milestone_order = milestone_order
)

payload_json <- jsonlite::toJSON(payload, auto_unbox = TRUE, digits = 12)

html <- paste0(
  "<!doctype html>\n",
  "<html lang=\"en\">\n",
  "<head>\n",
  "  <meta charset=\"utf-8\" />\n",
  "  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\" />\n",
  "  <title>Milestones — interactive</title>\n",
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
  "    /* If the (auto) disabled checkbox is ON, make that state obvious */\n",
  "    label:has(input:disabled:checked) { opacity: 1; }\n",
  "    label:has(input:disabled:checked) input { accent-color: #000; }\n",
  "    label:has(input:disabled:checked) .sexsample { outline: 2px solid #000; outline-offset: 2px; border-radius: 3px; }\n",
  "    label:has(input:disabled:checked) .sexline { border-top-color: #000; }\n",
  "    label:has(input:disabled:checked) .sexmark { background: #000; }\n",
  "    #plot { width: 100%; height: 740px; }\n",
  "    #tooltip { position: absolute; z-index: 20; padding: 8px 10px; border: 1px solid #ddd; background: rgba(255,255,255,0.95); font-size: 12px; color: #111; display: none; max-width: none; width: max-content; pointer-events: none; }\n",
  "    #tooltip .row { display: flex; align-items: baseline; justify-content: space-between; gap: 10px; white-space: nowrap; }\n",
  "    #tooltip .left { white-space: nowrap; display: inline-flex; align-items: baseline; gap: 6px; }\n",
  "    #tooltip .sym { width: 1.1em; display: inline-block; text-align: center; }\n",
  "    #tooltip .age { font-variant-numeric: tabular-nums; }\n",
  "    .hint { color: #555; font-size: 12px; margin: 6px 0 0; }\n",
  "    .key { display: grid; grid-template-columns: 90px 1fr; column-gap: 12px; row-gap: 6px; align-items: center; }\n",
  "    .keyitems { display: flex; flex-wrap: wrap; gap: 10px 16px; align-items: center; }\n",
  "    .swatch { display: inline-flex; align-items: center; gap: 8px; }\n",
  "    .chip { width: 14px; height: 14px; background: #000; display: inline-block; }\n",
  "    .sexsample { position: relative; width: 34px; height: 14px; display: inline-block; }\n",
  "    .sexline { position: absolute; left: 0; right: 0; top: 50%; transform: translateY(-50%); border-top: 3px solid #111; }\n",
  "    .sexline--dash { border-top-style: dashed; }\n",
  "    .sexmark { position: absolute; left: 50%; top: 50%; width: 10px; height: 10px; background: #111; transform: translate(-50%, -50%); }\n",
  "    .sexmark--diamond { transform: translate(-50%, -50%) rotate(45deg); }\n",
  "  </style>\n",
  "</head>\n",
  "<body>\n",
  "  <div class=\"wrap\">\n",
  "    <div class=\"controls\">\n",
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
  "      <fieldset>\n",
  "        <legend>Uncertainty</legend>\n",
  "        <label><input type=\"checkbox\" id=\"show-ribbons\" /> Show 95% CI ribbons</label>\n",
  "      </fieldset>\n",
  "    </div>\n",
  "    <div style=\"position: relative\">\n",
  "      <div id=\"plot\"></div>\n",
  "      <div id=\"tooltip\"></div>\n",
  "    </div>\n",
  "    <p class=\"hint\">Tip: uncheck treatments/sexes to compare specific trajectories. Turn on ribbons only when you’ve reduced clutter.</p>\n",
  "  </div>\n",
  "  <script>\n",
  "    const payload = ", payload_json, ";\n",
  "    const gd = document.getElementById('plot');\n",
  "    const treatments = payload.treat_levels;\n",
  "    const sexes = payload.sex_levels;\n",
  "    const milestoneOrder = payload.milestone_order;\n",
  "    const palette = payload.palette || {};\n",
  "    const tooltip = document.getElementById('tooltip');\n",
  "\n",
  "    const layout = payload.layout;\n",
  "\n",
  "    // Build checkboxes dynamically (generic; no hardcoded DC/SC/DT etc.)\n",
  "    function safeId(v) {\n",
  "      return String(v).replace(/[^a-zA-Z0-9_-]/g, '_');\n",
  "    }\n",
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
  "    const trBox = document.getElementById('treatment-boxes');\n",
  "    treatments.forEach(tr => {\n",
  "      const col = palette[tr] || '#111';\n",
  "      addCheckboxRow(trBox, `tr-${safeId(tr)}`, tr, `<span class=\\\"chip\\\" style=\\\"background:${col}\\\"></span> ${tr}`, true);\n",
  "    });\n",
  "\n",
  "    const sxBox = document.getElementById('sex-boxes');\n",
  "    function sexSample(label) {\n",
  "      if (label === 'M') return `<span class=\\\"sexsample\\\"><span class=\\\"sexline sexline--dash\\\"></span><span class=\\\"sexmark sexmark--diamond\\\"></span></span> M`;\n",
  "      if (label === 'F') return `<span class=\\\"sexsample\\\"><span class=\\\"sexline\\\"></span><span class=\\\"sexmark\\\" style=\\\"border-radius:50%\\\"></span></span> F`;\n",
  "      return `<span class=\\\"sexsample\\\"><span class=\\\"sexline\\\"></span><span class=\\\"sexmark\\\"></span></span> Averaged`;\n",
  "    }\n",
  "    sexes.forEach(sx => {\n",
  "      addCheckboxRow(sxBox, `sx-${safeId(sx)}`, sx, sexSample(sx), true, false);\n",
  "    });\n",
  "    // Averaged is auto-managed (indicator), not user-clickable.\n",
  "    addCheckboxRow(sxBox, 'sx-ALL', 'ALL', sexSample('ALL'), false, true);\n",
  "\n",
  "    // Compute and lock y-range from all trace y-values (keeps axes stable when toggling ribbons)\n",
  "    function computeYRange() {\n",
  "      let lo = Infinity, hi = -Infinity;\n",
  "      for (const t of payload.traces) {\n",
  "        if (!t || !Array.isArray(t.y)) continue;\n",
  "        for (const v of t.y) {\n",
  "          const n = Number(v);\n",
  "          if (!Number.isFinite(n)) continue;\n",
  "          lo = Math.min(lo, n);\n",
  "          hi = Math.max(hi, n);\n",
  "        }\n",
  "      }\n",
  "      if (!Number.isFinite(lo) || !Number.isFinite(hi)) return null;\n",
  "      lo = Math.floor(lo);\n",
  "      hi = Math.ceil(hi);\n",
  "      return [lo, hi];\n",
  "    }\n",
  "    const yRange = computeYRange();\n",
  "    if (yRange) {\n",
  "      layout.yaxis = Object.assign({}, layout.yaxis || {}, {\n",
  "        range: yRange,\n",
  "        autorange: false,\n",
  "        dtick: 5,\n",
  "        showgrid: true,\n",
  "        ticks: 'outside',\n",
  "        ticklen: 6,\n",
  "        minor: { dtick: 1, ticks: 'outside', ticklen: 3, showgrid: false }\n",
  "      });\n",
  "    }\n",
  "    // Add a bit of categorical padding left/right of first/last milestone.\n",
  "    if (Array.isArray(milestoneOrder) && milestoneOrder.length) {\n",
  "      const n = milestoneOrder.length;\n",
  "      // For category axes, Plotly uses category indices for range.\n",
  "      layout.xaxis = Object.assign({}, layout.xaxis || {}, { range: [-0.4, (n - 1) + 0.4] });\n",
  "    }\n",
  "\n",
  "    const baseTraces = payload.traces;\n",
  "    const baseLayout = layout;\n",
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
  "\n",
  "      // If user selects any sex, disable Averaged.\n",
  "      if (anySexChecked) {\n",
  "        cbA.checked = false;\n",
  "        return;\n",
  "      }\n",
  "\n",
  "      // If none are selected, automatically switch to Averaged.\n",
  "      cbA.checked = true;\n",
  "    }\n",
  "\n",
  "    function getActiveSexes() {\n",
  "      // Returns Set of allowed meta.sex values based on UI.\n",
  "      const cbA = document.getElementById('sx-ALL');\n",
  "      if (cbA && cbA.checked) return new Set(['ALL']);\n",
  "      const sel = getSelected('sx', sexes);\n",
  "      // If empty, syncSexSelection should have turned on Averaged, but fallback anyway.\n",
  "      return sel.size ? sel : new Set(['ALL']);\n",
  "    }\n",
  "\n",
  "    function buildTraces({showRibbons, selTr, selSx}) {\n",
  "      const ribbons = [];\n",
  "      const lines = [];\n",
  "      const useAllSex = selSx.has('ALL');\n",
  "      for (const t of baseTraces) {\n",
  "        const meta = t.meta || {};\n",
  "        const ok = selTr.has(meta.treatment) && (useAllSex ? (meta.sex === 'ALL') : selSx.has(meta.sex));\n",
  "        if (!ok) continue;\n",
  "        if ((meta.kind || '') === 'ribbon') {\n",
  "          // keep ribbons off unless explicitly enabled\n",
  "          if (!showRibbons) {\n",
  "            // hide ribbon entirely in this mode\n",
  "            continue;\n",
  "          }\n",
  "          ribbons.push(t);\n",
  "          continue;\n",
  "        }\n",
  "        lines.push(t);\n",
  "      }\n",
  "      // Ensure uncertainty ribbons are always behind lines.\n",
  "      return ribbons.concat(lines);\n",
  "    }\n",
  "\n",
  "    function render() {\n",
  "      const selTr = getSelected('tr', treatments);\n",
  "      syncSexSelection();\n",
  "      const selSx = getActiveSexes();\n",
  "      const showRibbons = document.getElementById('show-ribbons').checked;\n",
  "      const traces = buildTraces({showRibbons, selTr, selSx});\n",
  "      Plotly.react(gd, traces, baseLayout, {responsive: true, displayModeBar: false});\n",
  "    }\n",
  "\n",
  "    Plotly.newPlot(gd, buildTraces({showRibbons:false, selTr:new Set(treatments), selSx:new Set(sexes)}), baseLayout, {responsive: true, displayModeBar: false});\n",
  "\n",
  "    document.querySelectorAll('input[type=checkbox]').forEach(cb => {\n",
  "      cb.addEventListener('change', render);\n",
  "    });\n",
  "    render();\n",
  "\n",
  "    function traceColor(t) {\n",
  "      return (((t||{}).line||{}).color) || (((t||{}).marker||{}).color) || '#111';\n",
  "    }\n",
  "\n",
  "    function fmtRow(color, tr, sx, age) {\n",
  "      const sexSymbol = (sx === 'ALL') ? '■' : ((sx === 'M') ? '◆' : '●');\n",
  "      const sym = `<span class=\\\"sym\\\" style=\\\"color:${color}; font-weight: 900\\\">${sexSymbol}</span>`;\n",
  "      const label = `<span>${tr} · ${sx}</span>`;\n",
  "      const ageInt = Math.round(age);\n",
  "      return `<div class=\\\"row\\\"><span class=\\\"left\\\">${sym}${label}</span><span class=\\\"age\\\">${ageInt} d</span></div>`;\n",
  "    }\n",
  "\n",
  "    // Continuous crosshair line following the cursor's y position.\n",
  "    let rafPending = false;\n",
  "    let lastMouseEvent = null;\n",
  "    function setCrosshair(yVal) {\n",
  "      if (!Number.isFinite(yVal)) return;\n",
  "      Plotly.relayout(gd, {\n",
  "        shapes: [{\n",
  "          type: 'line',\n",
  "          xref: 'paper', x0: 0, x1: 1,\n",
  "          yref: 'y', y0: yVal, y1: yVal,\n",
  "          line: { color: 'rgba(120,120,120,0.25)', width: 1 }\n",
  "        }]\n",
  "      });\n",
  "    }\n",
  "\n",
  "    function clearCrosshair() {\n",
  "      Plotly.relayout(gd, { shapes: [] });\n",
  "    }\n",
  "\n",
  "    function onMouseMove(e) {\n",
  "      lastMouseEvent = e;\n",
  "      if (rafPending) return;\n",
  "      rafPending = true;\n",
  "      window.requestAnimationFrame(() => {\n",
  "        rafPending = false;\n",
  "        if (!lastMouseEvent) return;\n",
  "        const full = gd._fullLayout;\n",
  "        if (!full || !full.yaxis || !full._size) return;\n",
  "        const rect = gd.getBoundingClientRect();\n",
  "        const yClient = lastMouseEvent.clientY - rect.top;\n",
  "        const plotY = yClient - full._size.t;\n",
  "        if (plotY < 0 || plotY > full._size.h) return;\n",
  "        const yVal = full.yaxis.p2l(plotY);\n",
  "        setCrosshair(yVal);\n",
  "      });\n",
  "    }\n",
  "\n",
  "    gd.addEventListener('mousemove', onMouseMove);\n",
  "    gd.addEventListener('mouseleave', () => { clearCrosshair(); });\n",
  "\n",
  "    gd.on('plotly_hover', (ev) => {\n",
  "      const selTr = getSelected('tr', treatments);\n",
  "      syncSexSelection();\n",
  "      const selSx = getActiveSexes();\n",
  "      const useAllSex = selSx.has('ALL');\n",
  "      const pts = (ev && ev.points) ? ev.points : [];\n",
  "      const e = ev && ev.event ? ev.event : null;\n",
  "      const rows = [];\n",
  "      let yHover = null;\n",
  "      for (const p of pts) {\n",
  "        const t = p.data;\n",
  "        const meta = (t && t.meta) ? t.meta : {};\n",
  "        if ((meta.kind || '') !== 'line') continue;\n",
  "        if (!selTr.has(meta.treatment) || !(useAllSex ? (meta.sex === 'ALL') : selSx.has(meta.sex))) continue;\n",
  "        if (t && (t.visible === false || t.visible === 'legendonly')) continue;\n",
  "        const age = Number(p.y);\n",
  "        if (!Number.isFinite(age)) continue;\n",
  "        if (yHover === null) yHover = age;\n",
  "        rows.push({ age, html: fmtRow(traceColor(t), meta.treatment, meta.sex, age) });\n",
  "      }\n",
  "      rows.sort((a,b) => b.age - a.age);\n",
  "      if (!rows.length) { tooltip.style.display = 'none'; return; }\n",
  "      tooltip.innerHTML = rows.map(r => r.html).join('');\n",
  "      // Make visible briefly so we can measure size for clamping.\n",
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
  "    });\n",
  "    gd.on('plotly_unhover', () => {\n",
  "      tooltip.style.display = 'none';\n",
  "    });\n",
  "  </script>\n",
  "</body>\n",
  "</html>\n"
)

writeLines(html, con = paths$out_html, useBytes = TRUE)
message("Wrote: ", paths$out_html)

