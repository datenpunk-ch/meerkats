# Meerkats Against Patriarchy

This repository contains the project materials for **“Meerkats Against Patriarchy”** (data-journalism / thesis-based portfolio project).

## What’s in here
- **Project page (HTML)**: `meerkats-against-patriarchy.html`
- **Project text (Markdown)**: `meerkats-against-patriarchy.md`
- **Figure outputs**: `figures/` (includes `figures/CallProportions.png`)
- **Photos**: `pics/`
- **Data**: `data/` (RDS files used by the plotting script)
- **Scripts**: `scripts/`

## Preview locally
From the repository root:

```powershell
python -m http.server 8000
```

Then open:
- `http://localhost:8000/meerkats-against-patriarchy.html`

## Regenerate the call-proportions figure
The page includes `figures/CallProportions.png`, generated from a fitted `brms` model and `PROP_data`.

### 1) Requirements
- R installed
- A working Stan backend (this project used **CmdStanR** originally)

Install required R packages:

```r
install.packages(c("brms", "dplyr", "tidyr", "tidybayes", "ggplot2", "ggdist"))
```

### 2) Provide the fitted model locally
Large fitted model files are ignored by git:
- `models/*.rds` is in `.gitignore`

To regenerate the figure you need the model file available locally at:
- `models/B_prop.rds`

### 3) Run the script
From the repo root:

```powershell
Rscript scripts/b_prop_interactive.R
```

This will write:
- `figures/CallProportions.png`

