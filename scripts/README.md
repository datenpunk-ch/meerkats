# Model figure (B_prop)

This project includes a single script to reuse your fitted `brms` model (`models/B_prop.rds`) as a static figure referenced by:
- `meerkats-against-patriarchy.html`

## 1) Install R + packages
You need R available on your machine.

In an R session:

```r
install.packages(c("brms", "dplyr", "tidyr", "tidybayes", "ggplot2", "ggdist"))
```

> `brms` also requires a working Stan backend (CmdStanR or rstan). Use your existing thesis setup.

## 2) Export the interactive model figure
From the project root (the folder that contains `data/` and `models/`), run:

```powershell
Rscript scripts/b_prop_interactive.R
```

This creates:
- `figures/CallProportions.png`

## 3) View the HTML
Open:
- `meerkats-against-patriarchy.html`

The figure will appear automatically once the PNG exists.

### Tip: open via a local web server
This avoids path issues and makes it easy to preview assets.

From the project root you can run (Python example):

```powershell
python -m http.server 8000
```

Then visit:
- `http://localhost:8000/meerkats-against-patriarchy.html`

