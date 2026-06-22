# Convolutional-Rt

This repository contains the analysis code for the paper:

**Fast, Frequentist Estimation of Epidemic Reproduction Numbers**

The paper introduces **ConvRt**, a method for estimating the time-varying reproduction number R(t) from observed case or hospitalization counts. ConvRt uses a smoothing-spline parameterization of R(t) with penalized Poisson regression, cross-validated lambda selection, and a deconvolution step that imputes infectious prevalence from observations. The method supports both retrospective (end-of-season) and real-time (right-censored) estimation.

The ConvRt R package is available at: [https://github.com/jeremy-goldwasser/ConvRt](https://github.com/jeremy-goldwasser/ConvRt)

---

## Repository structure

```
Convolutional-Rt/
├── code/
│   ├── _paths.R               Auto-detects repo root; defines data_dir, results_dir, figures_dir
│   ├── install_packages.R     Install all R dependencies
│   ├── examples/
│   │   └── demo.R             Self-contained worked example
│   └── analyses/              One subdirectory per dataset / comparison
│       ├── covidestim/        COVID-19 comparison against covidestim (CA, LA, counties)
│       ├── flu/               U.S. flu hospitalization analysis (seasons 2022/23, 2023/24)
│       ├── rtestim/           Benchmark against rtestim on renewal- and SEIR-simulated data
│       ├── sim_flu/           Flu-based SEIR simulation benchmark (retro + real-time)
│       ├── tf/                Trend-filtering variant comparison
│       └── weekly/            Weekly-aggregated estimation proof-of-concept
├── Data/                      Input data: observed counts, ground-truth R(t), GI PMFs
├── results/                   Generated outputs: R(t) fits, CV diagnostics (.rds / .csv)
└── figures/                   Generated figures (.pdf)
```

## Quickstart

```r
# Install the ConvRt package (one-time)
devtools::install_github("jeremy-goldwasser/ConvRt")

# Install remaining R dependencies (one-time)
Rscript code/install_packages.R

# Retrospective fit
library(ConvRt)
res <- fit_convrt_retrospective(
  counts                 = counts,   # observed case/hospitalization counts
  pi_EY                  = pi_EY,   # generation-interval PMF
  mean_infectious_period = 3.5
)
res$rt_df  # data.frame with day, Rt_mean, Rt_lo, Rt_hi

# Real-time fit
res <- fit_convrt_realtime(
  counts = counts,
  pi_EY  = pi_EY,
  lam    = lam                       # or pass lam_grid for fresh CV
)
res$rt_df
```

See `code/examples/demo.R` for a self-contained worked example.

## Reproducing paper figures

Each subdirectory of `code/analyses/` is self-contained. Scripts auto-detect the repo root via `code/_paths.R` and write results/figures to the corresponding paths under `results/` and `figures/`.
