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
│   ├── _paths.R               Auto-detects repo root; defines code_dir, data_dir,
│   │                          results_dir, figures_dir — scripts are runnable from
│   │                          any working directory
│   ├── install_packages.R     Install all R dependencies
│   ├── R/                     ConvRt method library (modularized)
│   │   ├── wrappers.R         High-level entry points: fit_convrt_retrospective,
│   │   │                      fit_convrt_realtime
│   │   ├── design.R           Design-matrix construction (smoothing-spline basis)
│   │   ├── solve.R            Penalized IRLS + KKT solver
│   │   ├── lambda_select.R    K-fold CV lambda selection (5-fold, 1se, Poisson deviance)
│   │   ├── extract.R          R(t) extraction with Laplace pointwise CIs and simbands
│   │   ├── deconvolution.R    Deconvolve-then-convolve prevalence imputation (GCV-tuned)
│   │   ├── taper.R            Tapered penalty for real-time fits; forward-validation
│   │   │                      gamma tuning
│   │   ├── delays.R           Delay distribution utilities
│   │   ├── realtime.R         Real-time fitting helpers
│   │   ├── renewal.R          Renewal equation utilities
│   │   ├── conformal.R        Conformal prediction intervals
│   │   ├── tf.R               Trend-filtering utilities
│   │   └── weekly.R           Weekly aggregation helpers
│   ├── analyses/              Analysis scripts organized by dataset/comparison
│   │   ├── covidestim/        CA/LA/county retrospective fits; GI sensitivity analysis
│   │   ├── flu/               U.S. flu hospitalization pipelines (seasons 1 and 2)
│   │   ├── rtestim/           Comparison against rtestim on renewal/SEIR scenarios
│   │   ├── sim_flu/           Flu-based simulation benchmark (SEIR + hospitalization)
│   │   ├── tf/                Trend-filtering vs. smoothing-spline comparison
│   │   └── weekly/            Weekly-data estimation experiments
│   └── examples/
│       └── demo.R             Minimal worked example
├── Data/                      Input data: observed counts, ground-truth R(t), GI PMFs
│   ├── real/                  Real-data inputs (covidestim, flu hospitalizations)
│   └── sim/                   Simulated-data inputs
├── results/                   Generated outputs: R(t) fits, CV diagnostics (.rds / .csv)
│   ├── real/
│   └── sim/
└── figures/                   Generated figures (.pdf)
    ├── real/
    ├── sim/
    └── extra_methods/
```

## Quickstart

```r
# Install R dependencies (one-time)
Rscript code/install_packages.R

# Source the ConvRt library
source("code/R/wrappers.R")

# Retrospective fit
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
