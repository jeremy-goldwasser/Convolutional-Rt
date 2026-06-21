# code/

Analysis code for *Fast, Frequentist Estimation of Epidemic Reproduction Numbers*.
The ConvRt method is implemented in the [ConvRt](https://github.com/jeremy-goldwasser/ConvRt)
R package; install it once and the scripts here call it directly.

## Dependencies

```r
# Install the ConvRt package (one-time)
# devtools::install_github("jeremy-goldwasser/ConvRt")

# Install remaining R dependencies
Rscript code/install_packages.R
```

## Layout

```
code/
├── _paths.R               Auto-detects the repo root from any working directory;
│                          defines data_dir, results_dir, figures_dir.
├── install_packages.R     Installs all CRAN dependencies.
├── examples/
│   └── demo.R             Self-contained worked example: simulate from the model,
│                          fit ConvRt, plot the recovered R(t).
└── analyses/              One subdirectory per dataset / comparison.
    ├── covidestim/        COVID-19 comparison against covidestim (CA, LA, counties).
    ├── flu/               U.S. flu hospitalization analysis (seasons 2022/23, 2023/24).
    ├── rtestim/           Benchmark against rtestim on renewal- and SEIR-simulated data.
    ├── sim_flu/           Flu-based SEIR simulation benchmark (retro + real-time).
    ├── tf/                Trend-filtering variant comparison.
    └── weekly/            Weekly-aggregated estimation proof-of-concept.
```

Each subdirectory contains a `_common.R` that loads ConvRt, detects the repo
root, and defines dataset-specific helpers shared across that subdirectory's
scripts. Analysis scripts are standalone: run them with `Rscript` from any
directory and they write results to `results/` and figures to `figures/`.

## Model summary

ConvRt estimates the time-varying reproduction number R(t) from observed case or
hospitalization counts. The observation model is:

```
Y_t ~ Poisson( rho_t * omega_{t mod 7} * sum_{s<t} R_s * X_{t,s} )
X_{t,s} = Lambda_s * pi_EY[t-s]
Lambda_s = sum_{k>=1} g_k * X_exposure[s-k]
```

- `g` — generation interval PMF (supplied by the user; e.g. a discrete Gamma).
- `pi_EY` — exposure-to-observation reporting delay PMF.
- `X_exposure` — latent exposure incidence, recovered by deconvolving observations
  against `pi_EY`.
- `R(t)` — parameterized as a natural cubic smoothing spline; penalty strength
  `lambda` selected by K-fold cross-validation (Poisson deviance, 1se rule).
- Real-time fits add a CDF-tapered tail penalty (strength `gamma` tuned by forward
  validation) and split-conformal prediction intervals.

The only epidemiological inputs required are a generation interval and a reporting
delay — both as single PMFs.
