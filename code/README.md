# girt — Generation-Interval R_t

A regularized identity-link Poisson estimator of the time-varying reproduction
number R_t that fits a smoothing-spline R_t through the **renewal equation** and
predicts the **observed counts**. It is the non-mechanistic sibling of `mechrt`:
same observation model and same penalized-spline engine, but the per-day
transmission force is generated from a single **generation interval `g`** instead
of compartmental (E→I latent, I→R infectious) distributions.

## Model

```
Y_t ~ Poisson( mu_t ),   mu_t = rho_t * omega_(t mod 7) * sum_{s<t} R_s * X_{t,s}
X_{t,s} = Lambda_s * pi_EY[t-s],     Lambda_s = sum_{k>=1} g_k * X_exposure[s-k]
```

- `X_exposure` — exposure (infection) incidence, recovered by deconvolving the
  observations against the **reporting delay** `pi_EY` only (no SEIR kernel).
- `Lambda_s` — the renewal force: trailing exposures convolved with the
  generation interval `g`. The renewal equation is `E[X_s] = R_s * Lambda_s`.
- `pi_EY` — exposure→observed-outcome (e.g. E→hospitalization) reporting delay.
- `R_t = (B theta)_t` — a natural cubic smoothing spline (integrated 2nd-deriv
  penalty), `lambda` by K-fold CV. Real-time adds a CDF-tapered tail penalty with
  strength `gamma` tuned by forward validation, and split-conformal CIs.

The only epidemiological inputs are a generation interval and a reporting delay —
both are single PMFs (literature, external estimate, or, optionally,
`gi_from_compartmental()` for the GI implied by an SEIR model).

## Relationship to `mechrt`

With the SEIR-implied GI `g_k = zeta^EI_k / mu^IR` (= `P_still_I` normalized),
girt is **algebraically identical** to a mechanistic MechRt fit
(appendix3.tex, Prop. compartmental-generation): the renewal force `Lambda_s`
equals `prevalence_s / mu^IR`, so the design, the CV-λ, the FV-γ, and the fitted
R_t coincide to machine precision. girt only diverges from MechRt when you feed a
generation interval that the SEIR delays cannot reproduce. (This equivalence was
verified to ~1e-14 on flu3 retro + real-time; those mechrt-comparison scripts
live on the `main` branch, which still has the mechanistic codebase. This branch
keeps only girt.)

## Usage

```r
source("code/girt.R")

# generation interval (here: implied by an SEIR latent + infectious period;
# or supply any literature/external GI pmf)
g <- gi_from_compartmental(pi_lat = gi_discrete_gamma_delay(2.0, 1.2)$pmf,
                           pi_IR  = gi_discrete_gamma_delay(2.75, 1.0)$pmf)$g

# retrospective
retro <- fit_girt_retrospective(obs_inc = y, dates = dates, g = g,
                                mean_EY = 5.7, sd_EY = 2.3, severity = 0.015,
                                first_rt_date = as.Date("2022-07-01"),
                                likelihood_start_date = as.Date("2022-07-22"))
retro$rt_df    # date, day, Rt_mean, Rt_lo, Rt_hi

# real-time (right-censored vintage)
rt <- fit_girt_realtime(obs_inc = y_vintage, dates = dates_vintage, g = g,
                        mean_EY = 5.7, sd_EY = 2.3, severity = 0.015,
                        first_rt_date = as.Date("2022-07-01"),
                        likelihood_start_date = as.Date("2022-07-22"))
```

## Layout

```
girt.R              loader (source this)
R/
  delays.R          delay pmfs + gi_from_compartmental
  deconvolution.R   gi_deconvolve_exposures (reporting delay only)   [girt-native]
  renewal.R         gi_renewal_force                                  [girt-native]
  design.R          build_gi_design (renewal design)                 [girt-native]
  wrappers.R        fit_girt_retrospective / fit_girt_realtime        [girt-native]
  solve.R           gi_solve (KKT IRLS)                               [engine]
  lambda_select.R   gi_select_lambda_cv (+ _dow)                      [engine]
  extract.R         gi_extract_rt, gi_estimate_dispersion            [engine]
  taper.R           gi_build_tapered_penalty, .gi_solve_taper,        [engine]
                    gi_tune_gamma_fv
  realtime.R        severity / likelihood-start / rt_df helpers       [engine]
  conformal.R       split-conformal real-time CIs                     [engine]
examples/
  demo.R            self-contained: simulate from the model, recover R_t
```

`[engine]` files are the shared penalized-spline Poisson-GLM machinery, copied
from `mechrt` and renamed to the `gi_*` namespace so girt is fully self-contained
(it sources nothing from `mechrt`). `[girt-native]` files define the method.

> This is research code on the way to a package; names (incl. `girt`) may change.
