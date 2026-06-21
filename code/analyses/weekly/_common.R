# code/analyses/weekly/_common.R
#
# Shared setup for the girt WEEKLY proof-of-concept -- the girt port of the
# MechRt weekly_demo (code/mechrt/weekly_demo on `main`).  Weekly-aggregated
# girt: the daily renewal design is collapsed to weekly Poisson observations
# (gi_aggregate_design_weekly) and the same penalized smoothing-spline R_t is
# fit through the weekly likelihood.  No day-of-week (a 7-day sum collapses it).
#
# Nothing here is mechanistic: the only epidemiological inputs are a single
# generation interval `g` (a discrete-Gamma pmf) and a reporting delay pi_EY.
# There are NO E->I / I->R compartments or distributions anywhere.
#
# Two datasets, mirroring the MechRt demo:
#   sim  : flu "wiggly" simulation (raw Poisson obs, known truth rt_raw),
#          IDENTITY deconvolution, GI = discrete Gamma(3.2, 1.6), cv rule "min".
#   real : US flu hosps season 2 (2023/24), LOG deconvolution, GI = discrete
#          Gamma(3.2, 1.6) (the girt-flu convention), cv rule "1se".
#
# Both build the daily design with dow_dates = NULL.

.find_repo_root <- function() {
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) tryCatch(dirname(normalizePath(f[1])), error = function(e) getwd()) else getwd()
  for (i in 1:9) {
    if (all(c("data", "results", "figures", "code") %in% list.files(d))) return(normalizePath(d))
    p <- dirname(d); if (p == d) break; d <- p
  }
  stop("weekly/_common.R: repo root not found")
}
repo_root   <- .find_repo_root()
data_dir    <- file.path(repo_root, "data")
results_dir <- file.path(repo_root, "results")
figures_dir <- file.path(repo_root, "figures")
suppressPackageStartupMessages(library(ConvRt))
suppressPackageStartupMessages({ library(dplyr); library(tibble) })

# Output trees for the weekly PoC.
wk_results_dir <- file.path(results_dir, "extra_methods", "weekly")
wk_figures_dir <- file.path(figures_dir, "extra_methods", "weekly")

`%||%` <- function(a, b) if (is.null(a)) b else a

# --- Shared epidemiological inputs -------------------------------------------
wk_cfg <- list(
  N             = 332000000L,
  severity_rate = 0.015,
  mean_EH       = 5.7, sd_EH = 2.3,     # exposure -> hospitalization reporting delay
  knot_step     = 5L,
  level         = 0.95,
  # generation intervals (single discrete-Gamma pmf each; no compartments)
  sim_gi_mean  = 3.2,  sim_gi_sd  = 1.6,    # sim GI
  real_gi_mean = 3.2,  real_gi_sd = 1.6     # real-flu GI (matches EpiNow2)
)

pi_EH_obj <- function() gi_discrete_gamma_delay(wk_cfg$mean_EH, wk_cfg$sd_EH)

# Generation intervals: plain discrete-Gamma pmfs.  (The sim GI's moments match
# the generation interval the simulation transmits under; it is supplied directly
# as a GI, not reconstructed from any latent/infectious-period distributions.)
sim_gi  <- function() gi_discrete_gamma_delay(wk_cfg$sim_gi_mean,  wk_cfg$sim_gi_sd)$pmf
real_gi <- function() gi_discrete_gamma_delay(wk_cfg$real_gi_mean, wk_cfg$real_gi_sd)$pmf

# --- Build a DAILY girt design (no DoW) for one window -----------------------
# Deconvolve obs -> exposures, renewal design, likelihood-start + severity.
# Returns list(design, dates, day, X_hat).  `lik_date` defaults to beta_date.
wk_build_daily_design <- function(obs_inc, dates, g, beta_date, lik_date = NULL,
                                  severity = wk_cfg$severity_rate,
                                  knot_step = wk_cfg$knot_step, link = "identity") {
  dates   <- as.Date(dates)
  if (is.null(lik_date)) lik_date <- beta_date
  pi_EH   <- pi_EH_obj()$pmf
  n       <- length(obs_inc); end_t <- n - 1L; day <- 0:end_t
  X_hat   <- gi_deconvolve_exposures(obs_inc, pi_EY = pi_EH, severity_rate = severity,
                                     link = link, burn_in = 30L)$X_hat
  d <- suppressWarnings(build_gi_design(
    E_star = X_hat, g = g, pi_EY = pi_EH, max_time = end_t, inf_inc = obs_inc,
    y_min_count = 0, y_min_start = 0, y_min_end = 0, knot_step = knot_step,
    first_rt_index = beta_date, dates = dates, dow_dates = NULL))
  d <- gi_enforce_likelihood_start(d, day[which(d$valid_mask)],
                                   day[dates == as.Date(lik_date)])
  d <- gi_apply_severity_to_design(d, rep(severity, n))
  list(design = d, dates = dates, day = day, X_hat = X_hat)
}

load_flu_hosps <- function() {
  csv <- file.path(data_dir, "real", "flu", "hhs_flu_hosps_us_1d_clean.csv")
  d <- read.csv(csv); d$date <- as.Date(d$date)
  d[order(d$date), c("date", "hosps_1d")]
}

# Calendar anchors.
wk_sim <- list(
  sim_rds   = file.path(data_dir, "sim", "flu", "wiggly", "tuning_sim_results.rds"),
  obs_col   = "obs_cases_raw", truth_col = "rt_raw",
  beta_date = as.Date("2022-07-01"), lik_date = as.Date("2022-07-22"),
  rt_start  = as.Date("2022-07-01"), rt_end = as.Date("2023-04-30"),
  lam_grid  = 10^seq(-2, 7, length.out = 36)
)
wk_real_s2 <- list(
  beta_date  = as.Date("2023-07-01"), eos = as.Date("2024-04-26"),
  rt_start   = as.Date("2023-07-01"),
  lam_grid   = 10^seq(0, 7, length.out = 30),
  gamma_grid = 10^seq(-2, 5, length.out = 18),
  cutoffs    = as.Date("2024-01-28"),   # the cutoff the combined figure uses
  n_fv = 4L, cv_folds = 5L
)
