# code/analyses/flu/_common.R
#
# Shared helpers for reproducing the real-flu figures with girt.
#
# Generation interval: discrete Gamma(mean = 3.2, sd = 1.6), matching the GI
# given to EpiNow2 on this dataset.  (girt previously derived the GI from a
# compartmental E->I / I->R convolution; that path has been scrubbed -- the
# discrete Gamma has near-identical moments, is the same parametric family every
# other Rt method uses, and the fits were numerically indistinguishable.)
#
# Delays/severity/DoW/knots match realtime_config.R:
#   EH 5.7/2.3, severity 0.015, knot_step 5, DoW on.

# --- repo root + paths (no mechrt dependency) ---------------------------------
.find_repo_root <- function() {
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) tryCatch(dirname(normalizePath(f[1])), error = function(e) getwd()) else getwd()
  for (i in 1:9) {
    if (all(c("data", "results", "figures", "code") %in% list.files(d))) return(normalizePath(d))
    p <- dirname(d); if (p == d) break; d <- p
  }
  stop("could not locate repo root")
}
repo_root   <- .find_repo_root()
data_dir    <- file.path(repo_root, "data")
results_dir <- file.path(repo_root, "results")
figures_dir <- file.path(repo_root, "figures")
source(file.path(repo_root, "code", "girt.R"))
suppressPackageStartupMessages({ library(dplyr) })

flu_data_dir    <- file.path(data_dir, "real", "flu")
flu_results_dir <- file.path(results_dir, "real", "flu")
flu_figures_dir <- file.path(figures_dir, "real", "flu")

# --- config (matches code/real/flu/realtime_config.R) -------------------------
flu_cfg <- list(
  N = 332000000L, severity_rate = 0.015,
  mean_hosp = 5.7, sd_hosp = 2.3,
  knot_step = 5L, ns_lambda_grid = 10^seq(3, 8, length.out = 30),
  gamma_grid = 10^seq(-2, 5, length.out = 18), cv_folds = 5L, n_fv = 7L,
  deconv_link = "log",   # real flu uses LOG deconvolution (girt default is identity)
  csv = file.path(flu_data_dir, "hhs_flu_hosps_us_1d_clean.csv"))
flu_seasons <- list(
  s1 = list(beta = as.Date("2022-07-01"), eos = as.Date("2023-05-01")),
  s2 = list(beta = as.Date("2023-07-01"), eos = as.Date("2024-04-26")))

ce_gamma_pmf <- function(mean, sd) gi_discrete_gamma_delay(mean, sd)$pmf

# Generation interval: discrete Gamma matching EpiNow2's GI (mean = 3.2, sd = 1.6).
flu_gi <- function() gi_discrete_gamma_delay(3.2, 1.6)$pmf

load_flu_hosps <- function(csv = flu_cfg$csv) {
  d <- read.csv(csv); d$date <- as.Date(d$date)
  d[order(d$date), c("date", "hosps_1d")]
}

# Build the girt design for one (trimmed) flu window: deconvolve -> renewal design
# (DoW on) -> likelihood-start + severity.  Returns design + X_hat + dates + days.
flu_build_design_girt <- function(hosp_inc, dates, beta_date, lik_date,
                                  pi_EH, g, severity, knot_step = 5L,
                                  link = flu_cfg$deconv_link) {
  n <- length(hosp_inc); end_t <- n - 1L; day <- 0:end_t
  X_hat <- gi_deconvolve_exposures(hosp_inc, pi_EY = pi_EH, severity_rate = severity,
                                   link = link, burn_in = 30L)$X_hat
  d <- suppressWarnings(build_gi_design(
    E_star = X_hat, g = g, pi_EY = pi_EH, max_time = end_t, inf_inc = hosp_inc,
    y_min_count = 0, y_min_start = 0, y_min_end = 0, knot_step = knot_step,
    first_rt_index = beta_date, dates = dates, dow_dates = dates))
  d <- gi_enforce_likelihood_start(d, day[which(d$valid_mask)], day[dates == as.Date(lik_date)])
  d <- gi_apply_severity_to_design(d, rep(severity, n))
  list(design = d, X_hat = X_hat, day = day)
}
