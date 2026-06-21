# code/analyses/covidestim/_common.R
#
# Shared helpers for reproducing the covidestim comparison figures with girt.
# The mechanistic covidestim pipeline (on the `main` branch) fit MechRt with two
# conventions: a LAG-1 prevalence (pre_I = c(0, head(P_hat, -1))) and the NOMINAL
# mean infectious period as 1/gamma.  In the generation-interval view these are
# captured exactly by the generation interval
#
#     g_cov = c(0, zeta^EI / mean_infectious_nominal)
#
# where zeta^EI = pi_lat * surv_IR is the infectious kernel.  The leading 0 is the
# lag-1; dividing by the nominal mean infectious period (instead of sum(zeta^EI))
# matches MechRt's 1/gamma scaling.  With this g, girt's design X_{t,s} equals
# MechRt's term-for-term, so the fits coincide (verified per figure).  Everything
# else (deconvolution, CV-lambda, DoW, time-varying severity, extraction) is the
# girt engine.

# --- repo root + paths (no dependency on mechrt / _paths.R) --------------------
.find_repo_root <- function() {
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) tryCatch(dirname(normalizePath(f[1])), error = function(e) getwd()) else getwd()
  for (i in 1:9) {
    if (all(c("data", "results", "figures", "code") %in% list.files(d))) return(normalizePath(d))
    p <- dirname(d); if (p == d) break; d <- p
  }
  stop("could not locate repo root (needs data/results/figures/code)")
}
repo_root   <- .find_repo_root()
data_dir    <- file.path(repo_root, "data")
results_dir <- file.path(repo_root, "results")
figures_dir <- file.path(repo_root, "figures")
suppressPackageStartupMessages(library(ConvRt))

ce_data_dir <- file.path(data_dir, "real", "covidestim")
ce_res_dir  <- file.path(results_dir, "real", "covidestim")
ce_fig_dir  <- file.path(figures_dir, "real", "covidestim")
dir.create(ce_fig_dir, recursive = TRUE, showWarnings = FALSE)

# --- delay pmf (same midpoint-binned discrete Gamma as the mechrt scripts) ------
ce_gamma_pmf <- function(mean, sd) gi_discrete_gamma_delay(mean, sd)$pmf

# --- ascertainment: cases.fitted / (infections convolved with E->case delay) ---
# Clamped to [0.03, 0.8].  Computed from covidestim's raw cases.fitted +
# infections (depends only on pi_EC, so it is shared across GI variants).
ce_ascertainment <- function(cases_fitted, infections, pi_EC) {
  n <- length(cases_fitted); d_max <- length(pi_EC)
  eff <- vapply(seq_len(n), function(t) {
    k <- seq_len(d_max); s <- t - k; ok <- s >= 1
    if (any(ok)) sum(infections[s[ok]] * pi_EC[k[ok]]) else 0 }, numeric(1))
  pmin(pmax(cases_fitted / pmax(eff, 1e-6), 0.03), 0.8)
}

# Old hand-rolled helper for the covidestim-convention GI (lag-1 + nominal-mip).
# Now produced by gi_from_compartmental(pi_lat, pi_IR, mean_infectious_nominal = mip,
# lag_one = TRUE)$g .  Kept as a comment so the caller's intent is documented.
# ce_gi_lag1_nominal <- function(pi_lat, pi_IR, mean_infectious_nominal) {
#   surv_IR <- c(1, 1 - cumsum(pi_IR))
#   kernel  <- .gi_infectious_kernel(pi_lat, surv_IR)
#   c(0, kernel / mean_infectious_nominal)
# }

# --- one covidestim girt fit (girt engine; takes a GI directly) ---------------
# Prior signature took (pi_lat, pi_IR, mean_infectious, N) and reconstructed the
# GI inside.  Now callers build g themselves -- typically via
#   g <- gi_from_compartmental(pi_lat, pi_IR, mean_infectious_nominal = mip,
#                              lag_one = TRUE)$g
# and pass it in.  Returns design + fit + X_hat + lambda; callers run
# gi_extract_rt as needed.
fit_covidestim_girt <- function(cases, ascertain, dates, pi_EC, g,
                                beta_start_date, likelihood_start,
                                knot_step = 5L, tail_type = "none",
                                lambda_grid = 10^seq(-3, 6, length.out = 45),
                                nfold = 5L, cv_rule = "1se", use_dow = FALSE) {
  n <- length(cases); end_t <- n - 1L; day <- 0:end_t
  sev_ref  <- median(ascertain, na.rm = TRUE)
  y_scaled <- pmax(0, round(cases * sev_ref / ascertain))

  X_hat <- gi_deconvolve_exposures(y_scaled, pi_EY = pi_EC, severity_rate = sev_ref,
                                   burn_in = 30L)$X_hat

  beta_start_d <- day[dates == as.Date(beta_start_date)]
  lik_start_d  <- day[dates == as.Date(likelihood_start)]

  d <- suppressWarnings(build_gi_design(
    E_star = X_hat, g = g, pi_EY = pi_EC, max_time = end_t, inf_inc = cases,
    y_min_count = 0, y_min_start = lik_start_d, y_min_end = beta_start_d,
    knot_step = knot_step, dates = dates,
    dow_dates = if (use_dow) dates else NULL))    # g may not sum to 1 (covidestim nominal-mip convention)
  d <- gi_enforce_likelihood_start(d, day[which(d$valid_mask)], lik_start_d)
  d <- gi_apply_severity_to_design(d, ascertain)

  sel <- gi_select_lambda_cv(d, lambda_grid, tail = tail_type, nfold = nfold,
                             cv_select_rule = cv_rule, error_measure = "deviance")
  lam <- if (cv_rule == "1se" && !is.null(sel$best_lam_1se) && is.finite(sel$best_lam_1se))
    sel$best_lam_1se else if (!is.null(sel$best_lam_min) && is.finite(sel$best_lam_min))
    sel$best_lam_min else sel$best_lam
  fit <- gi_solve(d, lam, tail = tail_type)

  list(design = d, fit = fit, X_hat = X_hat, lambda = lam,
       lam_min = sel$best_lam_min, lam_1se = sel$best_lam_1se, sel = sel,
       beta_start_d = beta_start_d, lik_start_d = lik_start_d)
}
