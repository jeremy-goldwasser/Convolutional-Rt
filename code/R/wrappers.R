# girt :: wrappers.R
#
# High-level one-call entry points for the generation-interval Rt method.
#
#   fit_girt_retrospective(...)  end-of-season fit: deconvolve -> renewal design
#                                -> CV-tuned smoothing spline -> Rt + pointwise CI.
#   fit_girt_realtime(...)       right-censored fit: same, plus a tapered tail
#                                penalty whose strength gamma is tuned by forward
#                                validation; optional split-conformal CIs need
#                                daily vintages (see R/conformal.R).
#
# Inputs are a generation interval g, a reporting delay (mean/sd of the
# exposure->observation Gamma), the observed counts, and the severity rho.  No
# compartmental quantities.  (To match a specific SEIR model, build g with
# gi_from_compartmental().)
#
# Source girt.R before use.

suppressPackageStartupMessages({ library(tibble); library(dplyr) })

# ------------------------------------------------------------------------------
# Retrospective
# ------------------------------------------------------------------------------
fit_girt_retrospective <- function(
    obs_inc, dates, g,
    mean_EY, sd_EY,                 # exposure -> observed-outcome (reporting) delay
    severity        = 1,            # rho (scalar or length-n vector)
    first_rt_date   = NULL,         # burn-in anchor (Date or 0-based index); default = series start
    likelihood_start_date = NULL,   # first day entering the likelihood
    knot_step       = 5L,
    lam_grid        = 10^seq(-2, 8, length.out = 30),
    cv_select_rule  = c("min", "1se"),
    error_measure   = c("deviance", "mse", "mae"),
    nfold           = 5L,
    dow_dates       = NULL,
    deconv_link     = c("identity", "log"),  # deconvolution link (default identity)
    level           = 0.95,
    overdispersion  = TRUE,
    X_hat           = NULL          # supply to skip the internal deconvolution
) {
  cv_select_rule <- match.arg(cv_select_rule); error_measure <- match.arg(error_measure)
  deconv_link <- match.arg(deconv_link)
  dates <- as.Date(dates); n <- length(obs_inc); max_time <- n - 1L
  pi_EY <- gi_discrete_gamma_delay(mean_EY, sd_EY)$pmf
  sev   <- gi_resolve_severity(n, severity)

  if (is.null(X_hat))
    X_hat <- gi_deconvolve_exposures(obs_inc, pi_EY, severity_rate = mean(sev), link = deconv_link)$X_hat

  if (is.null(first_rt_date)) first_rt_date <- 0L
  d <- build_gi_design(E_star = X_hat, g = g, pi_EY = pi_EY, max_time = max_time,
                       inf_inc = obs_inc, knot_step = knot_step,
                       first_rt_index = first_rt_date, dates = dates, dow_dates = dow_dates)
  day <- 0:max_time
  if (!is.null(likelihood_start_date)) {
    ls_day <- day[dates == as.Date(likelihood_start_date)]
    d <- gi_enforce_likelihood_start(d, day[which(d$valid_mask)], ls_day)
  }
  d <- gi_apply_severity_to_design(d, sev)

  sel <- gi_select_lambda_cv(d, lam_grid, tail = "linear", nfold = nfold,
                             cv_select_rule = cv_select_rule, error_measure = error_measure)
  lam <- if (cv_select_rule == "min" && is.finite(sel$best_lam_min)) sel$best_lam_min else sel$best_lam
  fit <- gi_solve(d, lam, tail = "linear")
  rt  <- gi_rt_df_from_fit(fit, d, dates, lam, level = level, overdispersion = overdispersion)
  phi <- tryCatch(gi_estimate_dispersion(d, fit, lam), error = function(e) NA_real_)

  list(rt_df = rt, lam = lam, fit = fit, design = d, X_hat = X_hat,
       g = g, pi_EY = pi_EY, sel_lam = sel, phi_hat = phi)
}

# ------------------------------------------------------------------------------
# Real-time (right-censored) with forward-validated tapered tail penalty
# ------------------------------------------------------------------------------
fit_girt_realtime <- function(
    obs_inc, dates, g,
    mean_EY, sd_EY,
    severity        = 1,
    first_rt_date   = NULL,
    likelihood_start_date = NULL,
    knot_step       = 5L,
    lam_grid        = 10^seq(-2, 8, length.out = 30),
    cv_select_rule  = c("min", "1se"),
    error_measure   = c("deviance", "mse", "mae"),
    nfold           = 5L,
    use_taper       = TRUE,
    gamma_grid      = 10^seq(-4, 5, length.out = 18),
    n_fv            = 7L,
    dow_dates       = NULL,
    deconv_link     = c("identity", "log"),  # deconvolution link (default identity)
    level           = 0.95,
    overdispersion  = TRUE,
    X_hat           = NULL
) {
  cv_select_rule <- match.arg(cv_select_rule); error_measure <- match.arg(error_measure)
  deconv_link <- match.arg(deconv_link)
  dates <- as.Date(dates); n <- length(obs_inc); max_time <- n - 1L
  ey    <- gi_discrete_gamma_delay(mean_EY, sd_EY); pi_EY <- ey$pmf
  sev   <- gi_resolve_severity(n, severity); rho_scalar <- mean(sev)

  if (is.null(X_hat))
    X_hat <- gi_deconvolve_exposures(obs_inc, pi_EY, severity_rate = rho_scalar, link = deconv_link)$X_hat
  if (is.null(first_rt_date)) first_rt_date <- 0L
  day <- 0:max_time
  ls_day <- if (!is.null(likelihood_start_date)) day[dates == as.Date(likelihood_start_date)] else NULL

  build_one <- function(X_v, y_v, dts_v, dow_v) {
    mt <- length(X_v) - 1L
    dd <- build_gi_design(E_star = X_v, g = g, pi_EY = pi_EY, max_time = mt,
                          inf_inc = y_v, knot_step = knot_step,
                          first_rt_index = first_rt_date, dates = dts_v, dow_dates = dow_v)
    vday <- (0:mt)[which(dd$valid_mask)]
    if (!is.null(likelihood_start_date)) {
      lsd <- (0:mt)[dts_v == as.Date(likelihood_start_date)]
      dd <- gi_enforce_likelihood_start(dd, vday, lsd)
    }
    gi_apply_severity_to_design(dd, gi_resolve_severity(length(X_v), severity))
  }

  d <- build_one(X_hat, obs_inc, dates, dow_dates)

  sel <- gi_select_lambda_cv(d, lam_grid, tail = "linear", nfold = nfold,
                             cv_select_rule = cv_select_rule, error_measure = error_measure)
  lam <- if (cv_select_rule == "min" && is.finite(sel$best_lam_min)) sel$best_lam_min else sel$best_lam
  fit_main <- gi_solve(d, lam, tail = "linear")

  gamma <- 0; fit <- fit_main; P_taper <- NULL
  if (isTRUE(use_taper)) {
    fold_builder <- function(s) {
      idx <- seq_len(s + 1L)
      dd  <- tryCatch(build_one(X_hat[idx], obs_inc[idx], dates[idx],
                                if (is.null(dow_dates)) NULL else dow_dates[idx]),
                      error = function(e) NULL)
      if (is.null(dd)) return(NULL)
      ts_s <- tibble(S = 0, I = dd$gi_force, pre_S = 0, pre_I = dd$gi_force)
      list(design = dd, time_series = ts_s, pi_E = pi_EY)
    }
    ts_full <- tibble(S = 0, I = d$gi_force, pre_S = 0, pre_I = d$gi_force)
    sg <- gi_tune_gamma_fv(
      time_series = ts_full, inf_inc = obs_inc,
      report_shape = ey$shape, report_rate = ey$rate,   # reporting-delay Gamma (for per-fold taper recon)
      t_current = max_time, lam_tuned = lam, retro_theta = fit_main$theta,
      gamma_grid = gamma_grid, n_fv = n_fv, knot_spacing = knot_step,
      y_min_count = 0, tail = "linear", verbose = FALSE,
      spline_degree = 3L,
      basis = "natural", knot_step = knot_step, severity_rate = sev,
      dow_effects = fit_main$dow_effects, dow_dates = dow_dates,
      fold_builder = fold_builder)
    gamma <- if (!is.null(sg$best_gamma_min) && is.finite(sg$best_gamma_min)) sg$best_gamma_min else sg$best_gamma
    if (gamma > 0) {
      P_taper <- gi_build_tapered_penalty(d, pi_EY, quiet = TRUE)$P_taper
      fit <- .gi_solve_taper(d, lam, lam_taper = gamma, P_taper = P_taper,
                             tail = "linear", theta_init = fit_main$theta)
    }
  }

  rt <- gi_rt_df_from_fit(fit, d, dates, lam, lam_taper = gamma, P_taper = P_taper,
                          level = level, overdispersion = overdispersion)
  phi <- tryCatch(gi_estimate_dispersion(d, fit_main, lam), error = function(e) NA_real_)

  list(rt_df = rt, lam = lam, gamma = gamma, fit = fit, fit_main = fit_main,
       design = d, X_hat = X_hat, g = g, pi_EY = pi_EY, sel_lam = sel, phi_hat = phi)
}
