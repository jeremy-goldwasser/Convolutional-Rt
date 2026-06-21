# girt :: weekly.R
#
# Weekly-aggregated girt: collapse the daily renewal design to weekly Poisson
# observations (Y_w ~ Poisson(sum-of-7-days mu)), then fit the same penalized
# smoothing-spline R_t.  This is the girt port of the MechRt weekly_demo
# (code/mechrt/weekly_demo on the `main` branch): aggregate_design_weekly,
# fit_mechrt_weekly_realtime and its walk-forward gamma tuner .weekly_fv_gamma,
# re-expressed against the gi_* engine.
#
# The spline basis B (hence R_t = B theta) lives on the DAILY day axis and is
# untouched by aggregation -- weekly aggregation only sums the likelihood rows
# (Z, Y, off).  So R_t is still extracted per day via gi_extract_rt /
# gi_rt_df_from_fit; only the observation likelihood becomes weekly.
#
# Day-of-week is unidentifiable after a 7-day sum and is forced off (build the
# daily design with dow_dates = NULL).
#
# Public:
#   gi_aggregate_design_weekly(design, week_size, dates_valid, drop_partial)
#   fit_girt_weekly_retrospective(design, lam_grid, ...)   -> retro weekly fit
#   fit_girt_weekly_realtime(design, lam_grid, gamma_grid, pi_E, ...) -> taper RT
# Lower-level: .gi_weekly_fv_gamma

suppressPackageStartupMessages({ library(tibble); library(dplyr) })

# ------------------------------------------------------------------------------
# Row-aggregate a DAILY gi design to weekly Poisson observations.
# Mirrors mechrt::aggregate_design_weekly.  Sums Z_full / Y_valid / off_valid
# over consecutive `week_size`-day buckets of the valid rows; collapses DoW;
# sets n_obs = n_weeks so the gi_solve / gi_extract_rt lambda normalization
# (lam_eff = lam * n_obs / n_penalty) rebalances to the weekly loss.
# ------------------------------------------------------------------------------
gi_aggregate_design_weekly <- function(
    design,
    week_size    = 7L,
    dates_valid  = NULL,    # Date per valid (daily) row; attaches a rep date/week
    drop_partial = TRUE     # drop weeks with fewer than week_size valid days
) {
  week_size <- as.integer(week_size)
  if (!is.finite(week_size) || week_size < 2L)
    stop("gi_aggregate_design_weekly: week_size must be an integer >= 2.")
  if (is.null(design$valid_mask))
    stop("gi_aggregate_design_weekly: design has no $valid_mask.")

  valid_days <- which(design$valid_mask) - 1L   # 0-based day indices
  n_valid    <- length(valid_days)
  if (n_valid != nrow(design$Z_full))
    stop(sprintf("gi_aggregate_design_weekly: nrow(Z_full)=%d but |valid_mask|=%d.",
                 nrow(design$Z_full), n_valid))
  if (n_valid < week_size)
    stop("gi_aggregate_design_weekly: fewer valid days than one week.")

  bucket        <- (valid_days - valid_days[1L]) %/% week_size
  bucket_counts <- tabulate(bucket + 1L)

  keep_buckets <- if (drop_partial) which(bucket_counts == week_size) - 1L
                  else sort(unique(bucket))
  if (length(keep_buckets) == 0L)
    stop(sprintf("gi_aggregate_design_weekly: no complete %d-day weeks in range.",
                 week_size))

  keep_rows <- which(bucket %in% keep_buckets)
  n_weeks   <- length(keep_buckets)

  # Dense aggregator M (n_weeks x n_valid); small enough for one season.
  M <- matrix(0.0, n_weeks, n_valid)
  for (i in seq_along(keep_buckets))
    M[i, which(bucket == keep_buckets[i])] <- 1.0

  d_w <- design
  d_w$Z_full    <- M %*% design$Z_full
  d_w$Y_valid   <- as.numeric(M %*% design$Y_valid)
  d_w$off_valid <- as.numeric(M %*% design$off_valid)
  d_w$use_dow   <- FALSE
  d_w$dow_valid <- NULL
  d_w$n_obs     <- n_weeks   # drives lambda normalization in gi_solve/gi_extract_rt

  new_mask <- rep(FALSE, length(design$valid_mask))
  new_mask[valid_days[keep_rows] + 1L] <- TRUE
  d_w$valid_mask <- new_mask

  d_w$weekly_aggregator <- M
  d_w$weekly_bucket_id  <- bucket[keep_rows]
  d_w$weekly_kept_days  <- valid_days[keep_rows]   # 0-based day indices
  d_w$weekly_size       <- week_size

  if (!is.null(dates_valid) && length(dates_valid) == n_valid) {
    rep_idx <- vapply(keep_buckets,
                      function(b) max(which(bucket == b)), integer(1L))
    d_w$weekly_dates <- as.Date(dates_valid[rep_idx])
  } else {
    d_w$weekly_dates <- NULL
  }

  d_w
}

# ------------------------------------------------------------------------------
# Retrospective weekly fit: aggregate -> CV lambda -> solve -> Rt.
# Returns rt_df (chosen rule), plus rt_at_lam_min / rt_at_lam_1se.  `dates` is
# the FULL daily Date vector (length max_time+1) so Rt days map back to dates.
# ------------------------------------------------------------------------------
fit_girt_weekly_retrospective <- function(
    design,                       # DAILY design (no DoW), severity applied
    lam_grid,
    dates,                        # full daily Date vector, length max_time+1
    week_size      = 7L,
    dates_valid    = NULL,
    cv_select_rule = c("min", "1se"),
    error_measure  = c("deviance", "mse", "mae"),
    nfold          = 5L,
    tail           = "linear",
    level          = 0.95,
    overdispersion = TRUE,
    verbose        = FALSE
) {
  cv_select_rule <- match.arg(cv_select_rule)
  error_measure  <- match.arg(error_measure)
  if (isTRUE(design$use_dow))
    stop("fit_girt_weekly_retrospective: build the daily design with ",
         "dow_dates = NULL -- DoW is unidentifiable after weekly aggregation.")

  d_w <- gi_aggregate_design_weekly(design, week_size = week_size,
                                    dates_valid = dates_valid, drop_partial = TRUE)

  sel <- gi_select_lambda_cv(d_w, lam_grid, tail = tail,
                             nfold = min(nfold, nrow(d_w$Z_full) - 1L),
                             cv_select_rule = cv_select_rule,
                             error_measure = error_measure, verbose = verbose)
  lam_min <- sel$best_lam_min
  lam_1se <- sel$best_lam_1se
  lam     <- if (cv_select_rule == "1se") lam_1se else lam_min

  rt_df_at <- function(lx) {
    fit <- gi_solve(d_w, lx, tail = tail)
    gi_rt_df_from_fit(fit, d_w, dates, lx, level = level,
                      overdispersion = overdispersion)
  }
  rt_min <- rt_df_at(lam_min)
  rt_1se <- rt_df_at(lam_1se)
  rt_df  <- if (cv_select_rule == "1se") rt_1se else rt_min

  list(rt_df = rt_df, rt_at_lam_min = rt_min, rt_at_lam_1se = rt_1se,
       lam = lam, lam_min = lam_min, lam_1se = lam_1se,
       sel_lam = sel, design_weekly = d_w, design_daily = design)
}

# ------------------------------------------------------------------------------
# Real-time weekly fit (right-censored vintage).  Port of
# mechrt::fit_mechrt_weekly_realtime: keep the partial right-edge week, tune
# lambda by CV, then a CDF-tapered tail penalty whose strength gamma is tuned by
# walk-forward validation over WEEKS (.gi_weekly_fv_gamma).
# ------------------------------------------------------------------------------
fit_girt_weekly_realtime <- function(
    design,                       # DAILY design, truncated through cutoff, no DoW
    dates,                        # full daily Date vector, length max_time+1
    lam            = NULL,
    lam_grid       = NULL,
    use_taper      = TRUE,
    pi_E           = NULL,        # exposure -> observed-outcome PMF (taper kernel)
    gamma          = NULL,
    gamma_grid     = NULL,
    week_size      = 7L,
    dates_valid    = NULL,
    n_fv           = 4L,
    cv_select_rule = c("min", "1se"),
    error_measure  = c("deviance", "mse", "mae"),
    nfold          = 5L,
    tail           = "linear",
    fv_error       = c("deviance", "mse", "mae"),
    level          = 0.95,
    overdispersion = TRUE,
    verbose        = FALSE
) {
  cv_select_rule <- match.arg(cv_select_rule)
  error_measure  <- match.arg(error_measure)
  fv_error       <- match.arg(fv_error)

  if (isTRUE(design$use_dow))
    stop("fit_girt_weekly_realtime: build the daily design with dow_dates = NULL.")
  if (use_taper && is.null(pi_E))
    stop("fit_girt_weekly_realtime: pi_E required when use_taper = TRUE.")
  if (is.null(lam) && is.null(lam_grid))
    stop("fit_girt_weekly_realtime: provide either lam or lam_grid.")

  # 1. Weekly aggregation -- KEEP the partial right-edge week (real-time signal).
  d_w <- gi_aggregate_design_weekly(design, week_size = week_size,
                                    dates_valid = dates_valid, drop_partial = FALSE)

  # 2. Lambda by CV (unless supplied).
  sel_lam <- NULL
  if (is.null(lam)) {
    sel_lam <- gi_select_lambda_cv(
      d_w, lam_grid, tail = tail,
      nfold = min(nfold, nrow(d_w$Z_full) - 1L),
      cv_select_rule = cv_select_rule, error_measure = error_measure,
      verbose = verbose)
    lam <- if (cv_select_rule == "min" && is.finite(sel_lam$best_lam_min))
             sel_lam$best_lam_min else sel_lam$best_lam
  }

  # 3. Initial no-taper fit.
  fit_main <- gi_solve(d_w, lam, tail = tail)

  # 4. Gamma by walk-forward FV (unless supplied).
  P_taper <- NULL; fv_info <- NULL
  if (isTRUE(use_taper)) {
    P_taper <- gi_build_tapered_penalty(d_w, pi_E, quiet = !verbose)$P_taper
    if (is.null(gamma)) {
      if (is.null(gamma_grid))
        stop("fit_girt_weekly_realtime: with use_taper and gamma = NULL, ",
             "provide gamma_grid.")
      fv_info <- .gi_weekly_fv_gamma(
        design_w = d_w, lam = lam, P_taper = P_taper, gamma_grid = gamma_grid,
        tail = tail, n_fv = n_fv, theta_init = fit_main$theta,
        error_measure = fv_error, verbose = verbose)
      gamma <- fv_info$best_gamma
    }
  } else {
    gamma <- 0
  }

  # 5. Final fit.
  fit <- if (isTRUE(use_taper) && gamma > 0)
    .gi_solve_taper(d_w, lam, lam_taper = gamma, P_taper = P_taper,
                    tail = tail, theta_init = fit_main$theta)
  else fit_main

  # 6. Rt extraction (daily axis).
  rt_df <- gi_rt_df_from_fit(fit, d_w, dates, lam,
                             lam_taper = if (use_taper) gamma else 0,
                             P_taper = P_taper, level = level,
                             overdispersion = overdispersion)

  list(fit = fit, fit_main = fit_main, design_daily = design, design_weekly = d_w,
       lam = lam, gamma = gamma, sel_lam = sel_lam, fv_info = fv_info,
       P_taper = P_taper, rt_df = rt_df,
       config = list(use_taper = isTRUE(use_taper), tail = tail,
                     n_fv = n_fv, week_size = week_size))
}

# ------------------------------------------------------------------------------
# Walk-forward gamma tuning on the weekly design (port of .weekly_fv_gamma).
# Fold k holds out the k-th most recent week, trains on the preceding weekly
# rows, and scores the held-out weekly count under Z[held] %*% theta + off.
# Z_full was built on the full daily horizon, so the held-out week's Z row
# already exists -- no spline rebuild needed.
# ------------------------------------------------------------------------------
.gi_weekly_fv_gamma <- function(
    design_w, lam, P_taper, gamma_grid,
    tail = "linear", n_fv = 4L, theta_init = NULL,
    error_measure = c("deviance", "mse", "mae"), verbose = FALSE
) {
  error_measure <- match.arg(error_measure)
  err_fun <- switch(error_measure,
    mse      = function(y, m) (y - m)^2,
    mae      = function(y, m) abs(y - m),
    deviance = function(y, m) {
      term <- ifelse(y > 0, y * log(y / pmax(m, 1e-12)) - (y - m), m); 2 * term
    })

  n_w <- nrow(design_w$Z_full)
  if (n_fv >= n_w)
    stop(".gi_weekly_fv_gamma: n_fv (", n_fv, ") must be < n_weeks (", n_w, ").")

  fold_err <- matrix(NA_real_, nrow = n_fv, ncol = length(gamma_grid))
  for (i in seq_along(gamma_grid)) {
    g    <- gamma_grid[i]; warm <- theta_init
    for (k in seq_len(n_fv)) {
      held <- n_w - k + 1L
      if (held <= 1L) next
      train <- seq_len(held - 1L)

      d_tr <- design_w
      d_tr$Y_valid   <- design_w$Y_valid[train]
      d_tr$Z_full    <- design_w$Z_full[train, , drop = FALSE]
      d_tr$off_valid <- design_w$off_valid[train]
      d_tr$n_obs     <- length(train)

      fit_k <- tryCatch(
        if (g > 0)
          .gi_solve_taper(d_tr, lam, lam_taper = g, P_taper = P_taper,
                          tail = tail, theta_init = warm)
        else gi_solve(d_tr, lam, tail = tail, theta_init = warm),
        error = function(e) NULL)
      if (is.null(fit_k) || !all(is.finite(fit_k$theta))) next
      warm <- fit_k$theta

      Z_h  <- design_w$Z_full[held, , drop = FALSE]
      mu_h <- pmax(as.numeric(Z_h %*% fit_k$theta) + design_w$off_valid[held], 1e-10)
      fold_err[k, i] <- err_fun(design_w$Y_valid[held], mu_h)
    }
    if (verbose)
      cat(sprintf("  [FV gamma %d/%d] gamma=%.3e  mean_err=%.4f\n",
                  i, length(gamma_grid), g, mean(fold_err[, i], na.rm = TRUE)))
  }

  fv_mean <- apply(fold_err, 2L, function(x) {
    x <- x[is.finite(x)]; if (length(x) == 0L) Inf else mean(x) })
  fv_se <- apply(fold_err, 2L, function(x) {
    x <- x[is.finite(x)]
    if (length(x) <= 1L) NA_real_ else sd(x) / sqrt(length(x)) })
  best_i <- which.min(fv_mean)
  list(gamma_grid = gamma_grid, fold_err = fold_err, fv_mean = fv_mean,
       fv_se = fv_se, best_i = best_i, best_gamma = gamma_grid[best_i], n_fv = n_fv)
}
