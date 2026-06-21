# code/analyses/flu/_realtime.R
#
# girt real-time fit for one weekly/daily vintage, reproducing MechRt's
# fit_one_combo canonical "natural / linear / taper / cv_1se" path with LOG
# deconvolution and g = zeta^EI/2.75:
#   deconvolve(log) -> renewal design (DoW) -> CV-1se-DoW lambda -> solve (DoW)
#   -> FV-min gamma (taper) -> DoW-prescale -> tapered solve -> extract rt_df.
# Source _common.R first.

# Full per-fold gamma FV reusing the full-horizon X_hat (no re-deconv per fold),
# DoW-prescaled by the supplied dow_effects (matches fv_tune_gamma_realtime).
.flu_fv_tune_gamma <- function(hosps_fit, X_hat_full, g, pi_EH, eh_shape, eh_rate,
                               lam, dow_effects, beta, lik, severity, knot_step,
                               end_t, gamma_grid, n_fv) {
  fold_builder <- function(s) {
    hs <- hosps_fit[hosps_fit$day <= s, ]
    if (!(beta %in% hs$date) || !(lik %in% hs$date)) return(NULL)
    X_s <- X_hat_full[seq_len(nrow(hs))]
    d_s <- tryCatch(build_gi_design(E_star = X_s, g = g, pi_EY = pi_EH, max_time = nrow(hs) - 1L,
            inf_inc = hs$hosps_1d, y_min_count = 0, y_min_start = 0, y_min_end = 0,
            knot_step = knot_step, first_rt_index = beta, dates = hs$date, dow_dates = hs$date),
          error = function(e) NULL)
    if (is.null(d_s)) return(NULL)
    vd <- (0:(nrow(hs) - 1L))[which(d_s$valid_mask)]
    d_s <- tryCatch(gi_enforce_likelihood_start(d_s, vd, (0:(nrow(hs) - 1L))[hs$date == lik]),
                    error = function(e) NULL)
    if (is.null(d_s)) return(NULL)
    d_s <- gi_apply_severity_to_design(d_s, rep(severity, nrow(hs)))
    if (!is.null(dow_effects) && isTRUE(d_s$use_dow)) {          # pre-scale DoW
      dt <- dow_effects[d_s$dow_valid]
      d_s$Z_full <- d_s$Z_full * dt; d_s$off_valid <- d_s$off_valid * dt
      d_s$use_dow <- FALSE; d_s$dow_valid <- NULL
    }
    list(design = d_s, time_series = tibble::tibble(S = 0, I = d_s$gi_force, pre_S = 0, pre_I = d_s$gi_force),
         pi_E = pi_EH)
  }
  gi_tune_gamma_fv(
    time_series = tibble::tibble(S = 0, I = gi_renewal_force(X_hat_full, g), pre_S = 0,
                                 pre_I = gi_renewal_force(X_hat_full, g)),
    inf_inc = hosps_fit$hosps_1d, report_shape = eh_shape, report_rate = eh_rate,
    t_current = end_t, lam_tuned = lam, retro_theta = NULL, gamma_grid = gamma_grid,
    n_fv = n_fv, knot_spacing = knot_step, y_min_count = 0, tail = "linear",
    verbose = FALSE, spline_degree = 3L, basis = "natural", knot_step = knot_step,
    severity_rate = severity, dow_effects = dow_effects, dow_dates = hosps_fit$date,
    fold_builder = fold_builder)
}

# Regression-spline GI design with constant/linear tail KKT constraints (for the
# "regression_constant" real-time comparison variant).  Same GI force X_{t,s} =
# Lambda_s * pi_EY[t-s] as build_gi_design, but a clamped cubic B-spline basis +
# derivative tail constraints, mirroring MechRt's regression branch.  gi_solve's
# KKT path enforces the constraints.
build_gi_design_rs <- function(E_star, g, pi_EY, max_time, inf_inc, fixed_knot_step,
                               dates = NULL, first_rt_index = NULL, dow_dates = NULL,
                               y_min_count = 0, spline_degree = 3L) {
  Lambda <- gi_renewal_force(E_star, g); d_max <- length(pi_EY); ord <- spline_degree + 1L
  if (!is.null(first_rt_index) && inherits(first_rt_index, "Date"))
    first_rt_index <- match(first_rt_index, as.Date(dates)) - 1L
  if (is.null(first_rt_index)) first_rt_index <- 0L
  first_count_index <- first_rt_index + d_max
  X_ts <- matrix(0.0, max_time + 1L, max_time + 1L)
  for (t in seq_len(max_time)) for (s in max(0L, t - d_max):(t - 1L))
    X_ts[t + 1L, s + 1L] <- Lambda[s + 1L] * pi_EY[t - s]
  Y_all <- as.numeric(inf_inc); rs_X <- rowSums(X_ts)
  burn <- seq_along(Y_all) > first_count_index
  valid_mask <- (rs_X > 0) & (Y_all > y_min_count) & burn
  if (!any(valid_mask)) valid_mask <- (rs_X > 0) & (Y_all > 0) & burn
  Xv <- X_ts[valid_mask, , drop = FALSE]; Y_valid <- Y_all[valid_mask]; off_valid <- numeric(sum(valid_mask))
  t_data_end <- max(which((rs_X > 0) & (Y_all > y_min_count))) - 1L
  fr <- c(1L, as.integer(t_data_end))
  sknots <- seq(fixed_knot_step, max_time - fixed_knot_step, by = fixed_knot_step)
  sknots <- sknots[sknots < fr[2]]; if (!length(sknots)) sknots <- fr[1]
  if (fr[2] - max(sknots) > 9L) sknots <- sort(unique(c(sknots, as.numeric(fr[2] - 7L))))
  t_tail <- max(sknots); int_knots <- sknots[sknots > fr[1] & sknots < fr[2]]
  if (!length(int_knots)) int_knots <- (fr[1] + fr[2]) / 2.0
  t_tail <- max(int_knots)
  all_knots <- c(rep(fr[1], ord), int_knots, rep(fr[2], ord))
  t_eval <- 0:max_time
  B <- splines::splineDesign(all_knots, t_eval, ord = ord, derivs = 0L, outer.ok = TRUE)
  n_basis <- ncol(B); np <- n_basis; Z_full <- Xv %*% B
  D2 <- diff(diag(n_basis), differences = 2L); P_pspline <- crossprod(D2); P_ridge <- diag(n_basis)
  t_c <- (t_tail + fr[2]) / 2.0
  d1 <- splines::splineDesign(all_knots, t_c, ord = ord, derivs = 1L, outer.ok = TRUE)[1L, ]
  d2 <- splines::splineDesign(all_knots, t_c, ord = ord, derivs = 2L, outer.ok = TRUE)[1L, ]
  linear_right <- matrix(d2, nrow = 1L)
  if (spline_degree >= 3L) linear_right <- rbind(linear_right,
    matrix(splines::splineDesign(all_knots, t_c, ord = ord, derivs = 3L, outer.ok = TRUE)[1L, ], nrow = 1L))
  empty_A <- matrix(0.0, 0L, np)
  tail_constraints <- list(none = empty_A, linear = linear_right,
                           constant = rbind(linear_right, matrix(d1, nrow = 1L)))
  use_dow <- !is.null(dow_dates) && length(dow_dates) == (max_time + 1L)
  dow_all <- if (use_dow) as.integer(format(as.Date(dow_dates), "%u")) else NULL
  dow_valid <- if (use_dow) dow_all[valid_mask] else NULL
  list(Z_full = Z_full, Y_valid = Y_valid, off_valid = off_valid, B = B, n_basis = n_basis, np = np,
       n_jump = 0L, jump_times_idx = integer(0), D2 = D2, P_pspline = P_pspline, P_ridge = P_ridge,
       tail_constraints = tail_constraints, tail_const_extra_M = NULL, tail_const_days = NA_integer_,
       head_const_extra_M = NULL, head_const_days = 0L, all_knots = all_knots, int_knots = int_knots,
       fr = fr, t_tail = t_tail, kink_vec = NULL, use_kink = FALSE, T_kink = NULL, valid_mask = valid_mask,
       first_rt_index = first_rt_index, first_count_index = first_count_index, max_time = max_time,
       spline_degree = spline_degree, knot_step = fixed_knot_step,
       use_dow = use_dow, dow_valid = dow_valid, dow_all = dow_all, basis = "regression",
       n_obs = sum(valid_mask), n_penalty = max(n_basis - 2L, 1L), gi_force = Lambda)
}

# Regression-constant real-time fit (notaper, CV-1se-DoW, tail=constant).
fit_flu_constant_girt <- function(hosps_all, week_end, beta, lik, g, pEH_pmf,
                                  severity = flu_cfg$severity_rate, fixed_knot_step = 5L,
                                  lam_grid = 10^seq(0, 5, length.out = 30), nfold = flu_cfg$cv_folds,
                                  link = flu_cfg$deconv_link) {
  win <- hosps_all[hosps_all$date >= beta - 30L & hosps_all$date <= week_end, ]
  win$day <- seq_len(nrow(win)) - 1L
  X_hat <- gi_deconvolve_exposures(win$hosps_1d, pi_EY = pEH_pmf, severity_rate = severity, link = link, burn_in = 30L)$X_hat
  d <- build_gi_design_rs(X_hat, g, pEH_pmf, nrow(win) - 1L, win$hosps_1d, fixed_knot_step,
                          dates = win$date, first_rt_index = beta, dow_dates = win$date)
  d <- gi_enforce_likelihood_start(d, win$day[which(d$valid_mask)], win$day[win$date == lik])
  d <- gi_apply_severity_to_design(d, rep(severity, nrow(win)))
  sel <- gi_select_lambda_cv_dow(d, win$date[which(d$valid_mask)], lam_grid, tail = "constant",
                                 nfold = nfold, cv_select_rule = "1se", error_measure = "deviance")
  lam <- if (!is.null(sel$best_lam_1se) && is.finite(sel$best_lam_1se)) sel$best_lam_1se else sel$best_lam
  fit <- gi_solve(d, lam, tail = "constant")
  rt <- gi_extract_rt(fit, d, lam = lam, level = 0.95, overdispersion = TRUE)
  rt$date <- win$date[rt$day + 1L]
  list(meta = list(week_end_date = as.Date(week_end), spline_type = "regression", tail_mode = "constant",
                   use_taper = FALSE, deconv_link = link, lambda = lam),
       rt_df = rt[, c("date","day","Rt_mean","Rt_lo","Rt_hi")], lambda_chosen = lam)
}

# Retrospective end-of-season fit (notaper, CV-1se-DoW) with pointwise +
# simultaneous band, matching run_season().  Returns list(rt_df[+sim band], lam_1se).
fit_flu_retro_girt <- function(hosps_all, eos, beta, lik, pi_EH_obj, g,
                               severity = flu_cfg$severity_rate, knot_step = flu_cfg$knot_step,
                               lam_grid = flu_cfg$ns_lambda_grid, nfold = flu_cfg$cv_folds,
                               level = 0.90, link = flu_cfg$deconv_link) {
  pi_EH <- pi_EH_obj$pmf
  win <- hosps_all[hosps_all$date >= beta - 30L & hosps_all$date <= eos, ]
  win$day <- seq_len(nrow(win)) - 1L; end_t <- nrow(win) - 1L
  X_hat <- gi_deconvolve_exposures(win$hosps_1d, pi_EY = pi_EH, severity_rate = severity,
                                   link = link, burn_in = 30L)$X_hat
  d <- suppressWarnings(build_gi_design(E_star = X_hat, g = g, pi_EY = pi_EH, max_time = end_t,
        inf_inc = win$hosps_1d, y_min_count = 0, y_min_start = 0, y_min_end = 0,
        knot_step = knot_step, first_rt_index = beta, dates = win$date, dow_dates = win$date))
  d <- gi_enforce_likelihood_start(d, win$day[which(d$valid_mask)], win$day[win$date == lik])
  d <- gi_apply_severity_to_design(d, rep(severity, nrow(win)))
  sel <- gi_select_lambda_cv_dow(d, win$date[which(d$valid_mask)], lam_grid, tail = "linear",
                                 nfold = nfold, cv_select_rule = "1se", error_measure = "deviance")
  lam <- if (!is.null(sel$best_lam_1se) && is.finite(sel$best_lam_1se)) sel$best_lam_1se else sel$best_lam
  fit <- gi_solve(d, lam, tail = "linear")
  rt <- gi_extract_rt_simband(fit, d, lam = lam,
                              level = level, overdispersion = TRUE)
  rt$date <- win$date[rt$day + 1L]
  list(rt_df = rt, lam_1se = lam, dow_effects = fit$dow_effects, hosps = win)
}

# One vintage fit -> list(rt_df, lambda_chosen, gamma_chosen, dow_effects, ...).
# `use_taper=FALSE` reproduces the notaper path; lam_retro_1se (if given) yields
# the rt_at_lam_retro_1se curve too.
fit_flu_combo_girt <- function(hosps_all, week_end, beta, lik,
                               pi_EH_obj, g,
                               severity = flu_cfg$severity_rate, knot_step = flu_cfg$knot_step,
                               lam_grid = flu_cfg$ns_lambda_grid, gamma_grid = flu_cfg$gamma_grid,
                               nfold = flu_cfg$cv_folds, n_fv = flu_cfg$n_fv,
                               use_taper = TRUE, lam_retro_1se = NULL, link = flu_cfg$deconv_link) {
  pi_EH <- pi_EH_obj$pmf
  win <- hosps_all[hosps_all$date >= beta - 30L & hosps_all$date <= week_end, ]
  win$day <- seq_len(nrow(win)) - 1L; end_t <- nrow(win) - 1L
  X_hat <- gi_deconvolve_exposures(win$hosps_1d, pi_EY = pi_EH, severity_rate = severity,
                                   link = link, burn_in = 30L)$X_hat
  d <- suppressWarnings(build_gi_design(E_star = X_hat, g = g, pi_EY = pi_EH, max_time = end_t,
        inf_inc = win$hosps_1d, y_min_count = 0, y_min_start = 0, y_min_end = 0,
        knot_step = knot_step, first_rt_index = beta, dates = win$date, dow_dates = win$date))
  d <- gi_enforce_likelihood_start(d, win$day[which(d$valid_mask)], win$day[win$date == lik])
  d <- gi_apply_severity_to_design(d, rep(severity, nrow(win)))
  valid_dates <- win$date[which(d$valid_mask)]

  sel <- gi_select_lambda_cv_dow(d, valid_dates, lam_grid, tail = "linear", nfold = nfold,
                                 cv_select_rule = "1se", error_measure = "deviance")
  lam <- if (!is.null(sel$best_lam_1se) && is.finite(sel$best_lam_1se)) sel$best_lam_1se else sel$best_lam
  fit_main <- gi_solve(d, lam, tail = "linear"); dow <- fit_main$dow_effects

  rt_df_from <- function(fit, des, lamx, gam = 0, Pt = NULL)
    { rt <- gi_extract_rt(fit, des, lam = lamx, lam_taper = gam,
                          P_taper = Pt, level = 0.95, overdispersion = TRUE)
      rt$date <- win$date[rt$day + 1L]; rt[, c("date","day","Rt_mean","Rt_lo","Rt_hi")] }

  gamma_chosen <- 0; fit_final <- fit_main; design_rt <- d; Ptaper <- NULL; gamma_info <- NULL
  if (isTRUE(use_taper)) {
    d_sc <- d
    if (!is.null(dow) && isTRUE(d_sc$use_dow)) { dt <- dow[d_sc$dow_valid]
      d_sc$Z_full <- d_sc$Z_full * dt; d_sc$off_valid <- d_sc$off_valid * dt
      d_sc$use_dow <- FALSE; d_sc$dow_valid <- NULL }
    gamma_info <- .flu_fv_tune_gamma(win, X_hat, g, pi_EH, pi_EH_obj$shape, pi_EH_obj$rate,
                    lam, dow, beta, lik, severity, knot_step, end_t, gamma_grid, n_fv)
    gamma_chosen <- gamma_info$best_gamma_min
    Ptaper <- gi_build_tapered_penalty(d_sc, pi_EH, quiet = TRUE)$P_taper
    fit_final <- ConvRt:::.gi_solve_taper(d_sc, lam, lam_taper = gamma_chosen, P_taper = Ptaper,
                                 tail = "linear", theta_init = fit_main$theta)
    design_rt <- d_sc
  }
  rt_df <- rt_df_from(fit_final, design_rt, lam, gamma_chosen, Ptaper)

  out <- list(meta = list(week_end_date = as.Date(week_end), beta_start_date = beta,
                          spline_type = "natural", tail_mode = "linear", use_taper = use_taper,
                          deconv_link = link, lambda = lam, gamma = gamma_chosen),
              rt_df = rt_df, lambda_chosen = lam, gamma_chosen = gamma_chosen,
              dow_effects = dow, gamma_info = gamma_info)
  if (!is.null(lam_retro_1se) && use_taper)
    out$rt_at_lam_retro_1se <- rt_df_from(
      ConvRt:::.gi_solve_taper(design_rt, lam_retro_1se, lam_taper = gamma_chosen, P_taper = Ptaper,
                      tail = "linear", theta_init = fit_main$theta),
      design_rt, lam_retro_1se, gamma_chosen, Ptaper)
  out
}
