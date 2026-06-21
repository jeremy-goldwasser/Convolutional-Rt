# girt :: tf.R
#
# Trend-filter variant of girt: instead of the natural-cubic-spline basis with
# integrated-2nd-derivative^2 penalty, use the IDENTITY basis (theta_s = R_s
# directly) with an L1 penalty on the k-th finite differences of theta:
#
#   minimize_theta   sum_t [mu_t - y_t log mu_t]  +  lam * || D^(k) theta ||_1
#   s.t.             mu_t = Z_t' theta + off_t,   theta >= 0,
#                    Z_{t,s} = Lambda_s * pi_EY[t-s]              (girt kernel)
#
# k = 4 -> cubic trend filter (R_t is piecewise cubic between data-driven knots).
# k = 2 -> piecewise linear.
#
# This is the GI analogue of MechRt's experimental TF variant; same convolution
# kernel as girt's spline path (girt's build_gi_design already builds the right
# Xv = Lambda * pi_EY before spline aggregation), just swap the basis & penalty.
# No CIs (no closed-form Wald story for the L1 selection bias).
#
# Solver: CVXR (ECOS, SCS fallback) with the lambda canonicalized once as a
# Parameter so refits across the lambda path are fast.  Non-negativity on theta
# is enforced -- the cubic-TF null space lets unconstrained theta swing negative
# on this convolution kernel.

if (!requireNamespace("CVXR", quietly = TRUE))
  stop("girt/tf.R: needs the `CVXR` package.")


# ------------------------------------------------------------------------------
# k-th forward-difference matrix (k_diff = 4 = cubic TF).
# ------------------------------------------------------------------------------
.gi_tf_diff_matrix <- function(n, k_diff) {
  D <- diag(n); for (j in seq_len(k_diff)) D <- diff(D); D
}


# ------------------------------------------------------------------------------
# build_gi_design_tf: same args as build_gi_design; returned design has Z_full
# overwritten with the un-aggregated convolution kernel Xv (n_valid x (max_time+1))
# instead of Xv %*% B.  Spline-specific fields are nulled out.
# ------------------------------------------------------------------------------
build_gi_design_tf <- function(E_star, g, pi_EY, max_time, inf_inc, ...) {
  d <- build_gi_design(E_star = E_star, g = g, pi_EY = pi_EY,
                       max_time = max_time, inf_inc = inf_inc, ...)
  # Reconstruct the un-aggregated kernel (matches build_gi_design L88-93).
  d_max <- length(pi_EY); Lambda <- d$gi_force
  X_ts <- matrix(0.0, max_time + 1L, max_time + 1L)
  for (t in seq_len(max_time)) for (s in max(0L, t - d_max):(t - 1L))
    X_ts[t + 1L, s + 1L] <- Lambda[s + 1L] * pi_EY[t - s]
  Xv <- X_ts[d$valid_mask, , drop = FALSE]
  d$Z_full   <- Xv
  d$n_basis  <- ncol(Xv); d$np <- ncol(Xv)
  d$basis    <- "tf_identity"
  d$B        <- NULL; d$P_pspline <- NULL; d$P_ridge <- NULL
  d$tail_constraints <- NULL; d$tail_const_extra_M <- NULL; d$head_const_extra_M <- NULL
  d$all_knots <- NULL; d$int_knots <- NULL; d$D2 <- NULL
  d
}


# ------------------------------------------------------------------------------
# Build the CVXR Problem once (canonicalize), parameterized by lambda so the
# inner CV loop just re-solves the same prob with different lambda values.
# ------------------------------------------------------------------------------
gi_build_tf_problem <- function(design, k_diff = 4L, active_cols = NULL,
                                tail = c("none", "linear", "constant"),
                                n_tail_const = 7L, taper_weights = NULL) {
  tail <- match.arg(tail)
  Z_full  <- design$Z_full
  y       <- as.numeric(design$Y_valid)
  off     <- as.numeric(design$off_valid)
  np_full <- ncol(Z_full)

  if (is.null(active_cols)) {
    col_sums <- colSums(abs(Z_full)); tol_col <- max(col_sums) * 1e-10
    active   <- which(col_sums > tol_col)
    if (length(active) == 0L) stop("gi_build_tf_problem: Z_full has no non-zero columns.")
    keep_cols <- min(active):max(active)
  } else keep_cols <- as.integer(active_cols)
  Z  <- Z_full[, keep_cols, drop = FALSE]; np <- ncol(Z)
  D  <- .gi_tf_diff_matrix(np, k_diff); m <- nrow(D)

  theta <- CVXR::Variable(np, nonneg = TRUE)
  mu    <- Z %*% theta + off
  lam_p <- CVXR::Parameter(value = 1, nonneg = TRUE)
  gam_p <- CVXR::Parameter(value = 0, nonneg = TRUE)
  nll   <- sum(mu) - sum(y * log(mu))                # Poisson NLL (drop constant)
  Dtheta <- D %*% theta
  pen   <- lam_p * CVXR::cvxr_norm(Dtheta, 1)
  if (!is.null(taper_weights)) {
    if (length(taper_weights) != m)
      stop(sprintf("gi_build_tf_problem: taper_weights must have length %d (rows of D).", m))
    pen <- pen + gam_p * CVXR::cvxr_norm(taper_weights * Dtheta, 1)
  }
  # Tail equality constraints (linear/constant) on the trailing block of theta.
  constraints <- list()
  if (tail == "constant" && n_tail_const > 1L) {
    n_tc <- min(as.integer(n_tail_const), np); base_idx <- np - n_tc + 1L
    for (i in (base_idx + 1L):np)
      constraints <- c(constraints, list(theta[i] == theta[base_idx]))
  } else if (tail == "linear" && n_tail_const >= 3L) {
    n_tc <- min(as.integer(n_tail_const), np); base_idx <- np - n_tc + 1L
    for (i in (base_idx + 2L):np)
      constraints <- c(constraints, list(theta[i] - 2*theta[i - 1L] + theta[i - 2L] == 0))
  }
  prob  <- CVXR::Problem(CVXR::Minimize(nll + pen), constraints)

  list(prob = prob, theta = theta, lam_param = lam_p, gamma_param = gam_p,
       Z = Z, D = D, y = y, off = off,
       np_full = np_full, np_active = np, active_cols = keep_cols, k_diff = k_diff,
       tail = tail, n_tail_const = if (tail == "constant") n_tc else NA_integer_,
       has_taper = !is.null(taper_weights), taper_weights = taper_weights)
}


# ------------------------------------------------------------------------------
# Per-D-row taper weights for the trailing edge.
#   gap(i) = m - i   (0 at the last D row, grows leftward)
#   w(i)   = 1/F_pi_EY(gap+1)  for gap < d_max,  else 0
# Caps the last row at 2x the penultimate to keep ECOS conditioning sane.
# ------------------------------------------------------------------------------
gi_build_tf_taper_weights <- function(design, pi_EY, k_diff = 4L,
                                      active_cols = NULL, cap_gap1_at_2x = TRUE) {
  if (is.null(active_cols)) {
    cs <- colSums(abs(design$Z_full))
    if (max(cs) <= 0) stop("gi_build_tf_taper_weights: empty Z.")
    a <- which(cs > max(cs) * 1e-10); active_cols <- min(a):max(a)
  }
  np <- length(active_cols); m <- np - as.integer(k_diff)
  if (m <= 0L) return(numeric(0))
  F_cdf <- cumsum(pi_EY); d_max <- which(F_cdf >= 1 - 1e-3)[1]
  if (is.na(d_max)) d_max <- length(pi_EY)
  w <- numeric(m)
  for (i in seq_len(m)) {
    gap <- m - i
    if (gap < d_max) {
      F_val <- F_cdf[min(gap + 1L, length(F_cdf))]
      w[i] <- if (F_val > 1e-12) 1 / F_val else 0
    }
  }
  if (isTRUE(cap_gap1_at_2x) && m >= 2L && w[m] > 2 * w[m - 1L] && w[m - 1L] > 0)
    w[m] <- 2 * w[m - 1L]
  w
}


# ------------------------------------------------------------------------------
# Solve at one lambda.  ECOS first; SCS fallback when ECOS chokes.  A sanity
# guard rejects diverged solutions (Rt > 5 or median < 0.3) -- so the fallback
# can only help, never silently corrupt.
# ------------------------------------------------------------------------------
gi_solve_tf <- function(prob_state, lam, gamma = 0, solver = "ECOS",
                        fallback_solver = "SCS", verbose = FALSE) {
  CVXR::value(prob_state$lam_param) <- lam
  if (!is.null(prob_state$gamma_param)) CVXR::value(prob_state$gamma_param) <- gamma
  MAX_RT_SANE <- 5; MIN_MEDIAN_RT <- 0.3

  try_solver <- function(slv) {
    r <- tryCatch(CVXR::solve(prob_state$prob, solver = slv, verbose = verbose),
                  error = function(e) list(status = paste0("error: ", conditionMessage(e))))
    th <- tryCatch(as.numeric(r$getValue(prob_state$theta)),
                   error = function(e) rep(NA_real_, prob_state$np_active))
    ok <- length(th) == prob_state$np_active && !any(is.na(th)) &&
          max(abs(th)) <= MAX_RT_SANE && stats::median(th) >= MIN_MEDIAN_RT
    list(res = r, th = th, ok = ok)
  }
  att <- try_solver(solver)
  if (!att$ok && !is.null(fallback_solver) && !identical(fallback_solver, solver)) {
    att2 <- try_solver(fallback_solver); if (att2$ok) att <- att2 else att$th <- rep(NA_real_, prob_state$np_active)
  } else if (!att$ok) att$th <- rep(NA_real_, prob_state$np_active)

  theta_full <- rep(NA_real_, prob_state$np_full)
  if (!any(is.na(att$th))) theta_full[prob_state$active_cols] <- att$th
  mu <- if (!any(is.na(att$th))) pmax(as.numeric(prob_state$Z %*% att$th) + prob_state$off, 1e-10) else NA_real_
  list(theta = theta_full, theta_active = att$th, mu = mu, lam = lam, status = att$res$status)
}


# ------------------------------------------------------------------------------
# 5-fold CV over the lambda grid.  Pins the active-column block on the FULL
# design so every fold + final refit solves a problem of the same dimension.
# Returns lam_min + lam_1se (largest lam within 1 SE of min, lam >= lam_min).
# ------------------------------------------------------------------------------
gi_select_lambda_cv_tf <- function(design, lam_grid, k_diff = 4L, nfold = 5L,
                                   seed = 1L, solver = "ECOS", verbose = FALSE) {
  Z <- design$Z_full; y <- as.numeric(design$Y_valid); off <- as.numeric(design$off_valid)
  col_sums <- colSums(abs(Z)); tol_col <- max(col_sums) * 1e-10
  active   <- which(col_sums > tol_col); if (!length(active)) stop("CV: empty Z.")
  active_cols <- min(active):max(active); np_act <- length(active_cols)
  n <- nrow(Z); set.seed(seed); fold_id <- sample(rep(seq_len(nfold), length.out = n))
  fold_dev <- matrix(NA_real_, nfold, length(lam_grid))

  for (kf in seq_len(nfold)) {
    test  <- which(fold_id == kf); train <- which(fold_id != kf)
    if (length(train) < np_act / 2L) next
    d_tr <- design; d_tr$Z_full <- Z[train, , drop = FALSE]
    d_tr$Y_valid <- y[train]; d_tr$off_valid <- off[train]
    prob_tr <- gi_build_tf_problem(d_tr, k_diff = k_diff, active_cols = active_cols)
    cache <- vector("list", length(lam_grid))
    for (li in seq_along(lam_grid)) {
      f <- gi_solve_tf(prob_tr, lam_grid[li], solver = solver, verbose = FALSE)
      if (!any(is.na(f$theta_active))) cache[[li]] <- f$theta_active
    }
    Zte <- Z[test, active_cols, drop = FALSE]; yte <- y[test]; offte <- off[test]
    for (li in seq_along(lam_grid)) {
      th <- cache[[li]]; if (is.null(th)) next
      mu_te <- pmax(as.numeric(Zte %*% th) + offte, 1e-10)
      dev_t <- 2 * ifelse(yte > 0, yte * log(yte / mu_te), 0) - 2 * (yte - mu_te)
      fold_dev[kf, li] <- mean(dev_t)
    }
    if (verbose) cat(sprintf("  TF CV fold %d/%d done\n", kf, nfold))
  }
  cv_mean <- colMeans(fold_dev, na.rm = TRUE)
  cv_sd   <- apply(fold_dev, 2L, sd, na.rm = TRUE)
  cv_se   <- cv_sd / sqrt(pmax(colSums(!is.na(fold_dev)), 1L))
  best_i_min <- which.min(cv_mean); best_lam_min <- lam_grid[best_i_min]
  cv_thresh  <- cv_mean[best_i_min] + cv_se[best_i_min]
  cand <- which(cv_mean <= cv_thresh & lam_grid >= best_lam_min)
  if (!length(cand)) { best_i_1se <- best_i_min; best_lam_1se <- best_lam_min
  } else { best_i_1se <- cand[which.max(lam_grid[cand])]; best_lam_1se <- lam_grid[best_i_1se] }
  list(lam_grid = lam_grid, cv_mean = cv_mean, cv_sd = cv_sd, cv_se = cv_se,
       best_i_min = best_i_min, best_lam_min = best_lam_min,
       best_i_1se = best_i_1se, best_lam_1se = best_lam_1se,
       cv_threshold_1se = cv_thresh, nfold = nfold, seed = seed, active_cols = active_cols)
}


# ------------------------------------------------------------------------------
# theta -> Rt data frame.  No CIs (would need a TF-specific inference story).
# ------------------------------------------------------------------------------
gi_extract_rt_tf <- function(fit, design = NULL) {
  np <- length(fit$theta)
  data.frame(day = seq.int(0L, np - 1L), Rt_mean = as.numeric(fit$theta),
             Rt_lo = NA_real_, Rt_hi = NA_real_)
}


# ------------------------------------------------------------------------------
# Closed-form single-pass DoW from anchor TF residuals.  d_j = sum(Y|j)/sum(mu|j);
# auto-balanced to mean(weighted) = 1 because the anchor matches total Y.
# Pre-applies the DoW to design rows so downstream CV / final fits run unchanged.
# ------------------------------------------------------------------------------
.gi_tf_estimate_dow <- function(design, lam_anchor, k_diff, solver = "ECOS") {
  if (is.null(design$dow_valid) || length(design$dow_valid) != nrow(design$Z_full))
    stop(".gi_tf_estimate_dow: design$dow_valid missing or wrong length.")
  prob <- gi_build_tf_problem(design, k_diff = k_diff)
  fit  <- gi_solve_tf(prob, lam_anchor, solver = solver)
  if (any(is.na(fit$theta_active))) stop(".gi_tf_estimate_dow: anchor solve failed.")
  mu0 <- as.numeric(fit$mu); y <- as.numeric(design$Y_valid)
  dow <- as.integer(design$dow_valid)
  d_vec <- vapply(seq_len(7L), function(j) {
    idx <- which(dow == j); if (!length(idx)) return(1.0)
    sm <- sum(mu0[idx]); if (!is.finite(sm) || sm <= 0) return(1.0)
    sum(y[idx]) / sm
  }, numeric(1))
  setNames(d_vec, c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"))
}


# ------------------------------------------------------------------------------
# Public retrospective TF wrapper.  Returns curves at lam_min and lam_1se.
# Optional single-pass DoW (use_dow = TRUE; one extra anchor solve + closed-form
# d_j, then pre-applied to design rows).
# ------------------------------------------------------------------------------
fit_girt_tf_retrospective <- function(design, lam_grid, k_diff = 4L,
                                      nfold = 5L, seed = 1L,
                                      use_dow = FALSE, dow_anchor_lam = NULL,
                                      solver = "ECOS", verbose = FALSE) {
  if (!identical(design$basis, "tf_identity"))
    stop("fit_girt_tf_retrospective: design must come from build_gi_design_tf().")

  dow_effects <- NULL
  if (isTRUE(use_dow)) {
    if (is.null(design$dow_valid)) stop("use_dow = TRUE needs dow_dates in build_gi_design_tf().")
    if (is.null(dow_anchor_lam)) dow_anchor_lam <- stats::median(lam_grid)
    dow_effects <- .gi_tf_estimate_dow(design, dow_anchor_lam, k_diff, solver)
    d_t <- as.numeric(dow_effects[design$dow_valid])
    design$Z_full <- design$Z_full * d_t; design$off_valid <- design$off_valid * d_t
    design$use_dow <- FALSE; design$dow_valid <- NULL
    if (verbose) cat(sprintf("  TF DoW (anchor %.3e): %s\n", dow_anchor_lam,
                             paste(sprintf("%s=%.2f", names(dow_effects), dow_effects), collapse = " ")))
  }
  sel <- gi_select_lambda_cv_tf(design, lam_grid, k_diff = k_diff, nfold = nfold,
                                seed = seed, solver = solver, verbose = verbose)
  prob_full <- gi_build_tf_problem(design, k_diff = k_diff, active_cols = sel$active_cols)
  fit_min <- gi_solve_tf(prob_full, sel$best_lam_min, solver = solver)
  fit_1se <- gi_solve_tf(prob_full, sel$best_lam_1se, solver = solver)
  list(rt_min = gi_extract_rt_tf(fit_min), rt_1se = gi_extract_rt_tf(fit_1se),
       lam_min = sel$best_lam_min, lam_1se = sel$best_lam_1se,
       sel = sel, dow_effects = dow_effects,
       config = list(basis = "tf_identity", k_diff = k_diff, nfold = nfold, seed = seed,
                     use_dow = isTRUE(use_dow), dow_anchor_lam = if (use_dow) dow_anchor_lam else NA_real_))
}


# ------------------------------------------------------------------------------
# Forward-validation gamma tuner for the realtime tapered-TF path.
#
# For each truncation s in {t_current - n_fv, ..., t_current - 1}:
#   1. build a TF design on data through day s (taper weights recomputed),
#   2. for each gamma in gamma_grid, fit at fixed lam, linearly extrapolate
#      the last two theta values to predict R_{s+1}, predict Y_{s+1}, score |.|.
# Inputs include the renewal Force at the prediction day, the severity, and
# the (pre-applied) DoW effect at t = s+1, so the predicted Y matches the
# observation likelihood the spline FV uses.
# Returns best_gamma_min (default rule) + best_gamma_1se.
# ------------------------------------------------------------------------------
gi_tune_gamma_tf_fv <- function(E_star, g, pi_EY, inf_inc, dates, t_current,
                                lam, gamma_grid, max_time,
                                first_rt_index = 0L, severity_rate = 1,
                                dow_dates = NULL, dow_effects = NULL,
                                n_fv = 7L, k_diff = 4L,
                                tail = c("linear", "none", "constant"),
                                n_tail_const = 7L, knot_step = 5L,
                                solver = "ECOS", verbose = FALSE) {
  tail <- match.arg(tail)
  s_vals <- seq.int(t_current - n_fv, t_current - 1L)
  designs  <- vector("list", length(s_vals))
  taper_ws <- vector("list", length(s_vals))
  for (j in seq_along(s_vals)) {
    s <- s_vals[j]
    if (s + 1L > length(E_star)) next
    E_s <- E_star[seq_len(s + 1L)]; inc_s <- as.numeric(inf_inc)[seq_len(s + 1L)]
    dts_s <- as.Date(dates)[seq_len(s + 1L)]
    dow_s <- if (!is.null(dow_dates)) as.Date(dow_dates)[seq_len(s + 1L)] else NULL
    d_s <- tryCatch(build_gi_design_tf(
      E_star = E_s, g = g, pi_EY = pi_EY, max_time = s, inf_inc = inc_s,
      knot_step = knot_step, first_rt_index = first_rt_index,
      dates = dts_s, dow_dates = dow_s), error = function(e) NULL)
    if (is.null(d_s)) next
    sev_by_day <- if (length(severity_rate) == 1L) rep(severity_rate, s + 1L) else severity_rate[seq_len(s + 1L)]
    d_s <- gi_apply_severity_to_design(d_s, sev_by_day)
    if (!is.null(dow_effects) && !is.null(d_s$dow_valid)) {
      d_t <- as.numeric(dow_effects[d_s$dow_valid])
      d_s$Z_full <- d_s$Z_full * d_t; d_s$off_valid <- d_s$off_valid * d_t
      d_s$use_dow <- FALSE; d_s$dow_valid <- NULL
    }
    designs[[j]]  <- d_s
    taper_ws[[j]] <- gi_build_tf_taper_weights(d_s, pi_EY, k_diff = k_diff)
  }
  # Prediction-side multipliers at t = s+1: renewal force D_s = sum_k g_k X_hat[s+1-k] * pi_EY[k]... actually
  # the per-day kernel sum at the prediction step.  Use the row of the un-aggregated kernel that the spline
  # FV uses: D_pred = sum_s Lambda_s * pi_EY[t-s] for the predicted day, i.e. the row-sum of the X_ts matrix
  # at row t = s+1 -- equivalently the renewal forecast on inf_inc.  We approximate using Lambda directly.
  Lambda <- gi_renewal_force(E_star, g); d_max <- length(pi_EY)
  D_vals  <- numeric(length(s_vals)); sev_vals <- numeric(length(s_vals)); dow_vals <- numeric(length(s_vals))
  Y_true  <- numeric(length(s_vals))
  for (j in seq_along(s_vals)) {
    s <- s_vals[j]; t_pred <- s + 1L
    if (t_pred + 1L > length(inf_inc)) { Y_true[j] <- NA_real_; next }
    Dv <- 0
    for (k in seq_len(min(d_max, t_pred)))
      if ((t_pred - k + 1L) >= 1L && (t_pred - k + 1L) <= length(Lambda))
        Dv <- Dv + Lambda[t_pred - k + 1L] * pi_EY[k]
    D_vals[j]   <- Dv
    sev_vals[j] <- if (length(severity_rate) == 1L) severity_rate else severity_rate[t_pred + 1L]
    dow_vals[j] <- if (!is.null(dow_effects) && !is.null(dow_dates) && (t_pred + 1L) <= length(dow_dates))
                     as.numeric(dow_effects[as.integer(format(as.Date(dow_dates[t_pred + 1L]), "%u"))])
                   else 1.0
    Y_true[j]   <- as.numeric(inf_inc[t_pred + 1L])
  }

  step_err <- matrix(NA_real_, length(s_vals), length(gamma_grid))
  for (j in seq_along(s_vals)) {
    d_s <- designs[[j]]; if (is.null(d_s)) next
    prob_s <- gi_build_tf_problem(d_s, k_diff = k_diff, tail = tail,
                                  n_tail_const = n_tail_const, taper_weights = taper_ws[[j]])
    for (gi in seq_along(gamma_grid)) {
      fit_sg <- gi_solve_tf(prob_s, lam = lam, gamma = gamma_grid[gi], solver = solver)
      th_full <- fit_sg$theta
      if (is.null(th_full) || all(is.na(th_full))) next
      idx_act <- which(!is.na(th_full)); if (length(idx_act) < 2L) next
      last2 <- utils::tail(idx_act, 2L)
      Rt_pred <- th_full[last2[2]] + (th_full[last2[2]] - th_full[last2[1]]) /
                  max(last2[2] - last2[1], 1L)
      Y_hat  <- Rt_pred * D_vals[j] * sev_vals[j] * dow_vals[j]
      step_err[j, gi] <- abs(Y_hat - Y_true[j])
    }
    if (verbose) cat(sprintf("  TF FV fold s=%d done\n", s_vals[j]))
  }
  fv_mean <- apply(step_err, 2L, function(x) { x <- x[is.finite(x)]; if (length(x)) mean(x) else NA_real_ })
  fv_n    <- apply(step_err, 2L, function(x) sum(is.finite(x)))
  fv_sd   <- apply(step_err, 2L, function(x) { x <- x[is.finite(x)]; if (length(x) > 1L) sd(x) else NA_real_ })
  fv_se   <- ifelse(fv_n > 0L, fv_sd / sqrt(fv_n), NA_real_)
  best_i_min <- if (any(is.finite(fv_mean))) which.min(fv_mean) else 1L
  thresh <- fv_mean[best_i_min] + ifelse(is.finite(fv_se[best_i_min]), fv_se[best_i_min], 0)
  cand <- which(is.finite(fv_mean) & fv_mean <= thresh)
  best_i_1se <- if (length(cand)) max(cand) else best_i_min
  list(gamma_grid = gamma_grid, fv_mean = fv_mean, fv_sd = fv_sd, fv_se = fv_se,
       step_errors = step_err,
       best_i_min = best_i_min, best_gamma_min = gamma_grid[best_i_min],
       best_i_1se = best_i_1se, best_gamma_1se = gamma_grid[best_i_1se],
       best_gamma = gamma_grid[best_i_min], s_vals = s_vals, Y_true = Y_true)
}
