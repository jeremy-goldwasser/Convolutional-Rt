# MechRt :: taper.R
#
# Real-time machinery: tapered smoothing penalty + null-space IRLS solver +
# forward-validation gamma tuning.
#
#   - gi_build_tapered_penalty   construct P_taper (CDF-based tail weights)
#   - .gi_solve_taper    null-space IRLS solver.  Use this whenever
#                             lam_taper > 0 (taper requires it for stable
#                             solves at large lam_taper).  See fit.R header
#                             for the choice between gi_solve and this.
#   - gi_tune_gamma_fv           forward-validation gamma tuning.
#
# fit_mechrt_retrospective lives in wrappers.R.  A flu-specific wrapper
# fv_tune_gamma_realtime lives in real/flu/realtime_lib.R.
# ==============================================================================
# TAPERED SMOOTHING PENALTY (tail regularisation)
#
# The tapered penalty penalizes weighted squared first differences of beta(t):
#   sum_t  w_t * (beta_{t+1} - beta_t)^2
# where w_t = 1/F(T - t) for t >= T - d, 0 otherwise.
# F is the CDF of the E-to-I delay distribution (pi_E).
# ==============================================================================

# -----------------------------------------------------------------------------
# Real-time tapered-smoothing penalty (CDF-based).
#
# Penalty acts on first differences of TOTAL beta(t):
#
#   P_taper = D1' diag(w) D1
#
# Weights w(gap), where gap = T_end - t (days from real-time edge):
#
#   w(gap) = 1 / F(gap)   for  1 <= gap <= d
#   w(gap) = 0            for  gap = 0 or gap > d
#
# where F = cumsum(pi_E) is the CDF of the E-to-observed-outcome delay
# (for hosp pipelines: pi_EH) and d = smallest k with F(k) >= 1 - 1e-3 is
# the effective support.
#
# IMPORTANT: pi_E must be the E-to-OBSERVED-OUTCOME distribution (e.g. pi_EH
# for flu hosp pipelines).  Callers that pass the E-to-I latent delay
# by mistake will construct a penalty with the wrong support.
# -----------------------------------------------------------------------------
gi_build_tapered_penalty <- function(design, pi_E, quiet = FALSE) {
  # Capping option dead behind this hardcode (see REMOVED_OPTIONS.md H)
  cap_gap1_at_2x <- FALSE
  nb       <- design$n_basis
  B        <- design$B
  use_kink <- isTRUE(design$use_kink)
  kink_vec <- design$kink_vec

  G    <- if (use_kink) cbind(B, kink_vec) else B
  D1_G <- diff(G)   # max_time x np

  F_cdf <- cumsum(pi_E)
  d <- which(F_cdf >= 1 - 1e-3)[1]
  if (is.na(d)) d <- length(pi_E)

  T_end <- design$fr[2]   # anchor taper at spline domain boundary, NOT max_time
  w <- numeric(nrow(D1_G))
  for (i in seq_len(nrow(D1_G))) {
    t_0based <- i - 1L
    gap <- T_end - t_0based
    # gap > 0 (not >=): exclude the D1_G row at gap=0 (B(T_end+1) - B(T_end));
    # B(T_end+1) lies outside the spline domain and would pull theta[last] to 0.
    if (gap > 0 && gap <= d) {
      F_val <- F_cdf[min(gap, length(F_cdf))]
      w[i]  <- if (F_val > 1e-12) 1.0 / F_val else 0.0
    }
  }

  # Optional: cap the gap=1 weight at 2 * gap=2 weight to avoid the 1/F(1)
  # outlier (which can be 100x+ the next weight) dominating the taper.
  if (isTRUE(cap_gap1_at_2x)) {
    i1 <- which(T_end - (seq_along(w) - 1L) == 1L)
    i2 <- which(T_end - (seq_along(w) - 1L) == 2L)
    if (length(i1) == 1L && length(i2) == 1L && w[i1] > 2 * w[i2]) {
      w[i1] <- 2 * w[i2]
    }
  }

  if (!quiet)
    cat(sprintf("  Tapered penalty (1/F, d=%d support, cap_gap1_at_2x=%s): %d non-zero weights (t = %d..%d)\n",
                d, isTRUE(cap_gap1_at_2x), sum(w > 0), T_end - d, T_end - 1L))

  P_taper <- crossprod(D1_G * sqrt(w))   # np x np

  list(P_taper = P_taper, weights = w, D1_G = D1_G)
}


# ==============================================================================
# CORE IRLS SOLVER — null-space method (no KKT zero block)
#
# Constraints A theta = 0 are handled by projecting into null(A):
#   theta = N alpha,  where N is a basis for null(A)
#   Solve: (N'MN) alpha = N' ZtWyo   — always positive definite, no zero block
#
# When lam_taper = 0 and P_taper = NULL, produces same results as gi_solve.
# ==============================================================================

.gi_solve_taper <- function(
    design,
    lam,
    lam_taper       = 0,
    P_taper         = NULL,
    tail            = "linear",
    theta_init      = NULL,
    max_iter        = 300L,
    tol             = 1e-8,
    verbose         = FALSE
) {
  # ridge-penalty branch dead behind this hardcode (see REMOVED_OPTIONS.md D)
  penalty <- "pspline"

  if (identical(design$basis, "natural") && tail == "constant") {
    stop("MechRt: basis = 'natural' with tail = 'constant' is deprecated. ",
         "Use basis = 'regression' with tail = 'constant' for an explicit ",
         "constant-tail constraint, or stick with the smoothing spline ",
         "(basis = 'natural') and leave the tail unconstrained.")
  }

  Z   <- design$Z_full
  y   <- design$Y_valid
  off <- design$off_valid
  np  <- design$np
  nb  <- design$n_basis

  n_obs_t     <- if (!is.null(design$n_obs))     design$n_obs     else length(y)
  n_penalty_t <- if (!is.null(design$n_penalty)) design$n_penalty else max(nb - 2L, 1L)
  P_base <- if (penalty == "pspline") design$P_pspline else design$P_ridge
  has_taper <- (lam_taper > 0) && !is.null(P_taper)
  P_total <- lam * (P_base / n_penalty_t)
  if (has_taper) P_total <- P_total + lam_taper * (P_taper / n_penalty_t)
  if (!is.null(design$head_const_extra_M) && isTRUE(design$head_const_days > 0L)) {
    P_total <- P_total + design$head_const_extra_M
  }

  A    <- design$tail_constraints[[tail]]
  n_eq <- nrow(A)

  # --- Null-space basis for constraints A theta = 0 ---------------------------
  # N is np x (np - n_eq), columns span null(A).
  # theta = N %*% alpha satisfies A theta = 0 automatically.
  if (n_eq > 0) {
    qr_A <- qr(t(A))
    N_null <- qr.Q(qr_A, complete = TRUE)[, (n_eq + 1L):np, drop = FALSE]
  } else {
    N_null <- diag(np)
  }
  n_free <- ncol(N_null)

  # --- Initialise theta -------------------------------------------------------
  if (!is.null(theta_init) && length(theta_init) == np) {
    theta <- theta_init
  } else {
    rs <- rowSums(Z[, seq_len(nb), drop = FALSE])
    rs_pos <- rs[rs > 0]
    ref <- if (length(rs_pos) > 0) median(rs_pos) else 1.0
    init_val <- max(1e-6, (mean(y) - mean(off)) / ref)
    theta <- c(rep(init_val, nb), rep(0.0, np - nb))
  }
  # Project initial theta onto null(A)
  # alpha_init = N' theta (least-squares projection), theta = N alpha
  alpha <- as.numeric(crossprod(N_null, theta))
  theta <- as.numeric(N_null %*% alpha)

  # Precompute Z_N = Z %*% N  (n_valid x n_free)
  Z_N <- Z %*% N_null
  # Precompute P_total_N = N' P_total N  (n_free x n_free)
  P_total_N <- crossprod(N_null, P_total %*% N_null)

  # --- IRLS iterations --------------------------------------------------------
  llik_prev <- -Inf
  for (iter in seq_len(max_iter)) {
    mu <- pmax(as.numeric(Z %*% theta) + off, 1e-10)
    w  <- 1.0 / mu   # Poisson working weights (identity link)

    # Reduced normal equations in null-space coordinates:
    #   (Z_N' W Z_N + P_total_N + eps I) alpha = Z_N' W (y - off)
    ZNtWZN  <- crossprod(Z_N * sqrt(w)) / n_obs_t    # n_free x n_free (normalized)
    ZNtWyo  <- as.numeric(crossprod(Z_N, w * (y - off))) / n_obs_t  # n_free x 1

    eps_ridge <- max(1e-10 * max(diag(ZNtWZN), 1.0), .Machine$double.eps)
    M_red <- ZNtWZN + P_total_N + eps_ridge * diag(n_free)

    alpha_new <- tryCatch(
      solve(M_red, ZNtWyo),
      error = function(e) {
        tryCatch(qr.solve(M_red, ZNtWyo), error = function(e2) alpha)
      }
    )
    theta_new <- as.numeric(N_null %*% alpha_new)

    # Backtracking line search: keep mu > 0
    step <- theta_new - theta
    ss   <- 1.0
    for (k in seq_len(30L)) {
      mu_try <- as.numeric(Z %*% (theta + ss * step)) + off
      if (all(mu_try > 0)) break
      ss <- ss * 0.5
    }
    theta <- theta + ss * step
    alpha <- as.numeric(crossprod(N_null, theta))

    # Convergence check
    mu_new <- pmax(as.numeric(Z %*% theta) + off, 1e-10)
    llik   <- sum(y * log(mu_new) - mu_new)
    delta  <- abs(llik - llik_prev) / (1.0 + abs(llik_prev))
    if (verbose) cat(sprintf("  iter %3d: llik = %10.4f  delta = %.2e\n", iter, llik, delta))
    if (iter > 5L && delta < tol) break
    llik_prev <- llik
  }

  mu_final   <- pmax(as.numeric(Z %*% theta) + off, 1e-10)
  llik_final <- sum(y * log(mu_final) - mu_final)

  if (verbose) cat(sprintf("  Converged in %d iterations.  llik = %.4f\n", iter, llik_final))

  list(
    theta           = theta,
    mu              = mu_final,
    llik            = llik_final,
    n_iter          = iter,
    penalty         = penalty,
    tail_constraint = tail,
    lam             = lam,
    lam_taper       = lam_taper
  )
}


# ==============================================================================
# fit_mechrt_retrospective lives in wrappers.R (smoothing-spline default,
# CV-tuned lambda, no tail constraint).
# ==============================================================================


# ==============================================================================
# FORWARD VALIDATION: tune gamma (lam_taper) via 1-step-ahead prediction
#
# For each candidate gamma, for s = (t-n_fv):(t-1):
#   1. Build design with data through s
#   2. Fit with (lambda, gamma) taper
#   3. Linearly extrapolate last two fitted betas to get beta_hat(s+1)
#   4. Predict Y_hat_{s+1} = beta_hat(s+1) * D_{s+1}
#   5. Compare to true Y_{s+1}
# FV score = MAE across all s values.  Pick gamma with minimum FV score.
# ==============================================================================

# Forward-validation gamma tuner for girt.  Takes the reporting-delay (E->Y)
# Gamma parameters; the per-fold D_{s+1} kernel is built from those when the
# fold_builder doesn't supply its own pi_E.
gi_tune_gamma_fv <- function(
    time_series,      # full time_series (NOT truncated)
    inf_inc,          # full infectious incidence vector
    report_shape,     # shape parameter of E-to-Y (reporting) Gamma delay
    report_rate,      # rate parameter of E-to-Y (reporting) Gamma delay
    t_current,        # current time (e.g. 110)
    lam_tuned,        # lambda from retrospective tuning
    retro_theta,      # retrospective theta for warm-starting
    gamma_grid,       # candidate gamma (lam_taper) values
    n_fv          = 7L,         # number of forward validation folds
    knot_spacing  = 10L,
    y_min_count   = 0,
    tail          = "linear",
    verbose       = TRUE,
    first_rt_index    = NULL,
    dates             = NULL,
    spline_degree          = 3L,    # cubic, matches build_design default
    basis                  = c("natural", "regression"),
    knot_step              = 5L,    # canonical natural-spline density (matches build_design)
    severity_rate          = 1.0,     # scalar or length-(t_current+2) vector
    dow_effects            = NULL,    # numeric(7), Mon..Sun, with prod == 1
    dow_dates              = NULL,    # Date vector covering at least t_current+2 entries
    # Optional callback for per-fold design construction.  When supplied,
    # fold_builder(s) must return a list with at least
    #   $design       : a build_design() output (possibly post-processed)
    #   $time_series  : a data.frame with $pre_S and $pre_I usable for the
    #                   D_{s+1} computation (length >= s+2)
    # This is the hook used by real-data pipelines that need per-fold state
    # re-imputation, severity rescaling, likelihood-start filtering, and
    # DoW pre-scaling before the gamma-taper FV fit.
    fold_builder           = NULL
) {
  basis <- match.arg(basis)

  # E0 bolus offset disabled (section F).  Hardcoded so the conditional
  # dead branch at line ~395 still parses; revivable from REMOVED_OPTIONS.md.
  E0_init <- 0

  .t0_fv <- Sys.time()
  k_fv <- n_fv  # number of 1-step-ahead predictions

  # s ranges from (t_current - k_fv) to (t_current - 1)
  s_vals <- seq(t_current - k_fv, t_current - 1L)
  cat(sprintf("  Forward validation: s = %d..%d, predicting Y_{s+1} = Y_%d..Y_%d\n",
              s_vals[1], s_vals[length(s_vals)],
              s_vals[1] + 1L, s_vals[length(s_vals)] + 1L))

  # Renewal force Lambda(t) for D_{s+1} computation (full time_series path).
  # time_series$pre_I = Lambda; pre_S column retained for back-compat but unused.
  I_vec_full <- as.numeric(time_series$pre_I)
  # (pi_E computed per-s below using report_shape/report_rate)

  # ---- Build designs for each holdout s (reused across gamma grid) -----------
  cat("  Building designs for each holdout s...\n")
  designs    <- list()
  fold_tss   <- list()   # per-fold time_series (for D_{s+1} computation)
  fold_pi_Es <- list()   # per-fold pi_E pmf, if fold_builder supplies one
  for (s in s_vals) {
    key <- as.character(s)
    if (!is.null(fold_builder)) {
      fold <- tryCatch(
        fold_builder(s),
        error = function(e) {
          cat(sprintf("    s=%d fold_builder failed: %s\n", s, conditionMessage(e)))
          NULL
        }
      )
      if (is.null(fold) || is.null(fold$design)) {
        designs[[key]]    <- NULL
        fold_tss[[key]]   <- NULL
        fold_pi_Es[[key]] <- NULL
      } else {
        designs[[key]]    <- fold$design
        fold_tss[[key]]   <- fold$time_series
        fold_pi_Es[[key]] <- fold$pi_E   # may be NULL: fall back below
      }
    } else {
      # girt always tunes gamma through a fold_builder (the GI per-fold design
      # is built from truncated exposures + the generation interval).  The
      # mechanistic build_design default path is intentionally removed here.
      stop("gi_tune_gamma_fv: a fold_builder must be supplied (girt has no ",
           "mechanistic default design path).")
    }
  }

  # ---- Compute D_{s+1} for each s (total exposure at s+1) -------------------
  # girt Rt mode: I_vec_j IS the renewal force Lambda(t) = sum_k g_k X(t-k), so
  #   D_{s+1} = sum_{k=1}^{d_max} Lambda_{s+1-k} * pi_EY[k].
  D_vals   <- numeric(length(s_vals))
  off_vals <- numeric(length(s_vals))
  sev_vals <- numeric(length(s_vals))
  dow_vals <- numeric(length(s_vals))
  for (j in seq_along(s_vals)) {
    s <- s_vals[j]
    t_pred <- s + 1L  # 0-based time we're predicting
    ts_j <- fold_tss[[as.character(s)]]
    I_vec_j <- if (is.null(ts_j)) I_vec_full else as.numeric(ts_j$pre_I)
    # Prefer per-fold pi_E if fold_builder supplied one; otherwise reconstruct
    # from report_shape / report_rate (sim path).  This ensures the D_{s+1}
    # convolution uses the EXACT same kernel that was used to build the design.
    pi_E_fold <- fold_pi_Es[[as.character(s)]]
    pi_E_pred <- if (!is.null(pi_E_fold)) as.numeric(pi_E_fold)
                 else diff(pgamma(0:(t_pred + 1L),
                                  shape = report_shape, rate = report_rate))
    d_max_s <- length(pi_E_pred)
    D <- 0.0
    for (k in seq_len(min(d_max_s, t_pred))) {
      idx <- t_pred - k + 1L  # 1-based index into I_vec_j (the renewal force)
      if (idx >= 1 && idx <= length(I_vec_j))
        D <- D + I_vec_j[idx] * pi_E_pred[k]
    }
    D_vals[j] <- D
    # E0 offset at t_pred
    off_vals[j] <- if (E0_init > 0 && t_pred <= d_max_s) E0_init * pi_E_pred[t_pred] else 0.0

    # Severity multiplier at t_pred (1-based -> vector index t_pred + 1)
    sev_vals[j] <- if (length(severity_rate) == 1L) {
      as.numeric(severity_rate)
    } else if ((t_pred + 1L) <= length(severity_rate)) {
      as.numeric(severity_rate[t_pred + 1L])
    } else 1.0

    # DoW multiplier at t_pred
    if (!is.null(dow_effects) && !is.null(dow_dates) &&
        (t_pred + 1L) <= length(dow_dates)) {
      s1_date <- as.Date(dow_dates[t_pred + 1L])
      if (!is.na(s1_date)) {
        dow_idx <- as.integer(format(s1_date, "%u"))  # 1..7 Mon..Sun
        dow_vals[j] <- if (dow_idx >= 1L && dow_idx <= length(dow_effects))
                         as.numeric(dow_effects[dow_idx]) else 1.0
      } else {
        dow_vals[j] <- 1.0
      }
    } else {
      dow_vals[j] <- 1.0
    }
  }

  Y_true <- inf_inc[s_vals + 2L]  # Y at s+1 (1-based indexing: inf_inc[s+1+1])

  # ---- Sweep gamma grid ------------------------------------------------------
  mae_vals <- rep(NA_real_, length(gamma_grid))
  # Per-step absolute errors: rows = s_vals (FV folds), cols = gamma_grid
  step_err <- matrix(NA_real_, nrow = length(s_vals), ncol = length(gamma_grid))

  for (gi in seq_along(gamma_grid)) {
    gam <- gamma_grid[gi]
    gam_label <- if (gam == 0) "0" else sprintf("%.1e", gam)
    if (verbose) cat(sprintf("  gamma = %s: ", gam_label))

    for (j in seq_along(s_vals)) {
      s <- s_vals[j]
      d_s <- designs[[as.character(s)]]
      if (is.null(d_s)) next

      # Fit with lambda + gamma taper
      if (gam == 0) {
        fit_s <- tryCatch(
          gi_solve(d_s, lam_tuned, tail = tail),
          error = function(e) NULL
        )
      } else {
        pi_E_s_taper <- diff(pgamma(0:(s + 1L), shape = report_shape, rate = report_rate))
        taper_s <- tryCatch(
          gi_build_tapered_penalty(d_s, pi_E_s_taper, quiet = TRUE),
          error = function(e) NULL
        )
        if (is.null(taper_s)) { fit_s <- NULL } else {
          fit_s <- tryCatch(
            .gi_solve_taper(d_s, lam_tuned, lam_taper = gam, P_taper = taper_s$P_taper,
                         tail = tail, theta_init = retro_theta),
            error = function(e) NULL
          )
        }
      }
      if (is.null(fit_s) || !all(is.finite(fit_s$theta))) next

      # Evaluate fitted curve at times 0:s.
      # In beta mode this is the beta(t) curve; in Rt mode it is R_t.
      G_s <- if (isTRUE(d_s$use_kink)) cbind(d_s$B, d_s$kink_vec) else d_s$B
      curve <- as.numeric(G_s %*% fit_s$theta)

      # Linearly extrapolate last two values to get curve(s+1).
      # curve is indexed 1:(s+1) for times 0:s.
      c_s   <- curve[s + 1L]   # value at s
      c_sm1 <- curve[s]         # value at s-1
      c_ext <- 2.0 * c_s - c_sm1

      # Predict Y_{s+1}.  Severity and DoW multipliers default to 1 so that
      # the simple beta-mode sim callers see the unchanged formula
      #   Y_hat = beta_ext * D + off
      Y_hat <- c_ext * D_vals[j] * sev_vals[j] * dow_vals[j] + off_vals[j]
      step_err[j, gi] <- abs(Y_hat - Y_true[j])
    }

    finite_e <- step_err[is.finite(step_err[, gi]), gi]
    mae <- if (length(finite_e) > 0) mean(finite_e) else NA_real_
    mae_vals[gi] <- mae
    if (verbose) cat(sprintf("MAE = %.2f\n", mae))
  }

  # ---- Per-gamma summaries ---------------------------------------------------
  fv_mean <- apply(step_err, 2L, function(x) {
    x <- x[is.finite(x)]; if (!length(x)) Inf else mean(x)
  })
  fv_sd <- apply(step_err, 2L, function(x) {
    x <- x[is.finite(x)]; if (length(x) <= 1L) NA_real_ else sd(x)
  })
  fv_n  <- apply(step_err, 2L, function(x) sum(is.finite(x)))
  fv_se <- ifelse(fv_n > 0L, fv_sd / sqrt(fv_n), NA_real_)

  best_i_min <- which.min(fv_mean)
  if (!is.finite(fv_mean[best_i_min])) {
    best_i_min <- which.min(mae_vals)
  }

  # 1se rule (more regularization = larger gamma): largest gamma whose
  # fv_mean is within 1 SE of the min
  fv_threshold_1se <- fv_mean[best_i_min] +
    ifelse(is.finite(fv_se[best_i_min]), fv_se[best_i_min], 0.0)
  eligible_1se <- which(is.finite(fv_mean) & fv_mean <= fv_threshold_1se)
  best_i_1se   <- if (length(eligible_1se) > 0L) max(eligible_1se) else best_i_min

  # 2se rule
  fv_threshold_2se <- fv_mean[best_i_min] +
    ifelse(is.finite(fv_se[best_i_min]), 2.0 * fv_se[best_i_min], 0.0)
  eligible_2se <- which(is.finite(fv_mean) & fv_mean <= fv_threshold_2se)
  best_i_2se   <- if (length(eligible_2se) > 0L) max(eligible_2se) else best_i_min

  best_i <- best_i_min  # backward-compat: original behaviour was min-rule
  best_gamma <- gamma_grid[best_i]
  cat(sprintf("  Best gamma: %.3e  (MAE = %.2f)\n", best_gamma, mae_vals[best_i]))

  list(
    gamma_grid       = gamma_grid,
    mae              = mae_vals,
    best_gamma       = best_gamma,
    best_i           = best_i,
    best_mae         = mae_vals[best_i],
    s_vals           = s_vals,
    Y_true           = Y_true,
    # ---- New (additive) fields ---------------------------------------------
    step_errors      = step_err,
    fv_mean          = fv_mean,
    fv_sd            = fv_sd,
    fv_se            = fv_se,
    best_i_min       = best_i_min,
    best_gamma_min   = gamma_grid[best_i_min],
    best_i_1se       = best_i_1se,
    best_gamma_1se   = gamma_grid[best_i_1se],
    fv_threshold_1se = fv_threshold_1se,
    best_i_2se       = best_i_2se,
    best_gamma_2se   = gamma_grid[best_i_2se],
    fv_threshold_2se = fv_threshold_2se,
    elapsed_sec      = as.numeric(Sys.time() - .t0_fv, units = "secs")
  )
}


