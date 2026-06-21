# MechRt :: lambda_select.R
#
# Lambda selection for the smoothing penalty.  Two routes:
#   - gi_select_lambda_cv:      K-fold CV with cv_select_rule "min" or "1se" and
#                            Poisson "deviance" (default) or "mse" / "mae".
#   - gi_select_lambda_cv_dow:  K-fold CV that estimates DoW effects per fold via
#                            alternating coordinate descent.  Used by the
#                            real-data flu pipeline; new code should prefer
#                            build_design(..., dow_dates=...) + gi_select_lambda_cv.
# Both refit candidate lambdas via gi_solve.
#
# Note: an AIC/BIC grid-search route (compute_ic + select_lambda) was removed.
# The deprecated source lives at code/mechrt/deprecated/aic_lambda_select.R
# (see also REMOVED_OPTIONS.md section A).
# ==============================================================================

gi_select_lambda_cv <- function(
    design,
    lam_grid,
    tail = c("linear", "constant", "none"),
    nfold = 5L,
    cv_select_rule = c("1se", "min"),
    error_measure = c("deviance", "mse", "mae"),  # CV fold-loss
    verbose = FALSE,
    fixed_dow = NULL,             # hold DoW effects fixed across folds
    ...                           # passed to gi_solve
) {
  tail <- match.arg(tail)
  cv_select_rule <- match.arg(cv_select_rule)
  error_measure <- match.arg(error_measure)
  err_fun <- switch(error_measure,
    mse      = function(y, m) (y - m)^2,
    mae      = function(y, m) abs(y - m),
    deviance = function(y, m) {
      # Poisson deviance (per-obs): 2 * (y log(y/m) - (y - m)); limit at y=0 is 2*m
      term <- ifelse(y > 0, y * log(y / pmax(m, 1e-12)) - (y - m), m)
      2 * term
    }
  )

  n <- nrow(design$Z_full)
  if (n < 2L) stop("Need at least 2 observations for CV lambda tuning.")

  nfold <- as.integer(nfold)
  if (!is.finite(nfold) || nfold < 2L) nfold <- 2L
  nfold <- min(nfold, n)

  # Evenly spaced holdouts: fold k contains rows k, k+nfold, k+2*nfold, ...
  fold_id <- ((seq_len(n) - 1L) %% nfold) + 1L

  fold_err <- matrix(NA_real_, nrow = nfold, ncol = length(lam_grid))

  for (i in seq_along(lam_grid)) {
    lam <- lam_grid[i]
    if (verbose) cat(sprintf("  [CV %d/%d] lam = %.3e\n", i, length(lam_grid), lam))

    for (k in seq_len(nfold)) {
      hold_idx <- which(fold_id == k)
      train_idx <- which(fold_id != k)
      if (length(train_idx) == 0L || length(hold_idx) == 0L) next

      d_tr <- design
      d_tr$Y_valid   <- design$Y_valid[train_idx]
      d_tr$off_valid <- design$off_valid[train_idx]
      d_tr$Z_full    <- design$Z_full[train_idx, , drop = FALSE]
      if (isTRUE(design$use_dow))
        d_tr$dow_valid <- design$dow_valid[train_idx]

      fit_k <- tryCatch(
        gi_solve(d_tr, lam, tail = tail,
                    fixed_dow = fixed_dow, ...),
        error = function(e) NULL
      )
      if (is.null(fit_k) || !all(is.finite(fit_k$theta))) next

      # Holdout prediction (apply training DoW effects to holdout rows)
      Z_ho  <- design$Z_full[hold_idx, , drop = FALSE]
      off_ho <- design$off_valid[hold_idx]
      if (!is.null(fit_k$dow_effects)) {
        d_ho   <- fit_k$dow_effects[design$dow_valid[hold_idx]]
        mu_hold <- pmax(d_ho * (as.numeric(Z_ho %*% fit_k$theta) + off_ho), 1e-10)
      } else {
        mu_hold <- pmax(as.numeric(Z_ho %*% fit_k$theta) + off_ho, 1e-10)
      }
      y_hold <- design$Y_valid[hold_idx]
      fold_err[k, i] <- mean(err_fun(y_hold, mu_hold))
    }
  }

  cv_mean <- apply(fold_err, 2L, function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) Inf else mean(x)
  })

  cv_sd <- apply(fold_err, 2L, function(x) {
    x <- x[is.finite(x)]
    if (length(x) <= 1L) NA_real_ else sd(x)
  })

  cv_n <- apply(fold_err, 2L, function(x) sum(is.finite(x)))
  cv_se <- ifelse(cv_n > 0L, cv_sd / sqrt(cv_n), NA_real_)

  best_i_min <- which.min(cv_mean)
  if (!is.finite(cv_mean[best_i_min])) {
    cat("!!! CV lambda selection failed for all values\n")
    return(list(
      cv_mean = cv_mean,
      cv_sd = cv_sd,
      cv_se = cv_se,
      fold_errors = fold_err,
      best_i_min = NA_integer_,
      best_lam_min = NA_real_,
      best_i_1se = NA_integer_,
      best_lam_1se = NA_real_,
      cv_threshold_1se = NA_real_,
      best_lam = NA_real_,
      best_i = NA_integer_,
      best_fit = NULL
    ))
  }

  cv_threshold_1se <- cv_mean[best_i_min] + ifelse(is.finite(cv_se[best_i_min]), cv_se[best_i_min], 0.0)
  eligible_1se <- which(is.finite(cv_mean) & (cv_mean <= cv_threshold_1se))
  best_i_1se <- if (length(eligible_1se) > 0L) max(eligible_1se) else best_i_min

  # 2se rule (more conservative): largest lambda whose cv_mean is within 2 SE of min
  cv_threshold_2se <- cv_mean[best_i_min] + ifelse(is.finite(cv_se[best_i_min]), 2.0 * cv_se[best_i_min], 0.0)
  eligible_2se <- which(is.finite(cv_mean) & (cv_mean <= cv_threshold_2se))
  best_i_2se <- if (length(eligible_2se) > 0L) max(eligible_2se) else best_i_min

  best_i <- if (cv_select_rule == "1se") best_i_1se else best_i_min

  cat(sprintf("  Best lambda (CV): %.4e  (MSE = %.4f)\n", lam_grid[best_i], cv_mean[best_i]))

  best_fit <- gi_solve(design, lam_grid[best_i], tail = tail,
                          fixed_dow = fixed_dow, ...)

  list(
    cv_mean = cv_mean,
    cv_sd = cv_sd,
    cv_se = cv_se,
    fold_errors = fold_err,
    best_i_min = best_i_min,
    best_lam_min = lam_grid[best_i_min],
    best_i_1se = best_i_1se,
    best_lam_1se = lam_grid[best_i_1se],
    cv_threshold_1se = cv_threshold_1se,
    best_i_2se = best_i_2se,
    best_lam_2se = lam_grid[best_i_2se],
    cv_threshold_2se = cv_threshold_2se,
    best_lam = lam_grid[best_i],
    best_i = best_i,
    best_fit = best_fit
  )
}


# ==============================================================================
# 9.  LAMBDA SELECTION WITH DAY-OF-WEEK EFFECTS  (K-fold CV)
#
# Same CV structure as gi_select_lambda_cv, but each training fold estimates
# DoW effects jointly with theta via alternating descent.  Holdout predictions
# use the training-fold DoW scaling.
#
# Returns a list compatible with gi_select_lambda_cv's return format
# (best_lam_min, best_lam_1se, etc.).
# ==============================================================================

gi_select_lambda_cv_dow <- function(
    design,
    dates_valid,           # Date vector, length = nrow(design$Z_full)
    lam_grid,
    tail           = c("linear", "constant", "none"),
    nfold          = 5L,
    cv_select_rule = c("1se", "min"),
    error_measure  = c("deviance", "mse", "mae"),   # CV fold-loss
    max_outer      = 15L,
    dow_tol        = 1e-3,
    verbose        = FALSE,
    ...
) {
  tail           <- match.arg(tail)
  error_measure  <- match.arg(error_measure)
  err_fun        <- switch(error_measure,
    mse      = function(y, m) (y - m)^2,
    mae      = function(y, m) abs(y - m),
    deviance = function(y, m) 2 * ifelse(y > 0,
                 y * log(y / pmax(m, 1e-12)) - (y - m),
                 m)
  )
  cv_select_rule <- match.arg(cv_select_rule)

  n       <- nrow(design$Z_full)
  dow_idx <- as.integer(format(as.Date(dates_valid), "%u"))
  nfold   <- min(as.integer(nfold), n)
  fold_id <- ((seq_len(n) - 1L) %% nfold) + 1L

  fold_err <- matrix(NA_real_, nrow = nfold, ncol = length(lam_grid))

  for (i in seq_along(lam_grid)) {
    lam <- lam_grid[i]
    if (verbose) cat(sprintf("  [DoW CV %d/%d] lam = %.3e\n", i, length(lam_grid), lam))

    for (k in seq_len(nfold)) {
      train_idx <- which(fold_id != k)
      hold_idx  <- which(fold_id == k)
      if (length(train_idx) == 0L || length(hold_idx) == 0L) next

      # Training sub-design.  IMPORTANT: this outer loop manages DoW manually
      # via the alternating log_d coordinate descent below, so we explicitly
      # disable DoW inside gi_solve.  Otherwise gi_solve would read the
      # full-length design$dow_valid against the trimmed (|train_idx|-row)
      # Z_full and trigger recycling warnings / wrong fits.
      d_tr           <- design
      d_tr$Y_valid   <- design$Y_valid[train_idx]
      d_tr$off_valid <- design$off_valid[train_idx]
      d_tr$Z_full    <- design$Z_full[train_idx, , drop = FALSE]
      d_tr$use_dow   <- FALSE
      d_tr$dow_valid <- NULL
      dow_tr         <- dow_idx[train_idx]

      # Alternating DoW + theta on training set
      log_d <- rep(0.0, 7L)
      fit_k <- NULL
      for (outer in seq_len(max_outer)) {
        d_t_tr            <- exp(log_d[dow_tr])
        d_tr_sc           <- d_tr
        d_tr_sc$Z_full    <- d_tr$Z_full    * d_t_tr
        d_tr_sc$off_valid <- d_tr$off_valid * d_t_tr

        fit_k <- tryCatch(
          gi_solve(d_tr_sc, lam, tail = tail, ...),
          error = function(e) NULL
        )
        if (is.null(fit_k) || !all(is.finite(fit_k$theta))) { fit_k <- NULL; break }

        m_tr <- pmax(as.numeric(d_tr$Z_full %*% fit_k$theta) + d_tr$off_valid, 1e-10)
        log_d_new <- vapply(seq_len(7L), function(j) {
          idx <- which(dow_tr == j)
          if (length(idx) == 0L) return(0.0)
          log(sum(d_tr$Y_valid[idx]) / sum(m_tr[idx]))
        }, numeric(1L))
        log_d_new <- log_d_new - mean(log_d_new)

        delta <- max(abs(log_d_new - log_d))
        log_d <- log_d_new
        if (outer > 2L && delta < dow_tol) break
      }
      if (is.null(fit_k)) next

      # Holdout prediction using training DoW effects
      d_t_ho <- exp(log_d[dow_idx[hold_idx]])
      Z_ho   <- design$Z_full[hold_idx, , drop = FALSE]
      off_ho <- design$off_valid[hold_idx]
      mu_ho  <- pmax(d_t_ho * (as.numeric(Z_ho %*% fit_k$theta) + off_ho), 1e-10)
      y_ho   <- design$Y_valid[hold_idx]
      fold_err[k, i] <- mean(err_fun(y_ho, mu_ho))
    }
  }

  # Lambda selection  (identical logic to gi_select_lambda_cv)
  cv_mean <- apply(fold_err, 2L, function(x) {
    x <- x[is.finite(x)]; if (!length(x)) Inf else mean(x)
  })
  cv_sd <- apply(fold_err, 2L, function(x) {
    x <- x[is.finite(x)]; if (length(x) <= 1L) NA_real_ else sd(x)
  })
  cv_n  <- apply(fold_err, 2L, function(x) sum(is.finite(x)))
  cv_se <- ifelse(cv_n > 0L, cv_sd / sqrt(cv_n), NA_real_)

  best_i_min <- which.min(cv_mean)
  if (!is.finite(cv_mean[best_i_min])) {
    cat("!!! DoW CV lambda selection failed for all values\n")
    return(list(
      cv_mean = cv_mean, cv_sd = cv_sd, cv_se = cv_se,
      fold_errors = fold_err,
      best_i_min = NA_integer_, best_lam_min = NA_real_,
      best_i_1se = NA_integer_, best_lam_1se = NA_real_,
      cv_threshold_1se = NA_real_,
      best_i_2se = NA_integer_, best_lam_2se = NA_real_,
      cv_threshold_2se = NA_real_,
      best_lam = NA_real_, best_i = NA_integer_
    ))
  }

  cv_threshold_1se <- cv_mean[best_i_min] +
    ifelse(is.finite(cv_se[best_i_min]), cv_se[best_i_min], 0.0)
  eligible_1se <- which(is.finite(cv_mean) & cv_mean <= cv_threshold_1se)
  best_i_1se   <- if (length(eligible_1se) > 0L) max(eligible_1se) else best_i_min

  # 2se rule (more conservative): largest lambda whose cv_mean is within 2 SE of min
  cv_threshold_2se <- cv_mean[best_i_min] +
    ifelse(is.finite(cv_se[best_i_min]), 2.0 * cv_se[best_i_min], 0.0)
  eligible_2se <- which(is.finite(cv_mean) & cv_mean <= cv_threshold_2se)
  best_i_2se   <- if (length(eligible_2se) > 0L) max(eligible_2se) else best_i_min

  best_i       <- if (cv_select_rule == "1se") best_i_1se else best_i_min

  cat(sprintf("  [DoW CV] Best lambda: %.4e  (CV-MSE = %.4f)\n",
              lam_grid[best_i], cv_mean[best_i]))

  list(
    method           = "cv_dow",
    cv_mean          = cv_mean, cv_sd = cv_sd, cv_se = cv_se,
    fold_errors      = fold_err,
    best_i_min       = best_i_min, best_lam_min = lam_grid[best_i_min],
    best_i_1se       = best_i_1se, best_lam_1se = lam_grid[best_i_1se],
    cv_threshold_1se = cv_threshold_1se,
    best_i_2se       = best_i_2se, best_lam_2se = lam_grid[best_i_2se],
    cv_threshold_2se = cv_threshold_2se,
    best_lam         = lam_grid[best_i], best_i = best_i
  )
}

