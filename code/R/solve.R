# MechRt :: fit.R
#
# gi_solve — penalised IRLS + KKT solver.  This is the standard MechRt
# fitter; use it for retrospective fits AND for real-time fits with no
# tapered penalty.
#
# When to use which solver:
#   * gi_solve              KKT solve, supports DoW profiling.  Used for
#                               retrospective fits and tapered-OFF real-time
#                               fits.
#   * .gi_solve_taper       Null-space IRLS, supports tapered penalty.
#                               Use whenever lam_taper > 0; the KKT path in
#                               gi_solve is ill-conditioned for large
#                               tapers.  No DoW profiling — pre-apply DoW to
#                               the design first.  Defined in taper.R.
#
# ==============================================================================
# FIT BETA VIA PENALISED IRLS + KKT
#
# Solves:  max  sum_t [y_t log(mu_t) - mu_t]  -  lam/2 * theta' P theta
#          s.t. A theta = 0
#
# where mu_t = Z_t' theta + off_t  (identity-link Poisson).
#
# Each IRLS step is a weighted least-squares problem:
#   min  (y - off - Z theta)' W (y - off - Z theta)  +  lam * theta' P theta
#   s.t. A theta = 0
# where W = diag(1/mu).  Solved via the KKT system:
#   [Z'WZ + lam P   A'] [theta]   [Z'W(y - off)]
#   [A               0] [nu  ] = [0            ]
# ==============================================================================

gi_solve <- function(
    design,
    lam,
    tail            = c("linear", "constant", "none"),
    max_iter        = 300L,
    tol             = 1e-8,
    verbose         = FALSE,
    fixed_dow       = NULL,   # named length-7 vector of DoW effects to hold fixed (skip profiling)
    theta_init      = NULL    # optional warm-start for theta (skips fresh init)
) {
  # Deprecated args, hardcoded to canonical defaults.  See
  # code/mechrt/deprecated/REMOVED_OPTIONS.md (sections B, D) for the
  # ridge-penalty and constrain_mu QP code paths these dead-but-recoverable
  # branches preserve.
  penalty      <- "pspline"  # ridge dead behind this hardcode (section D)
  constrain_mu <- FALSE      # QP path dead behind this hardcode  (section B)

  tail    <- match.arg(tail)

  if (identical(design$basis, "natural") && tail == "constant") {
    stop("MechRt: basis = 'natural' with tail = 'constant' is deprecated. ",
         "Use basis = 'regression' with tail = 'constant' for an explicit ",
         "constant-tail constraint, or stick with the smoothing spline ",
         "(basis = 'natural') and leave the tail unconstrained.")
  }

  Z_orig <- design$Z_full
  y      <- design$Y_valid
  off_orig <- design$off_valid
  np     <- design$np
  nb     <- design$n_basis

  # Normalization factors: divide loss by n_obs, penalty by n_penalty
  n_obs     <- if (!is.null(design$n_obs))     design$n_obs     else length(y)
  n_penalty <- if (!is.null(design$n_penalty)) design$n_penalty else max(nb - 2L, 1L)

  P <- if (penalty == "pspline") design$P_pspline else design$P_ridge
  A <- design$tail_constraints[[tail]]
  n_eq <- nrow(A)
  b_kkt <- rep(0.0, n_eq)  # all derivative constraints are homogeneous

  # --- Day-of-week profiling setup ------------------------------------------
  # When design$use_dow is TRUE, the model is:
  #   Y_t ~ Pois( d_{dow(t)} * mu_t ),   mu_t = Z_t' theta + off_t
  #   d_j >= 0,  prod(d_j) = 1  (7 multiplicative effects, geometric mean = 1)
  #
  # At each IRLS step, d_j is updated in closed form (profile likelihood):
  #   d_j = sum(Y | dow=j) / sum(mu | dow=j),  then normalise to prod = 1.
  # The IRLS step then uses the scaled design: Z_scaled = diag(d_t) Z,
  # off_scaled = d_t * off.
  use_dow <- isTRUE(design$use_dow)
  dow_idx <- design$dow_valid   # integer 1..7 per valid row (NULL if no DoW)

  # If fixed_dow is supplied, use those effects and skip profiling in IRLS
  dow_is_fixed <- FALSE
  if (!is.null(fixed_dow) && use_dow) {
    stopifnot(length(fixed_dow) == 7L)
    d_vec <- as.numeric(fixed_dow)
    dow_is_fixed <- TRUE
  } else {
    d_vec <- rep(1.0, 7L)       # initialise: all DoW effects = 1
  }

  # --- Initialise theta to give mu ≈ mean(y) (or accept caller's warm start) -
  if (!is.null(theta_init) && length(theta_init) == np &&
      all(is.finite(theta_init))) {
    theta <- as.numeric(theta_init)
  } else {
    rs <- rowSums(Z_orig[, seq_len(nb), drop = FALSE])
    rs_pos <- rs[rs > 0]
    ref <- if (length(rs_pos) > 0) median(rs_pos) else 1.0
    init_val <- max(1e-6, (mean(y) - mean(off_orig)) / ref)
    theta <- c(rep(init_val, nb), rep(0.0, np - nb))
  }

  # Pre-compute projection onto constraint surface.
  # Homogeneous (A theta = 0): proj(v) = v - A'(AA')^{-1} A v
  # Inhomogeneous (A theta = b): proj(v) = v - A'(AA')^{-1}(A v - b)
  AAt_inv <- if (n_eq > 0)
    tryCatch(solve(tcrossprod(A)), error = function(e) NULL)
  else NULL

  project_constraint <- function(v) {
    if (is.null(AAt_inv)) return(v)
    v - as.numeric(t(A) %*% AAt_inv %*% (A %*% v - b_kkt))
  }

  theta <- project_constraint(theta)   # initial projection

  # --- IRLS iterations -------------------------------------------------------
  llik_prev <- -Inf
  for (iter in seq_len(max_iter)) {

    # -- DoW profile step: update d_j from current spline-based mu -----------
    if (use_dow) {
      if (!dow_is_fixed) {
        mu_underlying <- pmax(as.numeric(Z_orig %*% theta) + off_orig, 1e-10)
        for (j in seq_len(7L)) {
          idx_j <- which(dow_idx == j)
          if (length(idx_j) > 0L) {
            d_vec[j] <- sum(y[idx_j]) / sum(mu_underlying[idx_j])
          }
        }
        # Normalise: geometric mean = 1  =>  d_j <- d_j / prod(d_j)^(1/7)
        d_vec <- pmax(d_vec, 1e-10)
        d_vec <- d_vec / exp(mean(log(d_vec)))
      }

      # Scale design by current DoW effects
      d_t <- d_vec[dow_idx]
      Z   <- Z_orig * d_t
      off <- off_orig * d_t
    } else {
      Z   <- Z_orig
      off <- off_orig
    }

    # -- IRLS step on (possibly DoW-scaled) design ---------------------------
    mu <- as.numeric(Z %*% theta) + off
    mu <- pmax(mu, 1e-10)          # guard against log(0)

    w     <- 1.0 / mu              # Poisson working weights (identity link)
    ZtWZ  <- crossprod(Z * sqrt(w))                    # Z' diag(w) Z
    ZtWyo <- as.numeric(crossprod(Z, w * (y - off)))   # Z' diag(w) (y - off)

    # Normalize loss by n_obs and penalty by n_penalty for comparability
    ZtWZ_n  <- ZtWZ  / n_obs
    ZtWyo_n <- ZtWyo / n_obs
    P_n     <- P     / n_penalty

    # Regularisation: tiny ridge added for numerical stability
    eps_ridge <- max(1e-10 * max(diag(ZtWZ_n), 1.0), .Machine$double.eps)
    M <- ZtWZ_n + lam * P_n + eps_ridge * diag(np)
    # NS "constant head" augmented Lagrangian penalty (flat first few days)
    if (!is.null(design$head_const_extra_M) && isTRUE(design$head_const_days > 0L)) {
      M <- M + design$head_const_extra_M
    }

    # KKT linear solve.
    if (n_eq == 0L) {
      theta_new <- tryCatch(
        solve(M, ZtWyo_n),
        error = function(e) NULL
      )
    } else {
      KKT <- rbind(cbind(M,                          t(A)),
                   cbind(A, matrix(0.0, n_eq, n_eq)))
      rhs <- c(ZtWyo_n, b_kkt)
      sol <- tryCatch(solve(KKT, rhs), error = function(e) NULL)
      theta_new <- if (!is.null(sol)) sol[seq_len(np)] else NULL
    }

    # Fallback: ridge-only step if KKT failed; re-project to restore constraint
    if (is.null(theta_new) || !all(is.finite(theta_new))) {
      lam_fb <- max(lam, 1e-4)
      M_fb   <- ZtWZ_n + lam_fb * (design$P_ridge / n_penalty) + eps_ridge * diag(np)
      theta_new <- tryCatch(solve(M_fb, ZtWyo_n), error = function(e) theta)
    }
    # Always project theta_new onto the constraint surface
    theta_new <- project_constraint(theta_new)

    # Backtracking line search: keep mu > 0
    step <- theta_new - theta
    ss   <- 1.0
    for (k in seq_len(30L)) {
      mu_try <- as.numeric(Z %*% (theta + ss * step)) + off
      if (all(mu_try > 0)) break
      ss <- ss * 0.5
    }
    theta <- theta + ss * step

    # Convergence check (use DoW-scaled mu for log-likelihood)
    mu_new  <- pmax(as.numeric(Z %*% theta) + off, 1e-10)
    llik    <- sum(y * log(mu_new) - mu_new)
    delta   <- abs(llik - llik_prev) / (1.0 + abs(llik_prev))
    if (verbose) cat(sprintf("  iter %3d: llik = %10.4f  delta = %.2e\n",
                              iter, llik, delta))
    if (iter > 5L && delta < tol) break
    llik_prev <- llik
  }

  # Final projection: guarantees A theta = 0 to machine precision
  theta <- project_constraint(theta)

  # Final mu (with DoW scaling if applicable)
  if (use_dow) {
    d_t <- d_vec[dow_idx]
    Z   <- Z_orig * d_t
    off <- off_orig * d_t
  } else {
    Z   <- Z_orig
    off <- off_orig
  }
  mu_final <- pmax(as.numeric(Z %*% theta) + off, 1e-10)

  llik_final <- sum(y * log(mu_final) - mu_final)

  if (verbose) cat(sprintf("  Converged in %d iterations.  llik = %.4f\n",
                            iter, llik_final))

  # DoW effects: named vector (Mon..Sun), multiplicative, product = 1
  dow_names   <- c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
  dow_effects <- if (use_dow) setNames(d_vec, dow_names) else NULL

  if (use_dow && verbose) {
    cat("  DoW effects (multiplicative, product=1):\n")
    cat("    ", paste(sprintf("%s=%.3f", dow_names, d_vec), collapse = "  "), "\n")
  }

  list(
    theta           = theta,
    mu              = mu_final,
    llik            = llik_final,
    n_iter          = iter,
    penalty         = penalty,
    tail_constraint = tail,
    dow_effects     = dow_effects   # NULL when DoW not used
  )
}


