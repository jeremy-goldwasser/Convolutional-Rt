# MechRt :: rt_extract.R
#
# Convert a fitted theta into beta(t) and R_t(t), with pointwise Laplace CIs
# and (optionally) a sup-norm simultaneous confidence band drawn from the
# Laplace covariance.
# ==============================================================================
# beta(t) AND Rt(t) WITH LAPLACE CONFIDENCE INTERVALS
#
# Variance: Var(theta_hat) ≈ M_pen^{-1} Z'WZ M_pen^{-1}  (sandwich)
#           where M_pen = Z'WZ + lam*P
# For beta(t) = G(t)' theta:  SE^2(beta_hat(t)) = G(t)' Var(theta) G(t)
#
# Quasi-Poisson overdispersion correction (overdispersion = TRUE):
#   Under Poisson, Var(Y_t) = mu_t.  If the data are overdispersed the
#   true variance is phi * mu_t with phi > 1.  The point estimates are
#   unchanged (Poisson MLE is consistent regardless of phi), but the
#   sandwich covariance of theta_hat scales by phi:
#
#       V_quasi(theta) = phi * V_pois(theta)
#       SE_quasi       = sqrt(phi) * SE_pois
#
#   phi is estimated by the Pearson statistic:
#       phi_hat = sum( (y - mu)^2 / mu ) / (n - df_eff)
#   where df_eff = tr(H_pen^{-1} Z'WZ) is the effective degrees of freedom.
# ==============================================================================

#' Estimate the quasi-Poisson dispersion parameter phi from a fit.
#'
#' @param design  Design list (needs Z_full, Y_valid, off_valid, np, P_pspline/P_ridge).
#' @param fit     Fit list from solve_mechrt (needs theta, mu, penalty, dow_effects).
#' @param lam     Penalty weight used for the fit.
#' @return Scalar phi_hat (>= 1; clamped below at 1.0).
gi_estimate_dispersion <- function(design, fit, lam) {
  Z   <- design$Z_full
  y   <- design$Y_valid
  off <- design$off_valid
  np  <- design$np

  if (!is.null(fit$dow_effects) && isTRUE(design$use_dow)) {
    d_t <- fit$dow_effects[design$dow_valid]
    Z   <- Z * d_t
    off <- off * d_t
  }

  theta <- fit$theta
  mu    <- pmax(as.numeric(Z %*% theta) + off, 1e-10)

  # Pearson chi-squared
  pearson_chi2 <- sum((y - mu)^2 / mu)

  # Effective df (hat-matrix trace)
  P     <- if (fit$penalty == "pspline") design$P_pspline else design$P_ridge
  # Convert normalized lambda to effective lambda for unnormalized computation
  n_obs_d     <- if (!is.null(design$n_obs))     design$n_obs     else length(y)
  n_penalty_d <- if (!is.null(design$n_penalty)) design$n_penalty else max(design$n_basis - 2L, 1L)
  lam_eff <- lam * n_obs_d / n_penalty_d
  w     <- 1.0 / mu
  ZtWZ  <- crossprod(Z * sqrt(w))
  eps   <- max(1e-10 * max(diag(ZtWZ), 1.0), .Machine$double.eps)
  M_pen <- ZtWZ + lam_eff * P + eps * diag(np)
  # NS "constant tail" augmented Lagrangian penalty (NULL/no-op for RS)
  if (!is.null(design$head_const_extra_M) && isTRUE(design$head_const_days > 0L)) {
    M_pen <- M_pen + design$head_const_extra_M
  }
  Hi    <- tryCatch(solve(M_pen), error = function(e) diag(np) * 1e-6)
  df_eff <- sum(diag(Hi %*% ZtWZ))

  n <- length(y)
  phi <- pearson_chi2 / (n - df_eff)

  # Clamp: phi < 1 would mean underdispersion; don't shrink CIs below Poisson
  max(phi, 1.0)
}

# girt's R_t IS the spline; no S(t)/N or 1/gamma rescaling.
gi_extract_rt <- function(
    fit,
    design,
    lam,
    lam_taper = 0,  # tapered penalty weight (0 = no tapering)
    P_taper   = NULL,
    level = 0.95,
    overdispersion = TRUE
) {
  theta <- fit$theta
  if (is.null(theta) || !all(is.finite(theta))) {
    empty <- data.frame(day = integer(),
                        Rt_mean = numeric(), Rt_lo = numeric(), Rt_hi = numeric())
    return(empty)
  }

  np       <- design$np
  fr       <- design$fr
  max_time <- design$max_time
  B        <- design$B
  kink_vec <- design$kink_vec
  use_kink <- design$use_kink
  Z        <- design$Z_full
  off      <- design$off_valid

  # If DoW effects were estimated, apply them to Z and off for variance calc
  if (!is.null(fit$dow_effects) && isTRUE(design$use_dow)) {
    d_t <- fit$dow_effects[design$dow_valid]
    Z   <- Z * d_t
    off <- off * d_t
  }

  # Full basis matrix G: (max_time+1) x np
  G <- if (use_kink) cbind(B, kink_vec) else B

  # beta(t) = G(t)' theta
  beta_t <- as.numeric(G %*% theta)

  # Mask outside spline domain
  t_vals    <- 0:max_time
  in_domain <- (t_vals >= fr[1]) & (t_vals <= fr[2])
  beta_t[!in_domain] <- NA_real_

  # --- Sandwich variance of theta_hat --------------------------------------
  # Convert normalized lambda to effective lambda for unnormalized CI computation
  n_obs_rt     <- if (!is.null(design$n_obs))     design$n_obs     else nrow(Z)
  n_penalty_rt <- if (!is.null(design$n_penalty)) design$n_penalty else max(design$n_basis - 2L, 1L)
  lam_eff <- lam * n_obs_rt / n_penalty_rt
  lam_taper_eff <- lam_taper * n_obs_rt / n_penalty_rt

  P_base <- if (fit$penalty == "pspline") design$P_pspline else design$P_ridge
  has_taper <- (lam_taper > 0) && !is.null(P_taper)
  P_total <- lam_eff * P_base
  if (has_taper) P_total <- P_total + lam_taper_eff * P_taper

  if (!is.null(design$head_const_extra_M) && isTRUE(design$head_const_days > 0L)) {
    P_total <- P_total + design$head_const_extra_M
  }

  mu   <- pmax(as.numeric(Z %*% theta) + off, 1e-10)
  w    <- 1.0 / mu
  ZtWZ <- crossprod(Z * sqrt(w))

  if (has_taper) {
    # Null-space method: stable when lam_taper is large (avoids ill-conditioned full solve)
    tail_key <- if (!is.null(fit$tail_constraint)) fit$tail_constraint else "linear"
    A    <- design$tail_constraints[[tail_key]]
    n_eq <- nrow(A)
    if (n_eq > 0) {
      qr_A   <- qr(t(A))
      N_null <- qr.Q(qr_A, complete = TRUE)[, (n_eq + 1L):np, drop = FALSE]
    } else {
      N_null <- diag(np)
    }
    n_free    <- ncol(N_null)
    ZtWZ_N    <- crossprod(Z %*% N_null * sqrt(w))
    P_total_N <- crossprod(N_null, P_total %*% N_null)
    eps_N     <- max(1e-10 * max(diag(ZtWZ_N), 1.0), .Machine$double.eps)
    M_pen_N   <- ZtWZ_N + P_total_N + eps_N * diag(n_free)
    Hi_N      <- tryCatch(solve(M_pen_N), error = function(e) diag(n_free) * 1e-6)
    V_theta   <- N_null %*% (Hi_N %*% ZtWZ_N %*% Hi_N) %*% t(N_null)
  } else {
    eps     <- max(1e-10 * max(diag(ZtWZ), 1.0), .Machine$double.eps)
    M_pen   <- ZtWZ + P_total + eps * diag(np)
    Hi      <- tryCatch(solve(M_pen), error = function(e) diag(np) * 1e-6)
    V_theta <- Hi %*% ZtWZ %*% Hi
  }

  # SE of beta(t): sqrt(diag(G V_theta G'))  computed row-wise
  GV      <- G %*% V_theta
  se_beta <- sqrt(pmax(0.0, rowSums(GV * G)))
  se_beta[!in_domain] <- NA_real_

  # Quasi-Poisson overdispersion correction: inflate SE by sqrt(phi_hat)
  if (isTRUE(overdispersion)) {
    phi_hat <- gi_estimate_dispersion(design, fit, lam)
    se_beta <- se_beta * sqrt(phi_hat)
  }

  # girt's R_t IS the spline; no rescaling.
  scale <- rep(1.0, max_time + 1L)

  z_crit <- qnorm(1.0 - (1.0 - level) / 2.0)

  data.frame(
    day      = t_vals,
    Rt_mean  = beta_t * scale,
    Rt_lo    = (beta_t - z_crit * se_beta) * scale,
    Rt_hi    = (beta_t + z_crit * se_beta) * scale
  ) |> filter(!is.na(Rt_mean), .data$day >= fr[1])
}


# ==============================================================================
# extract_jumps: point estimates + Wald CIs for the Heaviside jump coefficients
# added by build_design(jump_times = ...) in the natural-spline branch.
#
# The jump coefficients are the LAST n_jump entries of fit$theta. Their
# variances are diag(V_theta)[(n_basis+1):np] from the same sandwich V_theta
# used in gi_extract_rt. Returns NULL if the design has no jump columns.
# ==============================================================================
extract_jumps <- function(
    fit,
    design,
    lam,
    lam_taper = 0,
    P_taper   = NULL,
    level     = 0.95,
    overdispersion = TRUE
) {
  n_jump <- if (!is.null(design$n_jump)) design$n_jump else 0L
  if (n_jump <= 0L) return(NULL)

  theta <- fit$theta
  if (is.null(theta) || !all(is.finite(theta))) return(NULL)

  np      <- design$np
  n_basis <- design$n_basis
  Z       <- design$Z_full
  off     <- design$off_valid

  if (!is.null(fit$dow_effects) && isTRUE(design$use_dow)) {
    d_t <- fit$dow_effects[design$dow_valid]
    Z   <- Z * d_t
    off <- off * d_t
  }

  n_obs_rt     <- if (!is.null(design$n_obs))     design$n_obs     else nrow(Z)
  n_penalty_rt <- if (!is.null(design$n_penalty)) design$n_penalty else max(n_basis - 2L, 1L)
  lam_eff       <- lam       * n_obs_rt / n_penalty_rt
  lam_taper_eff <- lam_taper * n_obs_rt / n_penalty_rt

  P_base <- if (fit$penalty == "pspline") design$P_pspline else design$P_ridge
  has_taper <- (lam_taper > 0) && !is.null(P_taper)
  P_total <- lam_eff * P_base
  if (has_taper) P_total <- P_total + lam_taper_eff * P_taper
  if (!is.null(design$head_const_extra_M) && isTRUE(design$head_const_days > 0L)) {
    P_total <- P_total + design$head_const_extra_M
  }

  mu   <- pmax(as.numeric(Z %*% theta) + off, 1e-10)
  w    <- 1.0 / mu
  ZtWZ <- crossprod(Z * sqrt(w))
  eps  <- max(1e-10 * max(diag(ZtWZ), 1.0), .Machine$double.eps)
  M_pen   <- ZtWZ + P_total + eps * diag(np)
  Hi      <- tryCatch(solve(M_pen), error = function(e) diag(np) * 1e-6)
  V_theta <- Hi %*% ZtWZ %*% Hi

  jump_idx_in_theta <- (n_basis + 1L):np
  beta_jump <- as.numeric(theta[jump_idx_in_theta])
  var_jump  <- diag(V_theta)[jump_idx_in_theta]
  se_jump   <- sqrt(pmax(0.0, var_jump))

  if (isTRUE(overdispersion)) {
    phi_hat <- gi_estimate_dispersion(design, fit, lam)
    se_jump <- se_jump * sqrt(phi_hat)
  }

  z_crit <- qnorm(1.0 - (1.0 - level) / 2.0)
  data.frame(
    jump_idx = seq_len(n_jump),
    day      = as.integer(design$jump_times_idx),
    beta     = beta_jump,
    se       = se_jump,
    lo       = beta_jump - z_crit * se_jump,
    hi       = beta_jump + z_crit * se_jump
  )
}


# ==============================================================================
# 5b.  SIMULTANEOUS CONFIDENCE BAND  (simulation-based sup-norm)
#
# 1. Draw theta_m ~ N(theta_hat, V_theta),  m = 1 .. n_sim
# 2. For each draw compute sup_t |Rt_m(t) - Rt_hat(t)| / SE_Rt(t)
# 3. c_alpha = quantile(..., level)
# 4. Simultaneous band: Rt_hat(t) ± c_alpha * SE_Rt(t)
#
# Returns additional columns sim_lo and sim_hi appended to the pointwise df.
# ==============================================================================

gi_extract_rt_simband <- function(fit, design, lam, level = 0.95, n_sim = 5000L,
                                   overdispersion = TRUE) {
  pw <- gi_extract_rt(fit, design, lam = lam,
                      level = level, overdispersion = overdispersion)
  if (nrow(pw) == 0) return(pw)

  nb       <- design$n_basis
  np       <- design$np
  fr       <- design$fr
  max_time <- design$max_time
  B        <- design$B
  kink_vec <- design$kink_vec
  use_kink <- design$use_kink
  Z        <- design$Z_full
  off      <- design$off_valid
  theta    <- fit$theta

  G <- if (use_kink) cbind(B, kink_vec) else B
  t_vals    <- 0:max_time
  in_domain <- (t_vals >= fr[1]) & (t_vals <= fr[2])

  P   <- if (fit$penalty == "pspline") design$P_pspline else design$P_ridge
  # Convert normalized lambda to effective lambda for unnormalized computation
  n_obs_sim     <- if (!is.null(design$n_obs))     design$n_obs     else nrow(Z)
  n_penalty_sim <- if (!is.null(design$n_penalty)) design$n_penalty else max(nb - 2L, 1L)
  lam_eff <- lam * n_obs_sim / n_penalty_sim
  mu  <- pmax(as.numeric(Z %*% theta) + off, 1e-10)
  w   <- 1.0 / mu
  ZtWZ <- crossprod(Z * sqrt(w))
  eps  <- max(1e-10 * max(diag(ZtWZ), 1.0), .Machine$double.eps)
  M_pen_sim <- ZtWZ + lam_eff * P + eps * diag(np)
  if (!is.null(design$head_const_extra_M) && isTRUE(design$head_const_days > 0L)) {
    M_pen_sim <- M_pen_sim + design$head_const_extra_M
  }
  Hi   <- tryCatch(solve(M_pen_sim),
                   error = function(e) diag(np) * 1e-6)
  V    <- Hi %*% ZtWZ %*% Hi

  # Scale covariance by phi if overdispersion correction is requested
  if (isTRUE(overdispersion)) {
    phi_hat <- gi_estimate_dispersion(design, fit, lam)
    V <- phi_hat * V
  }

  L <- tryCatch(
    t(chol(V)),
    error = function(e) {
      eig <- eigen(V, symmetric = TRUE)
      eig$vectors %*% diag(sqrt(pmax(eig$values, 0.0)))
    }
  )

  # girt's R_t IS the spline; no rescaling.
  scale <- rep(1.0, max_time + 1L)

  G_dom     <- G[in_domain, , drop = FALSE]
  scale_dom <- scale[in_domain]

  # SE_Rt at each in-domain time point
  GV_dom  <- G_dom %*% V
  se_rt   <- sqrt(pmax(0.0, rowSums(GV_dom * G_dom))) * scale_dom

  # Simulated standardised suprema
  Z_draws     <- matrix(rnorm(np * n_sim), np, n_sim)
  delta_beta  <- G_dom %*% (L %*% Z_draws)          # n_dom x n_sim
  delta_rt    <- delta_beta * scale_dom              # n_dom x n_sim
  se_rt_mat   <- matrix(se_rt, nrow = length(se_rt), ncol = n_sim)
  sup_std <- apply(abs(delta_rt) / pmax(se_rt_mat, 1e-12), 2, max)
  c_alpha <- quantile(sup_std, level, na.rm = TRUE)

  # Recover pointwise SE from the already-computed pw bounds
  z_crit   <- qnorm(1.0 - (1.0 - level) / 2.0)
  se_rt_pw <- (pw$Rt_mean - pw$Rt_lo) / z_crit   # pw$Rt_lo = Rt_mean - z*SE

  pw$sim_lo <- pw$Rt_mean - c_alpha * se_rt_pw
  pw$sim_hi <- pw$Rt_mean + c_alpha * se_rt_pw

  pw
}


