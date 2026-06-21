# girt :: design.R
#
# Build the design matrix Z and offset for the penalized Poisson regression that
# fits R_t through the renewal equation, predicting the observed counts.
#
# Model (identity-link Poisson, observation-level likelihood):
#
#   Y_t ~ Poisson( mu_t ),   mu_t = rho_t * omega_(t mod 7) * sum_{s<t} R_s * X_{t,s}
#   X_{t,s} = Lambda_s * pi_EY[t-s],   Lambda_s = sum_{k>=1} g_k * X_exposure[s-k]
#
# i.e. the renewal force Lambda (exposures convolved with the generation interval
# g) is generated at source day s, scaled by R_s, then pushed forward to the
# observation day t through the reporting delay pi_EY.  R_t = (B theta)_t is a
# natural cubic smoothing spline; severity rho and day-of-week omega are applied
# downstream (gi_apply_severity_to_design / gi_solve DoW profiling).
#
# This is the generation-interval analogue of MechRt's compartmental design:
# there is NO latent (E->I) or infectious-period (I->R) distribution and NO
# prevalence/recovery rate here -- only the generation interval g, the reporting
# delay pi_EY, and the exposures.  (With g = zeta^EI/mu^IR the resulting X_{t,s},
# hence the whole fit, is identical to MechRt -- Prop. compartmental-generation.)
#
# Natural cubic spline basis (linear tails) + integrated 2nd-derivative^2
# roughness penalty = a smoothing spline.  Returns the design list consumed by
# gi_solve / gi_select_lambda_cv / gi_extract_rt / gi_build_tapered_penalty.

library(splines)

build_gi_design <- function(
    E_star,                  # exposure incidence X (length max_time+1, days 0..max_time)
    g,                       # generation-interval pmf, g[k] = P(transmission k days after exposure)
    pi_EY,                   # exposure -> observed-outcome (reporting) delay pmf
    max_time,                # integer
    inf_inc,                 # observed counts Y (regression response), length max_time+1
    y_min_count    = 0,
    y_min_start    = NULL,
    y_min_end      = NULL,
    knot_step      = 5L,     # interior-knot spacing for the natural-spline basis
    spline_degree  = 3L,     # natural cubic spline (kept for API symmetry; always cubic)
    first_rt_index = NULL,   # first day (0-based) or Date where R_t is fit (burn-in anchor)
    dates          = NULL,   # Date vector length max_time+1 (to resolve Date first_rt_index / DoW)
    dow_dates      = NULL,   # Date vector length max_time+1 enabling multiplicative DoW effects
    head_const_days = 0L,    # force R_t flat over the first head_const_days (0 = off)
    jump_times     = NULL,   # 0-based day indices (or Dates) where R_t may jump (unpenalized Heaviside cols)
    pre_hist_rate  = 0,      # constant pre-day-0 exposure rate (survival-shaped offset); usually 0
    seed_rate      = 0       # pre-day-0 transmission seed: adds seed_rate * shifted survival of g
                            # to the renewal force (the GI-native form of MechRt's P_hat pre-history
                            # correction).  Pairs with pre_hist_rate (the offset).  Usually 0.
) {
  stopifnot(length(E_star) == max_time + 1L)
  if (abs(sum(g) - 1) > 1e-6)
    warning(sprintf("build_gi_design: g does not sum to 1 (sum = %.6f); it should be a normalized generation-interval pmf.", sum(g)))

  if (is.null(y_min_start) || length(y_min_start) == 0) y_min_start <- y_min_count
  if (is.null(y_min_end)   || length(y_min_end)   == 0) y_min_end   <- y_min_count
  knot_step <- as.integer(knot_step)
  if (knot_step < 1L) stop("knot_step must be >= 1")

  Lambda <- gi_renewal_force(E_star, g)   # the renewal force; this is the ONLY epi content
  if (seed_rate > 0) {                    # pre-day-0 seed: + seed_rate * shifted survival of g
    surv_g  <- rev(cumsum(rev(g))); shifted <- c(surv_g[-1L], 0)
    shifted <- c(shifted, rep(0, (max_time + 1L) - length(shifted)))[seq_len(max_time + 1L)]
    Lambda  <- Lambda + seed_rate * shifted
  }
  d_max  <- length(pi_EY)

  # --- Resolve first_rt_index (burn-in) --------------------------------------
  .date_to_idx <- function(dd) {
    if (inherits(dd, "Date")) {
      if (is.null(dates) || length(dates) != max_time + 1L)
        stop("first_rt_index is a Date but `dates` is missing / wrong length.")
      idx <- match(dd, as.Date(dates))
      if (is.na(idx)) stop(sprintf("first_rt_index = %s not in dates.", dd))
      return(idx - 1L)
    }
    as.integer(dd)
  }
  if (!is.null(first_rt_index)) {
    first_rt_index    <- .date_to_idx(first_rt_index)
    first_count_index <- first_rt_index + d_max
  } else {
    first_rt_index <- 0L; first_count_index <- 0L
  }
  if (first_count_index > max_time)
    stop(sprintf("first_count_index (%d) exceeds max_time (%d).", first_count_index, max_time))

  # --- Convolution matrix  X_{t,s} = Lambda_s * pi_EY[t-s] -------------------
  X_ts <- matrix(0.0, max_time + 1L, max_time + 1L)
  for (t in seq_len(max_time)) {
    s_lo <- max(0L, t - d_max)
    for (s in s_lo:(t - 1L)) {
      lag <- t - s
      X_ts[t + 1L, s + 1L] <- Lambda[s + 1L] * pi_EY[lag]
    }
  }

  # --- Offset (constant pre-day-0 exposures; usually zero) --------------------
  Y_all   <- as.numeric(inf_inc)
  off_all <- numeric(max_time + 1L)
  if (pre_hist_rate > 0) {
    surv_pi  <- rev(cumsum(rev(pi_EY)))
    surv_pad <- c(surv_pi, rep(0, max_time + 1L - length(surv_pi)))
    off_all  <- off_all + pre_hist_rate * c(0.0, surv_pad[seq_len(max_time)])
  }

  # --- Valid-row mask --------------------------------------------------------
  row_sums_X  <- rowSums(X_ts)
  burnin_mask <- seq_along(Y_all) > first_count_index
  valid_mask  <- (row_sums_X > 0) & (Y_all > y_min_count) & burnin_mask
  if (!any(valid_mask)) {
    valid_mask <- (row_sums_X > 0) & (Y_all > 0) & burnin_mask
    cat("  Warning: no points with Y > y_min_count; relaxed to Y > 0\n")
  }
  if (!any(valid_mask)) stop("No valid time points found.")

  Y_valid   <- Y_all[valid_mask]
  off_valid <- off_all[valid_mask]
  Xv        <- X_ts[valid_mask, , drop = FALSE]

  # --- Spline domain [fr[1], fr[2]] -----------------------------------------
  start_mask   <- (row_sums_X > 0) & (Y_all > y_min_start) & burnin_mask
  t_data_start <- if (any(start_mask)) as.integer(min(which(start_mask)) - 1L)
                  else as.integer(min(which(valid_mask)) - 1L)
  t_data_start <- max(t_data_start, first_rt_index)
  end_mask   <- (row_sums_X > 0) & (Y_all > y_min_end)
  t_data_end <- if (any(end_mask)) as.integer(max(which(end_mask)) - 1L)
                else as.integer(max(which(valid_mask)) - 1L)
  fr <- c(1L, t_data_end)
  if (fr[2] - fr[1] < 2L * knot_step)
    stop(sprintf("Spline domain [%d, %d] too narrow for knot_step=%d.", fr[1], fr[2], knot_step))

  spline_start <- fr[1]; spline_end <- fr[2]
  interior_knots <- seq(spline_start + knot_step, spline_end - knot_step, by = knot_step)
  t_eval <- 0:max_time
  ns_obj <- splines::ns(t_eval, knots = interior_knots,
                        Boundary.knots = c(spline_start, spline_end), intercept = TRUE)
  B <- as.matrix(ns_obj); n_basis <- ncol(B); np <- n_basis

  # --- Integrated 2nd-derivative^2 penalty (Omega) --------------------------
  n_fine <- max(1001L, as.integer(8L * (spline_end - spline_start) + 1L))
  t_fine <- seq(spline_start, spline_end, length.out = n_fine); dt <- t_fine[2] - t_fine[1]
  B_fine <- predict(ns_obj, t_fine)
  B_dd   <- (B_fine[1:(n_fine - 2L), , drop = FALSE] - 2 * B_fine[2:(n_fine - 1L), , drop = FALSE] +
             B_fine[3:n_fine, , drop = FALSE]) / dt^2
  Omega  <- dt * crossprod(B_dd)

  # --- optional value-jump (Heaviside) columns (unpenalized) ----------------
  jump_idx <- integer(0)
  if (!is.null(jump_times) && length(jump_times) > 0L) {
    jump_idx <- if (inherits(jump_times, "Date")) sort(unique(match(jump_times, as.Date(dates)) - 1L))
                else sort(unique(as.integer(jump_times)))
    jump_idx <- jump_idx[jump_idx > spline_start & jump_idx <= spline_end]
  }
  n_jump <- length(jump_idx)
  if (n_jump > 0L) {
    H <- vapply(jump_idx, function(xi) as.numeric(t_eval >= xi), numeric(length(t_eval)))
    B <- cbind(B, matrix(H, nrow = length(t_eval), ncol = n_jump)); np <- n_basis + n_jump
    Of <- matrix(0.0, np, np); Of[seq_len(n_basis), seq_len(n_basis)] <- Omega; Omega <- Of
  }

  Z_full       <- Xv %*% B
  P_pspline    <- Omega
  P_ridge      <- diag(c(rep(1.0, n_basis), rep(0.0, n_jump)), np, np)

  # --- Day-of-week (optional) ------------------------------------------------
  use_dow <- !is.null(dow_dates) && length(dow_dates) == (max_time + 1L)
  if (use_dow) {
    dow_all   <- as.integer(format(as.Date(dow_dates), "%u"))
    dow_valid <- dow_all[valid_mask]
  } else { dow_all <- NULL; dow_valid <- NULL }

  # --- Tail constraints: empty (natural spline => linear tails automatically) -
  empty_A <- matrix(0.0, nrow = 0L, ncol = np)
  tail_constraints <- list(none = empty_A, linear = empty_A, constant = empty_A)

  # --- Head-constant constraint (optional) -----------------------------------
  head_const_days <- as.integer(head_const_days)
  if (head_const_days > 0L) {
    anchor_row <- min(head_const_days + 1L, nrow(B))
    if (anchor_row >= 2L) {
      A_head <- sweep(B[1L:(anchor_row - 1L), , drop = FALSE], 2, as.numeric(B[anchor_row, ]), "-")
      head_const_extra_M <- 1e10 * crossprod(A_head)
    } else head_const_extra_M <- matrix(0.0, np, np)
  } else head_const_extra_M <- matrix(0.0, np, np)

  list(
    Z_full = Z_full, Y_valid = Y_valid, off_valid = off_valid,
    B = B, n_basis = n_basis, np = np, n_jump = n_jump, jump_times_idx = jump_idx,
    D2 = NULL, P_pspline = P_pspline, P_ridge = P_ridge,
    tail_constraints = tail_constraints,
    tail_const_extra_M = NULL, tail_const_days = NA_integer_,
    head_const_extra_M = head_const_extra_M, head_const_days = head_const_days,
    all_knots = NULL, int_knots = interior_knots, fr = fr, t_tail = NA_integer_,
    kink_vec = NULL, use_kink = FALSE, T_kink = NULL,
    valid_mask = valid_mask, t_data_start = t_data_start,
    first_rt_index = first_rt_index, first_count_index = first_count_index,
    max_time = max_time, spline_degree = 3L, knot_step = knot_step,
    use_dow = use_dow, dow_valid = dow_valid, dow_all = dow_all,
    basis = "natural", ns_obj = ns_obj,
    n_obs = sum(valid_mask), n_penalty = as.numeric(spline_end - spline_start),
    gi_force = Lambda
  )
}
