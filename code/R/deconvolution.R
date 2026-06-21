# girt :: deconvolution.R
#
# Recover daily exposure incidence X_t (a.k.a. infections) from observed counts
# by deconvolving ONLY the reporting delay:
#
#   Y_t ~ Poisson( rho * sum_k pi_EY[k] X_{t-k} )
#
# Two links (both mechanism-free: X_hat depends only on the reporting delay pi_EY
# and the severity rho -- NOT on any latent/infectious-period distribution):
#
#   link = "identity" (DEFAULT): X = B theta on a natural spline with a
#     magnitude-weighted integrated 2nd-derivative penalty.  The correctly
#     specified DGP link; the girt default.
#   link = "log": X = exp(B theta) on a cr regression spline with a uniform
#     2nd-derivative penalty (mgcv).  The legacy front-end; provided so girt can
#     reproduce log-deconvolution analyses (e.g. the canonical real-flu fits).
#
# Both reproduce MechRt's deconvolve_prevalence X_hat for the matching link
# (the identity path = .deconv_identity_weighted; the log path = .deconv_legacy,
# log), with the prevalence/SEIR-kernel convolution removed.  lambda by GCV.
# Returns list(X_hat, lambda_gcv).

library(splines)

gi_deconvolve_exposures <- function(obs_inc, pi_EY, severity_rate,
                                    link = c("identity", "log"),
                                    burn_in = 30L, knot_step = 5L, floor_frac = 0.03,
                                    K_basis = NULL, lambda_grid = NULL, verbose = FALSE) {
  link <- match.arg(link)
  if (link == "log")
    return(.gi_deconv_log(obs_inc, pi_EY, severity_rate, burn_in = burn_in,
                          K_basis = K_basis,
                          lambda_grid = if (is.null(lambda_grid)) 10^seq(-4, 6, length.out = 40) else lambda_grid,
                          verbose = verbose))
  .gi_deconv_identity(obs_inc, pi_EY, severity_rate, burn_in = burn_in,
                      knot_step = knot_step, floor_frac = floor_frac,
                      lambda_grid = if (is.null(lambda_grid)) 10^seq(-10, 4, length.out = 60) else lambda_grid,
                      verbose = verbose)
}

# --- identity-link, magnitude-weighted smoothing spline (DEFAULT) -------------
.gi_deconv_identity <- function(obs_inc, pi_EY, severity_rate, burn_in = 30L,
                                knot_step = 5L, floor_frac = 0.03,
                                lambda_grid = 10^seq(-10, 4, length.out = 60), verbose = FALSE) {
  y <- as.numeric(obs_inc); n <- length(y); rho <- severity_rate; tg <- seq_len(n)
  iknots <- seq(knot_step + 1L, n - knot_step, by = knot_step)
  ns_obj <- splines::ns(tg, knots = iknots, Boundary.knots = c(1, n), intercept = TRUE)
  B <- as.matrix(ns_obj); p <- ncol(B)
  n_fine <- max(1001L, as.integer(8L * (n - 1L) + 1L))
  t_fine <- seq(1, n, length.out = n_fine); dt <- t_fine[2] - t_fine[1]
  Bf <- predict(ns_obj, t_fine)
  B_dd <- (Bf[1:(n_fine - 2L), , drop = FALSE] - 2 * Bf[2:(n_fine - 1L), , drop = FALSE] +
           Bf[3:n_fine, , drop = FALSE]) / dt^2
  mid <- t_fine[2:(n_fine - 1L)]
  L_EY <- length(pi_EY); A <- matrix(0, n, n)
  for (t in 2:n) { kr <- 1:min(L_EY, t - 1); A[t, t - kr] <- pi_EY[kr] }
  Lbar <- round(sum(seq_along(pi_EY) * pi_EY))
  ma <- vapply(seq_len(n), function(t) mean(y[max(1L, t - 3L):min(n, t + 3L)]), numeric(1))
  xtil <- ma[pmin(seq_len(n) + Lbar, n)]
  pk <- max(xtil); v_daily <- (pk / pmax(xtil, floor_frac * pk))^2
  wf <- approx(seq_along(v_daily), v_daily, xout = mid, rule = 2)$y
  S_w <- dt * crossprod(sqrt(wf) * B_dd)
  v_idx <- (burn_in + 1):n; n_eff <- length(v_idx); yv <- y[v_idx]; AB <- A %*% B
  nll <- function(th, lam) { mu <- rho * as.numeric(AB %*% th); muv <- pmax(mu[v_idx], 1e-12)
    -sum(yv * log(muv) - muv) + lam * drop(crossprod(th, S_w %*% th)) }
  gr <- function(th, lam) { mu <- rho * as.numeric(AB %*% th); muv <- pmax(mu[v_idx], 1e-12)
    -rho * as.numeric(crossprod(AB[v_idx, , drop = FALSE], yv / muv - 1)) + 2 * lam * as.numeric(S_w %*% th) }
  th <- as.numeric(solve(crossprod(B) + 1e-3 * diag(p), crossprod(B, rep(mean(y[y > 0]) / rho, n))))
  gcvs <- numeric(length(lambda_grid)); thetas <- vector("list", length(lambda_grid))
  for (i in seq_along(lambda_grid)) {
    o <- optim(th, nll, gr, lam = lambda_grid[i], method = "BFGS",
               control = list(maxit = 800, reltol = 1e-10)); th <- o$par
    mu <- rho * as.numeric(AB %*% th); muv <- pmax(mu[v_idx], 1e-12)
    Jv <- rho * AB[v_idx, , drop = FALSE]; JtWJ <- crossprod(Jv / sqrt(muv)); Mpen <- JtWJ + lambda_grid[i] * S_w
    eps_r <- max(1e-8 * max(diag(Mpen), 1), .Machine$double.eps)
    edf <- tryCatch(sum(diag(solve(Mpen + eps_r * diag(p), JtWJ))), error = function(e) NA_real_)
    dev <- 2 * sum(ifelse(yv > 0, yv * log(yv / muv), 0) - (yv - muv))
    gcvs[i] <- dev * n_eff / (n_eff - edf)^2; thetas[[i]] <- th
  }
  i_best <- which.min(gcvs)
  X_hat <- pmax(as.numeric(B %*% thetas[[i_best]]), 0)
  if (verbose) cat(sprintf("gi_deconvolve_exposures[identity]: n=%d GCV lam=%.3g peak X=%.0f\n",
                           n, lambda_grid[i_best], max(X_hat)))
  list(X_hat = X_hat, lambda_gcv = lambda_grid[i_best])
}

# --- log-link, cr regression spline (LEGACY option; needs mgcv) ---------------
# Faithful copy of MechRt's .deconv_legacy(link = "log") X_hat path.
.gi_deconv_log <- function(obs_inc, pi_EY, severity_rate, burn_in = 30L,
                           K_basis = NULL, lambda_grid = 10^seq(-4, 6, length.out = 40),
                           verbose = FALSE) {
  if (!requireNamespace("mgcv", quietly = TRUE))
    stop("gi_deconvolve_exposures(link='log') requires the mgcv package")
  y <- as.numeric(obs_inc); n <- length(y); tg <- seq_len(n)
  if (is.null(K_basis)) K_basis <- max(10L, min(n, as.integer(ceiling(n / 5)) + 4L))
  sm <- mgcv::smoothCon(mgcv::s(tg, bs = "cr", k = K_basis),
                        data = data.frame(tg = tg), absorb.cons = FALSE, scale.penalty = FALSE)[[1]]
  B_sm <- sm$X; S_pen <- sm$S[[1]]; p <- ncol(B_sm)
  L_EY <- length(pi_EY); A <- matrix(0, n, n)
  for (t in 2:n) { kr <- 1:min(L_EY, t - 1); A[t, t - kr] <- pi_EY[kr] }
  v_idx <- (burn_in + 1):n; n_eff <- length(v_idx)
  X_from_theta <- function(theta) exp(as.numeric(B_sm %*% theta))
  nll <- function(theta, lambda) { X <- X_from_theta(theta); mu <- severity_rate * as.numeric(A %*% X)
    muv <- pmax(mu[v_idx], 1e-12); -sum(y[v_idx] * log(muv) - muv) + lambda * drop(t(theta) %*% S_pen %*% theta) }
  gr <- function(theta, lambda) { X <- X_from_theta(theta); mu <- severity_rate * as.numeric(A %*% X)
    r <- numeric(n); r[v_idx] <- severity_rate * (y[v_idx] / pmax(mu[v_idx], 1e-12) - 1)
    -as.numeric(crossprod(B_sm, X * as.numeric(crossprod(A, r)))) + 2 * lambda * as.numeric(S_pen %*% theta) }
  fit_at <- function(lambda, theta0) {
    fit <- optim(theta0, nll, gr = gr, lambda = lambda, method = "BFGS",
                 control = list(maxit = 800, reltol = 1e-10))
    th <- fit$par; X <- X_from_theta(th); mu <- severity_rate * as.numeric(A %*% X)
    G <- severity_rate * (A %*% (X * B_sm)); Gv <- G[v_idx, , drop = FALSE]; muv <- pmax(mu[v_idx], 1e-12)
    H <- crossprod(Gv / sqrt(muv)); eps_r <- max(1e-10 * max(diag(H), 1.0), .Machine$double.eps)
    edf <- sum(diag(solve(H + lambda * S_pen + eps_r * diag(nrow(H)), H)))
    yv <- y[v_idx]; dev <- 2 * sum(ifelse(yv > 0, yv * log(yv / muv), 0) - (yv - muv))
    list(theta = th, X = X, gcv = dev * n_eff / (n_eff - edf)^2)
  }
  theta_warm <- rep(log(mean(y[y > 0]) / severity_rate), p)
  gcvs <- numeric(length(lambda_grid)); thetas <- vector("list", length(lambda_grid))
  for (i in seq_along(lambda_grid)) {
    r <- fit_at(lambda_grid[i], theta_warm); theta_warm <- r$theta
    gcvs[i] <- r$gcv; thetas[[i]] <- r$theta
  }
  i_best <- which.min(gcvs)
  X_hat <- X_from_theta(thetas[[i_best]])
  if (verbose) cat(sprintf("gi_deconvolve_exposures[log]: n=%d GCV lam=%.3g peak X=%.0f\n",
                           n, lambda_grid[i_best], max(X_hat)))
  list(X_hat = X_hat, lambda_gcv = lambda_grid[i_best])
}
