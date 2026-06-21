# girt :: delays.R
#
# Delay-distribution utilities for the generation-interval Rt method.
#
# girt needs exactly two epidemiological inputs (besides the data):
#   * a generation interval g  (the renewal-equation transmission kernel), and
#   * a reporting delay pi_EY  (exposure -> observed outcome, e.g. E -> hosp).
# Neither is "mechanistic": g is a single PMF you can take from the literature
# or estimate externally, and pi_EY is the observation/reporting model that ANY
# observation-level Rt method needs.
#
# `gi_from_compartmental()` is an OPTIONAL convenience for users who want the
# generation interval IMPLIED by an SEIR-style compartmental model (latent E->I
# delay + infectious period).  It is the only place anything compartmental
# appears, and it produces a plain GI PMF -- the rest of girt never sees the
# compartments.

# Discrete Gamma delay PMF on {1, 2, ...} (midpoint-binned, renormalized).
.gi_gamma_shape_rate <- function(m, s) list(shape = (m / s)^2, rate = m / s^2)

gi_discrete_gamma_delay <- function(mean, sd) {
  sr  <- .gi_gamma_shape_rate(mean, sd)
  d   <- ceiling(qgamma(0.9995, shape = sr$shape, rate = sr$rate))
  up  <- seq_len(d) + 0.5
  lo  <- c(0, seq_len(d - 1L) + 0.5)
  pmf <- pgamma(up, shape = sr$shape, rate = sr$rate) -
         pgamma(lo, shape = sr$shape, rate = sr$rate)
  pmf <- pmf / sum(pmf)
  list(pmf = pmf, shape = sr$shape, rate = sr$rate)
}

# Infectious kernel zeta^EI_k = P(infectious k days after exposure)
#   = sum_j pi_lat[j] * P(W^IR > k - j),  P(W^IR > 0) = 1.
# surv_IR[a] = P(W^IR >= a) = P(W^IR > a-1); pass surv_IR = c(1, 1-cumsum(pi_IR)).
.gi_infectious_kernel <- function(pi_lat, surv_IR) {
  K <- numeric(length(pi_lat) + length(surv_IR) - 1L)
  for (j in seq_along(pi_lat))
    for (a in seq_along(surv_IR))
      K[j + a - 1L] <- K[j + a - 1L] + pi_lat[j] * surv_IR[a]
  K
}

# OPTIONAL: the generation interval implied by a two-stage (SEIR-style)
# compartmental model.  zeta^EI_k = P(infectious k days after exposure) is the
# infectious kernel; the GI is zeta^EI rescaled (and possibly lag-shifted) by
# the convention of whichever MechRt variant we're trying to reproduce.
#
#   pi_lat                  : E->I (latent) delay pmf
#   pi_IR                   : infectious-period pmf
#   mean_infectious_nominal : if NULL (default), divide by sum(zeta) = mu^IR so
#                             g sums to 1 -- the canonical MechRt-equivalent
#                             GI (g_k = P_still_I / mu^IR).  If a number, divide
#                             by that instead (matches the flu / rtestim / covidestim
#                             MechRt convention of dividing by the nominal 1/gamma;
#                             resulting g does NOT sum to 1 in general).
#   lag_one                 : if TRUE, prepend a leading zero so g starts at
#                             lag-1 (matches the covidestim MechRt convention).
#
# Returns: list(g, mu_IR, kernel).
gi_from_compartmental <- function(pi_lat, pi_IR,
                                  mean_infectious_nominal = NULL,
                                  lag_one = FALSE) {
  surv_IR <- c(1, 1 - cumsum(pi_IR))
  zeta    <- .gi_infectious_kernel(pi_lat, surv_IR)
  mu_IR   <- sum(zeta)
  scale   <- if (is.null(mean_infectious_nominal)) mu_IR else as.numeric(mean_infectious_nominal)
  g_core  <- zeta / scale
  g       <- if (isTRUE(lag_one)) c(0, g_core) else g_core
  list(g = g, mu_IR = mu_IR, kernel = zeta)
}
