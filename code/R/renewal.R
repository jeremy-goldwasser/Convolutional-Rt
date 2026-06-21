# girt :: renewal.R
#
# The renewal-equation transmission force.
#
#   Lambda_s = sum_{k>=1} g_k * X_{s-k}
#
# the convolution of trailing exposure incidence X (= deconvolved infections)
# with the generation interval g.  Lambda_s is the expected new exposures at day
# s per unit R_s: the renewal equation is  E[X_s] = R_s * Lambda_s.  girt fits a
# smoothing spline for R_s through this relationship, then pushes the predicted
# exposures R_s * Lambda_s through the reporting model to the observation
# likelihood (see design.R).
#
# Truncation (k <= s-1, no same-day transmission) matches the convention used to
# build the design's reporting convolution.

gi_renewal_force <- function(X, g) {
  n <- length(X); L <- length(g); lam <- numeric(n)
  for (t in 2:n) { kr <- 1:min(L, t - 1L); lam[t] <- sum(g[kr] * X[t - kr]) }
  lam
}
