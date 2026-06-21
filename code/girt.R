# =============================================================================
# girt :: top-level loader
#
# girt is a generation-interval Rt estimator: a regularized identity-link
# Poisson regression that fits a smoothing-spline R_t through the RENEWAL
# EQUATION and predicts the observed counts.  It is the non-mechanistic sibling
# of mechrt -- same observation model and same penalized-spline engine, but the
# per-day transmission force comes from a single generation interval g rather
# than from compartmental (E->I, I->R) distributions.
#
# Usage:
#   source("path/to/code/girt.R")
#
# -----------------------------------------------------------------------------
# PUBLIC API
#   gi_discrete_gamma_delay(mean, sd)      discrete Gamma delay pmf
#   gi_from_compartmental(pi_lat, pi_IR,   OPTIONAL: GI implied by an SEIR model.
#       mean_infectious_nominal = NULL,    Default sum-normalizes zeta^EI (g sums
#       lag_one = FALSE)                   to 1, the canonical MechRt-equivalent
#                                          form).  Pass mean_infectious_nominal to
#                                          divide by the nominal 1/gamma instead
#                                          (flu/rtestim MechRt convention);
#                                          add lag_one = TRUE for the covidestim
#                                          MechRt convention (prepend zero).
#   gi_deconvolve_exposures(...)           obs -> exposure incidence X (reporting
#                                          delay only; no SEIR kernel)
#   fit_girt_retrospective(...)            end-of-season Rt + pointwise CI
#   fit_girt_realtime(...)                 right-censored Rt + tapered tail + CI
#   gi_aggregate_design_weekly(...)        collapse daily design -> weekly Poisson
#   fit_girt_weekly_retrospective(...)     weekly-aggregated retro Rt (min + 1se)
#   fit_girt_weekly_realtime(...)          weekly-aggregated real-time Rt + taper
#
# Lower-level building blocks (used by the wrappers; call directly for custom
# pipelines): gi_renewal_force, build_gi_design, gi_solve, gi_select_lambda_cv,
#   gi_extract_rt, gi_build_tapered_penalty, gi_tune_gamma_fv, .gi_solve_taper,
#   gi_apply_severity_to_design, gi_enforce_likelihood_start, gi_rt_df_from_fit,
#   and the split-conformal helpers in R/conformal.R
#   (gi_conformal_q_multilevel, gi_apply_conformal_to_rt_df, ...).
#
# -----------------------------------------------------------------------------
# Module layout (sourced in order):
#   R/delays.R         delay pmfs + gi_from_compartmental
#   R/deconvolution.R  gi_deconvolve_exposures (exposures-only)
#   R/renewal.R        gi_renewal_force
#   R/design.R         build_gi_design (natural-spline renewal design)
#   R/solve.R          gi_solve (KKT IRLS)
#   R/lambda_select.R  gi_select_lambda_cv (+ _dow)
#   R/extract.R        gi_extract_rt, gi_estimate_dispersion, simband
#   R/taper.R          gi_build_tapered_penalty, .gi_solve_taper, gi_tune_gamma_fv
#   R/realtime.R       severity / likelihood-start / rt_df helpers
#   R/conformal.R      split-conformal real-time CIs
#   R/wrappers.R       fit_girt_retrospective, fit_girt_realtime
#   R/weekly.R         weekly aggregation + weekly retro/real-time wrappers
#
# The files under R/ named solve/lambda_select/extract/taper/conformal are the
# shared penalized-spline-Poisson-GLM engine, copied from mechrt and renamed to
# the gi_* namespace so girt is fully self-contained (sources nothing from
# mechrt).  The method-defining files (delays, deconvolution, renewal, design,
# wrappers) are girt-native.
# =============================================================================

.girt_dir <- (function() {
  this <- character(0)
  if (exists("sys.frames", mode = "function")) {
    frms <- sys.frames()
    for (i in rev(seq_along(frms))) {
      sf <- try(get0("ofile", envir = frms[[i]], inherits = FALSE), silent = TRUE)
      if (!inherits(sf, "try-error") && is.character(sf) && length(sf) == 1L &&
          grepl("girt\\.R$", sf)) { this <- sf; break }
    }
  }
  if (length(this) == 0L) {
    args <- commandArgs(trailingOnly = FALSE)
    this <- sub("^--file=", "", args[grep("^--file=", args)])
  }
  if (length(this) == 0L || identical(this, "")) {
    env <- Sys.getenv("GIRT_DIR", unset = "")
    if (nzchar(env)) return(normalizePath(env, mustWork = FALSE))
    return(normalizePath(file.path(getwd(), "girt"), mustWork = FALSE))
  }
  normalizePath(dirname(this), mustWork = FALSE)
})()

for (.mod in c("delays.R", "deconvolution.R", "renewal.R", "design.R",
               "solve.R", "lambda_select.R", "extract.R", "taper.R",
               "realtime.R", "conformal.R", "tf.R", "wrappers.R",
               "weekly.R")) {
  source(file.path(.girt_dir, "R", .mod), local = FALSE)
}
rm(.mod)
