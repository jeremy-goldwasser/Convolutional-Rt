# girt :: realtime.R
#
# Design post-processing + Rt-curve extraction helpers shared by the
# retrospective and real-time wrappers.  These operate purely on the design
# object and are method-agnostic (copied from the MechRt real-time helpers,
# renamed for the girt namespace).

suppressPackageStartupMessages({ library(tibble); library(dplyr) })

# Resolve a scalar-or-vector severity into a length-n vector.
gi_resolve_severity <- function(n, severity = 1) {
  if (length(severity) == 1L) return(rep(as.numeric(severity), n))
  p <- as.numeric(severity)
  if (length(p) != n) stop(sprintf("severity must have length 1 or %d", n))
  p
}

# Multiply the design's valid rows by the per-day severity rho (the observation
# scale): mu_t = rho_t * sum_s R_s X_{t,s}.
gi_apply_severity_to_design <- function(design, severity_by_day) {
  valid_rows <- which(design$valid_mask)
  p_valid    <- severity_by_day[valid_rows]
  design$Z_full    <- design$Z_full    * p_valid
  design$off_valid <- design$off_valid * p_valid
  design
}

# Drop design rows whose day < min_day (likelihood-window start filter).
gi_enforce_likelihood_start <- function(design, valid_days, min_day) {
  valid_rows <- which(design$valid_mask)
  keep <- valid_days >= min_day
  if (!any(keep)) stop("No valid design rows after likelihood start filter")
  keep_rows <- valid_rows[keep]
  design$Y_valid   <- design$Y_valid[keep]
  design$off_valid <- design$off_valid[keep]
  design$Z_full    <- design$Z_full[keep, , drop = FALSE]
  new_mask <- rep(FALSE, length(design$valid_mask)); new_mask[keep_rows] <- TRUE
  design$valid_mask    <- new_mask
  design$valid_row_idx <- keep_rows
  design$valid_days    <- valid_days[keep]
  if (isTRUE(design$use_dow) && !is.null(design$dow_valid))
    design$dow_valid <- design$dow_valid[keep]
  design
}

# Extract the R_t curve (+ pointwise CI) from a fit, attaching calendar dates.
gi_rt_df_from_fit <- function(fit, design, dates, lam,
                              lam_taper = 0, P_taper = NULL, level = 0.95,
                              overdispersion = TRUE) {
  rt <- gi_extract_rt(fit, design, lam = lam,
                      lam_taper = lam_taper, P_taper = P_taper,
                      level = level, overdispersion = overdispersion)
  if (nrow(rt) == 0) return(rt)
  rt$date <- as.Date(dates)[rt$day + 1L]
  rt[, c("date", "day", "Rt_mean", "Rt_lo", "Rt_hi")]
}
