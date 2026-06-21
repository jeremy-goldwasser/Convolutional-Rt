# =============================================================================
# code/mechrt/realtime_conformal.R
#
# CANONICAL real-time conformal CI module for MechRt.
#
# Recipe (split-conformal, absolute score, strictly real-time):
#
#   At a target vintage T, given a cache of daily MechRt fits {R^{W'}}_{W' <= T}:
#     1. For each calibration vintage W' < T and horizon b in {0,...,MAX_LAG},
#        define the reference date t' = W' - b.  Keep scores where
#        t' <= T - BUFFER_DAYS so that R^T(t') is a stable pseudo-truth.
#     2. Score:  s(W', b) = | R^{W'}_{t'} - R^T(t') |.
#     3. Per-horizon conformal quantile:
#          q_b(T, alpha) = ceil((1 - alpha) * (n_b + 1))-th order statistic
#        where n_b = #{W' : eligible}.  q is NA when k > n_b.
#     4. CI at the target at date T-b:  R^T_{T-b} +/- q_b(T, alpha).
#     5. For dates beyond MAX_LAG behind the edge, leave the Wald CI alone
#        (conformal has no calibration there; far from the edge the Wald CI
#        is already well-calibrated).
#
# Public API:
#   gi_load_daily_vintage_cache(dir, variant = NA)
#       -> tibble(date, day, Rt_mean, Rt_lo, Rt_hi, vintage, variant)
#
#   gi_compute_conformal_q(daily_vintages, target_W,
#                       buffer_days, max_lag, alpha)
#       -> tibble(d_to_edge, n_cal, q)
#
#   gi_apply_conformal_to_rt_df(rt_df, daily_vintages, target_W,
#                            ci_level = 0.90,
#                            buffer_days = 14, max_lag = 13)
#       -> rt_df with Rt_lo/Rt_hi overwritten where conformal is available;
#          original Wald CI preserved as Rt_lo_wald/Rt_hi_wald;
#          new columns conformal_q, conformal_n_cal, ci_source ("conformal"
#          or "wald") for transparency.
#
# Caller notes:
#   - ci_level here is the NOMINAL coverage of the produced interval (e.g.
#     0.90).  Internally this becomes alpha = 1 - ci_level for the
#     conformal quantile.
#   - The function is idempotent: re-running on an already-conformalized
#     rt_df will re-apply (reads Rt_lo_wald if present, else treats current
#     Rt_lo/Rt_hi as Wald).
#   - For rows where no conformal q is available (insufficient cal history
#     at that horizon, or d_to_edge > max_lag), the Wald CI is kept and
#     ci_source = "wald".
#
# Dependencies: dplyr, tibble.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tibble)
})

# -----------------------------------------------------------------------------
gi_load_daily_vintage_cache <- function(dir, variant = NA_character_) {
  files <- sort(list.files(dir,
                           pattern = "^daily_\\d{4}-\\d{2}-\\d{2}\\.rds$",
                           full.names = TRUE))
  if (!length(files))
    stop(sprintf("No daily_<date>.rds files under %s", dir))
  out <- lapply(files, function(f) {
    r  <- readRDS(f)
    df <- r$rt_df
    vd <- r$meta$vintage
    if (is.null(vd)) vd <- r$meta$vintage_date
    if (is.null(vd))
      stop(sprintf("Daily fit %s missing meta$vintage(_date)", f))
    df$vintage <- as.Date(vd)
    df$variant <- variant
    df$date    <- as.Date(df$date)
    df
  })
  bind_rows(out)
}

# -----------------------------------------------------------------------------
.gi_conformal_q_one <- function(s, alpha) {
  s <- s[is.finite(s)]
  n <- length(s)
  if (n == 0L) return(c(NA_real_, 0L))
  k <- ceiling((1 - alpha) * (n + 1))
  q <- if (k > n) NA_real_ else sort(s)[k]
  c(q, n)
}

gi_compute_conformal_q <- function(daily_vintages, target_W,
                                buffer_days, max_lag, alpha) {
  target_W <- as.Date(target_W)
  cal <- daily_vintages |>
    filter(vintage < target_W) |>
    mutate(d_to_edge = as.integer(vintage - date)) |>
    filter(d_to_edge >= 0L, d_to_edge <= as.integer(max_lag),
           date <= target_W - as.integer(buffer_days))

  # Anchor = R^T's daily fit, on dates <= target_W - buffer_days
  rtT <- daily_vintages |>
    filter(vintage == target_W) |>
    transmute(date, Rt_anchor = Rt_mean)
  if (!nrow(rtT) || !nrow(cal))
    return(tibble(d_to_edge = integer(0), n_cal = integer(0), q = numeric(0)))

  cal <- cal |>
    inner_join(rtT, by = "date") |>
    transmute(d_to_edge, score = abs(Rt_mean - Rt_anchor))

  cal |>
    group_by(d_to_edge) |>
    summarise(qn = list(.gi_conformal_q_one(score, alpha)),
              .groups = "drop") |>
    mutate(q     = vapply(qn, `[`, numeric(1), 1L),
           n_cal = vapply(qn, function(x) as.integer(x[2L]), integer(1))) |>
    select(d_to_edge, n_cal, q)
}

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Multi-level conformal half-widths.  Returns a long tibble
#   (d_to_edge, ci_level, n_cal, q)
# with q = conformal half-width at level (1 - alpha) for each horizon b.
# Useful for evaluation where you need per-alpha bands (e.g. CE).
gi_conformal_q_multilevel <- function(daily_vintages, target_W,
                                   buffer_days, max_lag, ci_levels) {
  target_W <- as.Date(target_W)
  rtT <- daily_vintages |>
    filter(vintage == target_W) |>
    transmute(date, Rt_anchor = Rt_mean)
  cal <- daily_vintages |>
    filter(vintage < target_W) |>
    mutate(d_to_edge = as.integer(vintage - date)) |>
    filter(d_to_edge >= 0L, d_to_edge <= as.integer(max_lag),
           date <= target_W - as.integer(buffer_days)) |>
    inner_join(rtT, by = "date") |>
    transmute(d_to_edge, score = abs(Rt_mean - Rt_anchor))
  if (!nrow(cal))
    return(tibble(d_to_edge = integer(0), ci_level = numeric(0),
                  n_cal = integer(0), q = numeric(0)))
  cal |>
    group_by(d_to_edge) |>
    summarise(scores = list(score), n_cal = dplyr::n(), .groups = "drop") |>
    tidyr::crossing(ci_level = ci_levels) |>
    rowwise() |>
    mutate(q = .gi_conformal_q_one(scores, 1 - ci_level)[1L]) |>
    ungroup() |>
    select(d_to_edge, ci_level, n_cal, q)
}

gi_apply_conformal_to_rt_df <- function(rt_df, daily_vintages, target_W,
                                     ci_level = 0.90,
                                     buffer_days = 14L,
                                     max_lag = 13L) {
  target_W <- as.Date(target_W)
  alpha    <- 1 - ci_level

  rt_df$date <- as.Date(rt_df$date)
  # Snapshot Wald CI on the first call (idempotent on later calls).
  if (!"Rt_lo_wald" %in% names(rt_df)) rt_df$Rt_lo_wald <- rt_df$Rt_lo
  if (!"Rt_hi_wald" %in% names(rt_df)) rt_df$Rt_hi_wald <- rt_df$Rt_hi

  q_df <- gi_compute_conformal_q(daily_vintages, target_W,
                              buffer_days = buffer_days,
                              max_lag     = max_lag,
                              alpha       = alpha)

  out <- rt_df |>
    mutate(d_to_edge = as.integer(target_W - date)) |>
    left_join(q_df, by = "d_to_edge") |>
    mutate(have_conf = is.finite(q) &
                       d_to_edge >= 0L & d_to_edge <= as.integer(max_lag),
           Rt_lo = ifelse(have_conf, Rt_mean - q, Rt_lo_wald),
           Rt_hi = ifelse(have_conf, Rt_mean + q, Rt_hi_wald),
           conformal_q     = q,
           conformal_n_cal = n_cal,
           ci_source = ifelse(have_conf, "conformal", "wald")) |>
    select(-q, -n_cal, -have_conf)
  out
}
