# code/analyses/weekly/run_s2_weekly_curves.R
#
# Real-time weekly-aggregated girt on US flu hosps, season 2.  At each cutoff:
#   1. build the daily renewal design truncated through the cutoff (no DoW, LOG
#      deconv, GI = discrete Gamma(3.2,1.6));
#   2. fit_convrt_weekly_realtime -> FV-tuned gamma;
#   3. refit at EVERY gamma in gamma_grid, extracting the Rt curve at each.
# Also a retrospective weekly fit over the full season for the overlay.
#
#   results/extra_methods/weekly/s2/weekly_combined_rt_curves_data.rds
#   results/extra_methods/weekly/s2/weekly_fit.rds   (retro weekly rt_df)

sd <- dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))
if (!length(sd) || is.na(sd)) sd <- "code/analyses/weekly"
source(file.path(sd, "_common.R"))
options(warn = 1)

cfg     <- wk_real_s2
out_dir <- file.path(wk_results_dir, "s2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

g       <- real_gi()
pi_EH   <- pi_EH_obj()$pmf
hosps   <- load_flu_hosps()

win_for <- function(end_date) {
  w <- hosps[hosps$date >= cfg$beta_date - 30L & hosps$date <= end_date, ]
  w$day <- seq_len(nrow(w)) - 1L; w
}

# --- Retrospective weekly fit over the full season (overlay) -----------------
cat("[s2] retro weekly fit (full season) ...\n")
wr  <- win_for(cfg$eos)
bdr <- wk_build_daily_design(wr$hosps_1d, wr$date, g, cfg$beta_date,
                             lik_date = cfg$beta_date, link = "log")
retro <- fit_convrt_weekly_retrospective(
  bdr$design, lam_grid = cfg$lam_grid, dates = wr$date,
  dates_valid = wr$date[which(bdr$design$valid_mask)],
  cv_select_rule = "1se", error_measure = "deviance",
  nfold = cfg$cv_folds, level = wk_cfg$level)
# Keep only the retro Rt curve the combined figure overlays.
saveRDS(list(rt_df = retro$rt_df, lam = retro$lam),
        file.path(out_dir, "weekly_fit.rds"))
cat(sprintf("[s2] retro weekly lam_1se=%.3e (%d weeks)\n",
            retro$lam, nrow(retro$design_weekly$Z_full)))

# --- Per-cutoff: FV gamma + full gamma sweep ---------------------------------
build_one_cutoff <- function(week_end) {
  w  <- win_for(week_end)
  bd <- wk_build_daily_design(w$hosps_1d, w$date, g, cfg$beta_date,
                              lik_date = cfg$beta_date, link = "log")
  d_d <- bd$design
  valid_dates <- w$date[which(d_d$valid_mask)]

  res <- fit_convrt_weekly_realtime(
    d_d, dates = w$date, lam_grid = cfg$lam_grid, use_taper = TRUE, pi_E = pi_EH,
    gamma = NULL, gamma_grid = cfg$gamma_grid, dates_valid = valid_dates,
    n_fv = cfg$n_fv, cv_select_rule = "1se", error_measure = "deviance",
    nfold = cfg$cv_folds, tail = "linear", fv_error = "deviance",
    level = wk_cfg$level, verbose = FALSE)

  lam_chosen <- res$lam; gamma_min <- res$gamma
  d_w <- res$design_weekly; P_taper <- res$P_taper

  rt_per_gamma <- vector("list", length(cfg$gamma_grid))
  warm <- res$fit_main$theta
  for (i in seq_along(cfg$gamma_grid)) {
    gm <- cfg$gamma_grid[i]
    fit_g <- tryCatch(
      if (gm > 0) ConvRt:::.gi_solve_taper(d_w, lam_chosen, lam_taper = gm, P_taper = P_taper,
                                  tail = "linear", theta_init = warm)
      else gi_solve(d_w, lam_chosen, tail = "linear", theta_init = warm),
      error = function(e) NULL)
    if (is.null(fit_g) || !all(is.finite(fit_g$theta))) next
    warm <- fit_g$theta
    rt_g <- gi_extract_rt(fit_g, d_w, lam = lam_chosen, lam_taper = gm,
                          P_taper = P_taper, level = wk_cfg$level, overdispersion = TRUE)
    if (!nrow(rt_g)) next
    # only the columns the figure uses (drop pointwise CI from the 18-curve fan)
    rt_per_gamma[[i]] <- tibble(date = w$date[rt_g$day + 1L], Rt_mean = rt_g$Rt_mean,
                                gamma = gm, gamma_index = i)
  }
  rt_per_gamma_df <- bind_rows(rt_per_gamma)

  cat(sprintf("[s2 | %s] lam=%.2e gamma_min=%.2e n_curves=%d\n",
              format(week_end), lam_chosen, gamma_min,
              length(unique(rt_per_gamma_df$gamma))))
  list(week_end_date = week_end, lam_chosen = lam_chosen, gamma_min = gamma_min,
       rt_per_gamma = rt_per_gamma_df, gamma_grid = cfg$gamma_grid)
}

cutoff_data <- lapply(cfg$cutoffs, build_one_cutoff)
names(cutoff_data) <- format(cfg$cutoffs)
saveRDS(cutoff_data, file.path(out_dir, "weekly_combined_rt_curves_data.rds"))
cat("[s2] done.\n")
