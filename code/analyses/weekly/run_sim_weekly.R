# code/analyses/weekly/run_sim_weekly.R
#
# Retrospective weekly-aggregated girt on the flu "wiggly" simulation (raw
# Poisson obs, known truth).  Also fits the canonical DAILY girt at lam_min for
# the daily-vs-weekly comparison.  Saves the fits the combined PDF consumes.
#
#   results/extra_methods/weekly/sim/weekly_fit.rds  (rt_at_lam_min, rt_at_lam_1se, ...)
#   results/extra_methods/weekly/sim/daily_fit.rds   (rt_daily with Rt_truth)

sd <- dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))
if (!length(sd) || is.na(sd)) sd <- "code/analyses/weekly"
source(file.path(sd, "_common.R"))
options(warn = 1)

out_dir <- file.path(wk_results_dir, "sim")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sim <- readRDS(wk_sim$sim_rds); sim$date <- as.Date(sim$date)
y     <- as.integer(round(pmax(0, sim[[wk_sim$obs_col]])))
truth <- as.numeric(sim[[wk_sim$truth_col]])
dates <- sim$date

g  <- sim_gi()
bd <- wk_build_daily_design(y, dates, g, wk_sim$beta_date, wk_sim$lik_date,
                            link = "identity")
design_d <- bd$design

attach_truth <- function(df) {
  df$Rt_truth <- truth[match(df$date, dates)]; df
}

# --- DAILY girt (lam_min) ----------------------------------------------------
cat("[sim] daily fit (cv min) ...\n")
sel_d <- gi_select_lambda_cv(design_d, wk_sim$lam_grid, tail = "linear",
                             nfold = 5L, cv_select_rule = "min",
                             error_measure = "deviance")
lam_d <- sel_d$best_lam_min
fit_d <- gi_solve(design_d, lam_d, tail = "linear")
rt_daily <- attach_truth(gi_rt_df_from_fit(fit_d, design_d, dates, lam_d,
                                           level = wk_cfg$level, overdispersion = TRUE))
saveRDS(list(rt_daily = rt_daily, lam = lam_d), file.path(out_dir, "daily_fit.rds"))
cat(sprintf("[sim] daily lam_min=%.3e\n", lam_d))

# --- WEEKLY girt (min + 1se) -------------------------------------------------
cat("[sim] weekly fit ...\n")
valid_dates <- dates[which(design_d$valid_mask)]
w <- fit_convrt_weekly_retrospective(
  design_d, lam_grid = wk_sim$lam_grid, dates = dates,
  dates_valid = valid_dates, cv_select_rule = "min",
  error_measure = "deviance", nfold = 5L, level = wk_cfg$level)
# Keep only the Rt curves the combined figure reads (drop heavy design objects).
saveRDS(list(rt_at_lam_min = attach_truth(w$rt_at_lam_min),
             rt_at_lam_1se = attach_truth(w$rt_at_lam_1se),
             lam_min = w$lam_min, lam_1se = w$lam_1se),
        file.path(out_dir, "weekly_fit.rds"))
cat(sprintf("[sim] weekly lam_min=%.3e lam_1se=%.3e  (%d weekly obs)\n",
            w$lam_min, w$lam_1se, nrow(w$design_weekly$Z_full)))
cat("[sim] done.\n")
