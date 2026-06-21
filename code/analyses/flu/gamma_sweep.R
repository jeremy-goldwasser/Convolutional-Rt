# code/analyses/flu/gamma_sweep.R  --  girt gamma-sweep data for the gamma figures.
#   Rscript gamma_sweep.R <season> <date> <mode>
#     mode: realtime_cv1se (->_1se) | retro_self_notaper (->_retrolam) | s1_retro (->_s1retrolam)
# Writes results/real/flu/girt/<season>/real_time/gamma/<date>/data<suffix>.rds
# matching the contract of realtime_gamma_sweep.R (girt-log).

sd <- dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))
if (!length(sd) || is.na(sd)) sd <- "code/analyses/flu"
source(file.path(sd, "_common.R")); source(file.path(sd, "_realtime.R"))
suppressPackageStartupMessages({ library(dplyr); library(tibble) })
args <- commandArgs(trailingOnly = TRUE)
season <- args[1]; WEEK_END <- as.Date(args[2]); MODE <- args[3]
suffix <- switch(MODE, realtime_cv1se = "_1se", retro_self_notaper = "_retrolam", s1_retro = "_s1retrolam")

pEH    <- gi_discrete_gamma_delay(flu_cfg$mean_hosp, flu_cfg$sd_hosp); pi_EH <- pEH$pmf
g      <- flu_gi()
hosps  <- load_flu_hosps(); beta <- flu_seasons[[season]]$beta
retro_lam <- function(s) readRDS(list.files(file.path(flu_results_dir, "girt", s, "retrospective", "results"), full.names = TRUE)[1])$lam_1se
retro_rt  <- function(s) { r <- readRDS(list.files(file.path(flu_results_dir, "girt", s, "retrospective", "results"), full.names = TRUE)[1])$rt_df
                           data.frame(date = as.Date(r$date), Rt_mean = r$Rt_mean) }

win <- hosps[hosps$date >= beta - 30L & hosps$date <= WEEK_END, ]; win$day <- seq_len(nrow(win)) - 1L; end_t <- nrow(win) - 1L
X_hat <- gi_deconvolve_exposures(win$hosps_1d, pi_EY = pi_EH, severity_rate = flu_cfg$severity_rate, link = "log", burn_in = 30L)$X_hat
d <- suppressWarnings(build_gi_design(E_star = X_hat, g = g, pi_EY = pi_EH, max_time = end_t, inf_inc = win$hosps_1d,
      y_min_count = 0, y_min_start = 0, y_min_end = 0, knot_step = flu_cfg$knot_step,
      first_rt_index = beta, dates = win$date, dow_dates = win$date))
d <- gi_enforce_likelihood_start(d, win$day[which(d$valid_mask)], win$day[win$date == beta])
d <- gi_apply_severity_to_design(d, rep(flu_cfg$severity_rate, nrow(win)))
valid_dates <- win$date[which(d$valid_mask)]

LAM_USED <- switch(MODE,
  realtime_cv1se = { s <- gi_select_lambda_cv_dow(d, valid_dates, flu_cfg$ns_lambda_grid, tail = "linear",
                          nfold = flu_cfg$cv_folds, cv_select_rule = "1se", error_measure = "deviance")
                     if (!is.null(s$best_lam_1se) && is.finite(s$best_lam_1se)) s$best_lam_1se else s$best_lam },
  retro_self_notaper = retro_lam(season),
  s1_retro = retro_lam("s1"))

fit_main <- gi_solve(d, LAM_USED, tail = "linear"); dow <- fit_main$dow_effects
d_sc <- d
if (!is.null(dow) && isTRUE(d_sc$use_dow)) { dt <- dow[d_sc$dow_valid]
  d_sc$Z_full <- d_sc$Z_full * dt; d_sc$off_valid <- d_sc$off_valid * dt; d_sc$use_dow <- FALSE; d_sc$dow_valid <- NULL }
P_taper <- gi_build_tapered_penalty(d_sc, pi_EH, quiet = TRUE)$P_taper
fv <- .flu_fv_tune_gamma(win, X_hat, g, pi_EH, pEH$shape, pEH$rate, LAM_USED, dow, beta, beta,
                         flu_cfg$severity_rate, flu_cfg$knot_step, end_t, flu_cfg$gamma_grid, flu_cfg$n_fv)
gamma_grid <- fv$gamma_grid

rrt <- retro_rt(season); rt_retro_at_cutoff <- rrt$Rt_mean[rrt$date == WEEK_END]
theta_init <- fit_main$theta; rows <- list()
for (gi in seq_along(gamma_grid)) {
  fg <- tryCatch(ConvRt:::.gi_solve_taper(d_sc, LAM_USED, lam_taper = gamma_grid[gi], P_taper = P_taper,
                                 tail = "linear", theta_init = theta_init), error = function(e) NULL)
  if (is.null(fg) || !all(is.finite(fg$theta))) next
  rt <- gi_extract_rt(fg, d_sc, lam = LAM_USED, lam_taper = gamma_grid[gi],
                      P_taper = P_taper, level = 0.90, overdispersion = TRUE)
  rt$date <- win$date[rt$day + 1L]; rt <- rt[rt$date <= WEEK_END, ]
  rt$gamma <- gamma_grid[gi]; rt$gamma_index <- gi; rows[[gi]] <- rt; theta_init <- fg$theta
}
rt_per_gamma <- bind_rows(rows)
rt_at_cutoff_by_gamma <- rt_per_gamma %>% filter(date == WEEK_END) %>%
  transmute(gamma_index, gamma, Rt_at_cutoff = Rt_mean, err_to_retro = abs(Rt_mean - rt_retro_at_cutoff))
idx_retro <- rt_at_cutoff_by_gamma$gamma_index[which.min(rt_at_cutoff_by_gamma$err_to_retro)]

out_dir <- file.path(flu_results_dir, "girt", season, "real_time", "gamma", format(WEEK_END))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(season = season, week_end_date = WEEK_END, lam_used = LAM_USED,
             retro_rt_at_cutoff = rt_retro_at_cutoff, gamma_grid = gamma_grid,
             fv_mean = fv$fv_mean, fv_se = fv$fv_se, best_gamma_min = fv$best_gamma_min,
             best_gamma_1se = fv$best_gamma_1se, closest_to_retro_gamma = gamma_grid[idx_retro],
             rt_per_gamma = rt_per_gamma, rt_at_cutoff_by_gamma = rt_at_cutoff_by_gamma,
             retro_rt_curve = rrt, lambda_mode = MODE),
        file.path(out_dir, sprintf("data%s.rds", suffix)))
cat(sprintf("wrote %s  (lam=%.3e gamma_min=%.3e gamma_1se=%.3e)\n",
            file.path(out_dir, sprintf("data%s.rds", suffix)), LAM_USED, fv$best_gamma_min, fv$best_gamma_1se))
