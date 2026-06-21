# code/analyses/rtestim/regen_baselines.R
#
# Regenerate the rtestim-SEIR baseline Rt fits (estimateR, EpiEstim, rtestim;
# EpiLPS appended below) LOCALLY, using the TRUE SEIR-implied generation
# interval Gamma(theor_mean, theor_sd) from seir_params.rds, instead of the
# rounded target Gamma(8.4, 3.8) used in the original cluster runs.
#
# The estimateR / EpiEstim / rtestim calls are copied verbatim from the cluster
# generator (Rt_scripts/.../rtestim_paper_retro_benchmark.r); the ONLY change is
# that delays$si_w is built at the true GI. estimateR (bootstrap) and rtestim
# (CV folds) are stochastic and were unseeded on the cluster, so we seed them
# here for reproducibility.
#
# Env vars (for validation):
#   GI_MEAN, GI_SD  -- override the GI (default: true theor values)
#   OUT_DIR         -- where to write fit_<method>_<tag>.rds (default: canonical
#                      results/sim/rtestim). Set to a staging dir to validate
#                      without overwriting the cached fits.

sd_dir <- dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))
if (!length(sd_dir) || is.na(sd_dir)) sd_dir <- "code/analyses/rtestim"
source(file.path(sd_dir, "_common.R"))   # discrete_gamma, auto_k_gamma, sei_par, rt_*_dir, rt_cfg
suppressPackageStartupMessages({ library(EpiEstim); library(estimateR); library(rtestim); library(EpiLPS); library(dplyr); library(tibble) })

BOOT_SEED <- 1L
sp  <- sei_par[sei_par$si_type == "fake", ]
gi_mean <- suppressWarnings(as.numeric(Sys.getenv("GI_MEAN", ""))); if (!is.finite(gi_mean)) gi_mean <- sp$theor_mean
gi_sd   <- suppressWarnings(as.numeric(Sys.getenv("GI_SD",   ""))); if (!is.finite(gi_sd))   gi_sd   <- sp$theor_sd
out_dir <- Sys.getenv("OUT_DIR", rt_results_dir); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("GI = Gamma(mean=%.6f, sd=%.6f)\nOUT_DIR = %s\n\n", gi_mean, gi_sd, out_dir))

# Delays exactly as the cluster benchmark builds them, but GI = (gi_mean, gi_sd).
p <- sp
delays <- list(
  pi_EY   = discrete_gamma(p$mean_EY, p$sd_EY, k = auto_k_gamma(p$mean_EY, p$sd_EY)),
  si_w    = discrete_gamma(gi_mean,   gi_sd,   k = auto_k_gamma(gi_mean,   gi_sd)),
  mean_EY = p$mean_EY, sd_EY = p$sd_EY
)

tag_for <- function(sc) sprintf("fake__%s", gsub(":? ", "_", sc))
sim     <- readRDS(file.path(rt_data_dir, "sim_combined.rds")); sim$date <- rt_cfg$date0 + (sim$time - 1L)
scenarios <- unique(sim$scenario)

# ---- METHOD: estimateR (LOESS + R-L deconv + EpiEstim sliding window) --------
fit_estimater <- function(dates, y, sc, delays) {
  case_delay <- list(name = "gamma",
                     shape = (delays$mean_EY / delays$sd_EY)^2,
                     scale = delays$sd_EY^2 / delays$mean_EY)
  si_distr   <- c(0, delays$si_w)
  set.seed(BOOT_SEED)
  boot <- get_block_bootstrapped_estimate(
    incidence_data = y, N_bootstrap_replicates = 50L,
    ref_date = dates[1], time_step = "day",
    smoothing_method = "LOESS", deconvolution_method = "Richardson-Lucy delay distribution",
    estimation_method = "EpiEstim sliding window", delay = case_delay,
    estimation_window = 3L, method = "non_parametric_si", si_distr = si_distr, mean_Re_prior = 1)
  rt <- as_tibble(boot) |> transmute(method = "estimateR", scenario = sc, si_type = "fake",
            date = as.Date(date), Rt_mean = Re_estimate, Rt_lo = CI_down_Re_estimate, Rt_hi = CI_up_Re_estimate)
  list(rt = rt, deconvolved = NULL)
}

# ---- METHOD: EpiEstim (raw reports, 7d sliding window, paper SI) -------------
fit_epiestim <- function(dates, y, sc, delays) {
  n <- length(y); t_start <- 2L:(n - 6L); t_end <- t_start + 6L
  si_distr <- c(0, delays$si_w)
  res <- suppressMessages(estimate_R(incid = y, method = "non_parametric_si",
    config = make_config(list(si_distr = si_distr, t_start = t_start, t_end = t_end,
                              mean_prior = 1, std_prior = 5))))
  tibble(method = "EpiEstim", scenario = sc, si_type = "fake", date = dates[t_end],
         Rt_mean = res$R$`Mean(R)`, Rt_lo = res$R$`Quantile.0.025(R)`, Rt_hi = res$R$`Quantile.0.975(R)`)
}

# ---- METHOD: rtestim (trend filter korder=3, CV-1se, paper SI) ---------------
fit_rtestim <- function(dates, y, sc, delays) {
  first_pos <- which(cumsum(y) > 0)[1]; if (is.na(first_pos)) first_pos <- 1L
  idx <- seq(first_pos, length(y)); y_fit <- as.integer(pmax(y[idx], 0L)); d_fit <- dates[idx]
  lam_grid_rt <- 10^seq(5, -2, length.out = 30)
  si_for_rt   <- c(0, delays$si_w)
  set.seed(BOOT_SEED)
  cv <- cv_estimate_rt(observed_counts = y_fit, korder = 3L, delay_distn = si_for_rt,
                       nfold = 3L, lambda = lam_grid_rt, maxiter = 1e7L)
  cb <- confband(cv, lambda = cv$lambda.1se, level = rt_cfg$level, type = "Rt")
  rt <- tibble(method = "rtestim", scenario = sc, si_type = "fake", date = d_fit,
               Rt_mean = as.numeric(cb$fit), Rt_lo = as.numeric(cb$`2.5%`), Rt_hi = as.numeric(cb$`97.5%`))
  list(rt = rt, lambda_1se = cv$lambda.1se, lambda_min = cv$lambda.min)
}

# ---- METHOD: EpiLPS (Bayesian P-splines, K=30, default priors, paper SI) -----
# NOTE: EpiLPS::estimR takes the SI starting at lag 1 (no leading 0), unlike the
# Cori-convention c(0, si_w) used by the other three. Verified to reproduce the
# cached fits exactly at the target GI.
fit_epilps <- function(dates, y, sc, delays) {
  priors <- Rmodelpriors(listcontrol = list(a_delta = 10, b_delta = 10, phi = 2,
                                             a_rho = 1e-4, b_rho = 1e-4))
  t0 <- Sys.time()
  res <- EpiLPS::estimR(incidence = as.integer(pmax(0L, round(y))), si = delays$si_w,
                        K = 30L, dates = dates, priors = priors)
  t1 <- Sys.time()
  rt <- tibble(method = "EpiLPS", scenario = sc, si_type = "fake",
               date = as.Date(res$RLPS$Time), Rt_mean = res$RLPS$R,
               Rt_lo = res$RLPS$Rq0.025, Rt_hi = res$RLPS$Rq0.975)
  list(rt = rt, elapsed_sec = as.numeric(difftime(t1, t0, units = "secs")))
}

for (sc in scenarios) {
  d <- sim[sim$scenario == sc & sim$si_type == "fake", ]; d <- d[order(d$time), ]
  y <- as.numeric(d$seir_reports); dates <- d$date; tag <- tag_for(sc)
  cat(sprintf("[%s]\n", sc))
  saveRDS(fit_estimater(dates, y, sc, delays), file.path(out_dir, sprintf("fit_estimater_%s.rds", tag)))
  cat("  estimateR done\n")
  saveRDS(fit_epiestim(dates, y, sc, delays),  file.path(out_dir, sprintf("fit_epiestim_%s.rds", tag)))
  cat("  EpiEstim done\n")
  saveRDS(fit_rtestim(dates, y, sc, delays),   file.path(out_dir, sprintf("fit_rtestim_%s.rds", tag)))
  cat("  rtestim done\n")
  saveRDS(fit_epilps(dates, y, sc, delays),    file.path(out_dir, sprintf("fit_epilps_%s.rds", tag)))
  cat("  EpiLPS done\n")
}
cat("\nDONE: estimateR/EpiEstim/rtestim/EpiLPS written to", out_dir, "\n")
