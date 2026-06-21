# code/analyses/flu/regen.R  --  regenerate all girt-log real-flu fits.
#
#   Rscript regen.R <season> <mode>
#     season: s1 | s2
#     mode  : retro | weekly | daily | all   (default all)
#
# Writes to results/real/flu/girt/<season>/{retrospective,real_time}/...
# Idempotent for daily/weekly (skips existing).  g = discrete Gamma GI, LOG deconv.

sd <- dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))
if (!length(sd) || is.na(sd)) sd <- "code/analyses/flu"
source(file.path(sd, "_common.R")); source(file.path(sd, "_realtime.R"))
suppressPackageStartupMessages(library(tibble))
args <- commandArgs(trailingOnly = TRUE)
season <- if (length(args) >= 1) args[1] else "s2"
mode   <- if (length(args) >= 2) args[2] else "all"
stopifnot(season %in% c("s1", "s2"))

pEH    <- gi_discrete_gamma_delay(flu_cfg$mean_hosp, flu_cfg$sd_hosp)
g      <- flu_gi()
hosps  <- load_flu_hosps()
beta   <- flu_seasons[[season]]$beta
eos    <- flu_seasons[[season]]$eos
gtree  <- file.path(flu_results_dir, "girt", season)
wk     <- function(a, b) { d <- seq(as.Date(a), as.Date(b), by = 7L); d[d %in% hosps$date] }
weekly_W <- if (season == "s1") wk("2022-10-01", "2023-02-25") else wk("2023-10-01", "2024-04-21")
daily_lo <- if (season == "s1") as.Date("2022-08-01") else as.Date("2023-08-01")
daily_W  <- { d <- seq(daily_lo, eos, by = "1 day"); d[d %in% hosps$date] }

retro_lam_path <- file.path(gtree, "retrospective", "results")
get_retro_lam <- function() {
  f <- list.files(retro_lam_path, pattern = "natural_linear_notaper_cv_1se.rds$", full.names = TRUE)
  if (length(f)) readRDS(f[1])$lam_1se else NA_real_
}

if (mode %in% c("retro", "all")) {
  cat(sprintf("[%s] retro EOS=%s ...\n", season, eos))
  dir.create(retro_lam_path, recursive = TRUE, showWarnings = FALSE)
  r <- fit_flu_retro_girt(hosps, eos, beta, beta, pEH, g)
  saveRDS(list(meta = list(season = season, eos = eos, deconv_link = "log"),
               rt_df = r$rt_df, lam_1se = r$lam_1se, dow_effects = r$dow_effects, hosps = r$hosps),
          file.path(retro_lam_path, sprintf("combo_00001__%s__natural_linear_notaper_cv_1se.rds", format(eos, "%Y-%m-%d"))))
  cat(sprintf("  retro lam_1se=%.4g done\n", r$lam_1se))
}

lam_retro <- get_retro_lam()

run_vintages <- function(W_vec, out_dir, namer, slim = FALSE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  for (i in seq_along(W_vec)) {
    W <- W_vec[i]; out <- file.path(out_dir, namer(i, W))
    if (file.exists(out)) next
    gf <- tryCatch(fit_flu_combo_girt(hosps, W, beta, beta, pEH, g,
                     use_taper = TRUE, lam_retro_1se = lam_retro),
                   error = function(e) { cat(sprintf("  [%s] ERR %s\n", W, conditionMessage(e))); NULL })
    if (is.null(gf)) next
    if (slim) gf <- list(meta = gf$meta, rt_df = gf$rt_df)
    saveRDS(gf, out)
    if (i %% 20 == 0 || i == length(W_vec)) cat(sprintf("  [%s %s] %d/%d\n", season, basename(out_dir), i, length(W_vec)))
  }
}

if (mode %in% c("weekly", "all")) {
  cat(sprintf("[%s] weekly: %d vintages\n", season, length(weekly_W)))
  run_vintages(weekly_W, file.path(gtree, "real_time", "results"),
               function(i, W) sprintf("combo_%05d__%s__natural_linear_taper_cv_1se.rds", i, format(W, "%Y-%m-%d")))
}
if (mode %in% c("constant", "all")) {
  cat(sprintf("[%s] regression-constant weekly: %d vintages\n", season, length(weekly_W)))
  out_dir <- file.path(gtree, "real_time", "results"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  for (i in seq_along(weekly_W)) {
    W <- weekly_W[i]
    out <- file.path(out_dir, sprintf("combo_%05d__%s__regression_constant_notaper_cv_1se.rds", i, format(W, "%Y-%m-%d")))
    if (file.exists(out)) next
    gf <- tryCatch(fit_flu_constant_girt(hosps, W, beta, beta, g, pEH$pmf),
                   error = function(e) { cat(sprintf("  [%s] ERR %s\n", W, conditionMessage(e))); NULL })
    if (!is.null(gf)) saveRDS(gf, out)
    if (i %% 10 == 0 || i == length(weekly_W)) cat(sprintf("  [%s constant] %d/%d\n", season, i, length(weekly_W)))
  }
}
if (mode %in% c("daily", "all")) {
  cat(sprintf("[%s] daily: %d vintages\n", season, length(daily_W)))
  run_vintages(daily_W, file.path(gtree, "real_time", "daily_fits"),
               function(i, W) sprintf("daily_%s.rds", format(W, "%Y-%m-%d")), slim = TRUE)
}
cat(sprintf("[%s] regen mode=%s complete.\n", season, mode))
