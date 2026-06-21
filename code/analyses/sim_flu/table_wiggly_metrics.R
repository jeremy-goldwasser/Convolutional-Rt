# code/analyses/sim_flu/table_wiggly_metrics.R
#   -> figures/sim/flu/wiggly/flu3_wiggly_metrics_identity.csv
#   (also writes table/{retro,realtime}_metrics.csv as intermediates)
# Headline metrics table for the wiggly scenario: MAE, WIS@90, CE
# (mean |emp - alpha| across alpha_grid).  Eval window: 2022-10-01 to 2023-01-31.
# All inputs are the cached canonical fit objects (post-conformal CIs already on
# the weekly real-time files); Gaussian-rescaling of (Rt_lo, Rt_hi) gives the
# per-alpha sigma for the CE rollup.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(tibble) })
sd_ <- dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))
if (!length(sd_) || is.na(sd_)) sd_ <- "code/analyses/sim_flu"
source(file.path(sd_, "_common.R"))
p <- sim_flu_paths("wiggly")
dir.create(p$table_dir, recursive = TRUE, showWarnings = FALSE)

eval_lo <- as.Date("2022-10-01"); eval_hi <- as.Date("2023-01-31")
alpha_grid <- c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95); z95 <- qnorm(0.975); z90 <- qnorm(0.95)

sim <- readRDS(p$sim_rds); sim$date <- as.Date(sim$date)
truth <- tibble(date = sim$date, Rt_truth = sim$rt_raw)

sigma_from <- function(lo, hi, level = 0.95) (hi - lo) / (2 * qnorm(1 - (1 - level)/2))
empirical_coverage <- function(rt) {
  rt <- rt[is.finite(rt$Rt_truth) & is.finite(rt$Rt_mean) & is.finite(rt$sigma), ]
  if (!nrow(rt)) return(rep(NA_real_, length(alpha_grid)))
  vapply(alpha_grid, function(a) {
    z <- qnorm(1 - (1 - a)/2)
    mean(rt$Rt_truth >= rt$Rt_mean - z*rt$sigma & rt$Rt_truth <= rt$Rt_mean + z*rt$sigma)
  }, numeric(1))
}
wis90 <- function(y, m, lo, hi) (1/1.5) * (0.5*abs(y-m) + 0.05*(hi-lo) + pmax(0, lo-y) + pmax(0, y-hi))

eval_band <- function(date, Rt_mean, Rt_lo, Rt_hi, level = 0.95, has_ci = TRUE) {
  ev <- tibble(date = as.Date(date), Rt_mean = Rt_mean, Rt_lo = Rt_lo, Rt_hi = Rt_hi) |>
    inner_join(truth, by = "date") |> filter(date >= eval_lo, date <= eval_hi)
  if (!nrow(ev)) return(tibble(n = 0L, MAE = NA_real_, WIS = NA_real_, CE = NA_real_))
  MAE <- mean(abs(ev$Rt_mean - ev$Rt_truth))
  if (!has_ci || all(!is.finite(ev$Rt_lo))) return(tibble(n = nrow(ev), MAE = MAE, WIS = NA_real_, CE = NA_real_))
  ev$sigma <- sigma_from(ev$Rt_lo, ev$Rt_hi, level)
  lo90 <- ev$Rt_mean - z90 * ev$sigma; hi90 <- ev$Rt_mean + z90 * ev$sigma
  WIS <- mean(wis90(ev$Rt_truth, ev$Rt_mean, lo90, hi90))
  CE  <- mean(abs(empirical_coverage(ev) - alpha_grid))
  tibble(n = nrow(ev), MAE = MAE, WIS = WIS, CE = CE)
}

# ================= RETRO =================
rt_all <- readRDS(file.path(p$retro_dir, "retro_sim_rt_estimates.rds"))
rt_all$date <- as.Date(rt_all$date)
retro_rows <- list()
for (m in c("MechRt", "estimateR", "EpiEstim", "rtestim", "EpiLPS")) {
  sub <- rt_all |> filter(as.character(method) == m, variant == "raw")
  if (!nrow(sub)) next
  r <- eval_band(sub$date, sub$Rt_mean, sub$Rt_lo, sub$Rt_hi, level = 0.95,
                 has_ci = any(is.finite(sub$Rt_lo)))
  retro_rows[[length(retro_rows)+1L]] <- cbind(tibble(fit = "retrospective", Method = m), r)
}
# EpiNow2 retro from epinow2_retro/
ep_files <- list.files(p$ep_retro, "__raw_gp__summary\\.rds$", full.names = TRUE)
if (length(ep_files)) {
  ep <- readRDS(ep_files[1])$rt_df
  rt_lo <- if (!is.null(ep$Rt_sd)) ep$Rt_mean - z95*ep$Rt_sd else NA_real_
  rt_hi <- if (!is.null(ep$Rt_sd)) ep$Rt_mean + z95*ep$Rt_sd else NA_real_
  r <- eval_band(ep$date, ep$Rt_mean, rt_lo, rt_hi, level = 0.95, has_ci = !is.null(ep$Rt_sd))
  retro_rows[[length(retro_rows)+1L]] <- cbind(tibble(fit = "retrospective", Method = "EpiNow2"), r)
}
retro_df <- do.call(rbind, retro_rows)
write_csv(retro_df, file.path(p$table_dir, "retro_metrics.csv"))

# ================= REALTIME (trailing 7 days, d_to_edge <= 6) =================
D_MAX <- 6L
man <- read.csv(p$manifest); man$week_end_date <- as.Date(man$week_end_date)
load_mech <- do.call(rbind, lapply(seq_len(nrow(man)), function(i) {
  f <- file.path(p$weekly_dir, man$result_filename[i]); if (!file.exists(f)) return(NULL)
  r <- readRDS(f); rt <- r$rt_df
  data.frame(date = as.Date(rt$date), week_end = man$week_end_date[i],
             Rt_mean = rt$Rt_mean, Rt_lo = rt$Rt_lo, Rt_hi = rt$Rt_hi, stringsAsFactors = FALSE)
})) |> mutate(d_to_edge = as.integer(week_end - date)) |> filter(d_to_edge >= 0L, d_to_edge <= D_MAX)
load_other <- function(method) {
  path <- file.path(p$other_rt, sprintf("%s__raw.rds", method)); if (!file.exists(path)) return(NULL)
  d <- readRDS(path)$rt |> mutate(date = as.Date(date), week_end = as.Date(week_end),
                                  d_to_edge = as.integer(week_end - date)) |>
    filter(d_to_edge >= 0L, d_to_edge <= D_MAX)
  # ensure Rt_lo/Rt_hi columns exist (some cache only Rt_mean)
  if (!"Rt_lo" %in% names(d)) d$Rt_lo <- NA_real_
  if (!"Rt_hi" %in% names(d)) d$Rt_hi <- NA_real_
  d
}
ep_rt_files <- list.files(p$ep_rt, "__raw_gp__summary\\.rds$", full.names = TRUE)
load_ep <- do.call(rbind, lapply(ep_rt_files, function(f) {
  r <- readRDS(f); rt <- r$rt_df
  Rt_lo <- if (!is.null(rt$Rt_sd)) rt$Rt_mean - z95*rt$Rt_sd else NA_real_
  Rt_hi <- if (!is.null(rt$Rt_sd)) rt$Rt_mean + z95*rt$Rt_sd else NA_real_
  data.frame(date = as.Date(rt$date), week_end = as.Date(r$meta$week_end_date),
             Rt_mean = rt$Rt_mean, Rt_lo = Rt_lo, Rt_hi = Rt_hi, stringsAsFactors = FALSE)
})) |> mutate(d_to_edge = as.integer(week_end - date)) |> filter(d_to_edge >= 0L, d_to_edge <= D_MAX)

realtime_rows <- list()
realtime_rows[[length(realtime_rows)+1L]] <- cbind(tibble(fit = "realtime", Method = "MechRt"),
  eval_band(load_mech$date, load_mech$Rt_mean, load_mech$Rt_lo, load_mech$Rt_hi, level = 0.95, has_ci = TRUE))
realtime_rows[[length(realtime_rows)+1L]] <- cbind(tibble(fit = "realtime", Method = "EpiNow2"),
  eval_band(load_ep$date, load_ep$Rt_mean, load_ep$Rt_lo, load_ep$Rt_hi, level = 0.95,
            has_ci = any(is.finite(load_ep$Rt_lo))))
for (m in c("estimateR", "EpiEstim", "rtestim", "EpiLPS")) {
  d <- load_other(m); if (is.null(d) || !nrow(d)) next
  realtime_rows[[length(realtime_rows)+1L]] <- cbind(tibble(fit = "realtime", Method = m),
    eval_band(d$date, d$Rt_mean, d$Rt_lo, d$Rt_hi, level = 0.95,
              has_ci = any(is.finite(d$Rt_lo))))
}
realtime_df <- do.call(rbind, realtime_rows)
write_csv(realtime_df, file.path(p$table_dir, "realtime_metrics.csv"))

# ================= FINAL combined CSV (figures/) =================
final_df <- bind_rows(retro_df, realtime_df) |> mutate(Runtime_sec = NA_real_) |>
  select(fit, Method, n, MAE, WIS, CE, Runtime_sec)
out <- file.path(p$fig_dir, "flu3_wiggly_metrics_identity.csv")
dir.create(p$fig_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(final_df, out); cat("wrote", out, "\n")
print(as.data.frame(final_df), digits = 3, row.names = FALSE)
