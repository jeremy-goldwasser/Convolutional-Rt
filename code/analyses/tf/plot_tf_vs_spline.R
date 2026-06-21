# code/analyses/tf/plot_tf_vs_spline.R
#   -> figures/extra_methods/tf/tf_vs_spline.pdf  (Season 1 2022-23)
# Reuses the flu analysis setup (flu_cfg, load_flu_hosps, flu_seasons, ...) by
# sourcing ../flu/_common.R.
# girt cubic trend-filter vs the cached girt natural-spline retro on the real
# flu Season 1.  Curves shown:
#   - MechRt (spline)              -- cached girt retro
#   - MechRt-TF (1se)              -- CV-1se lambda
#   - MechRt-TF (oracle vs spline) -- lambda in the grid minimising
#                                     |TF - spline| MAE in the season window
# (TF-min dropped per user preference; 1se is the default rule.)
# Season windows match plot_retro_both_seasons.R:
#   s1: 2022-09-01 .. 2023-02-01,  s2: 2023-09-01 .. 2024-04-12.
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(tibble); library(patchwork)
})
sd_ <- dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))
if (!length(sd_) || is.na(sd_)) sd_ <- "code/analyses/tf"
source(file.path(sd_, "..", "flu", "_common.R"))

tf_cfg <- list(k_diff = 4L, lam_grid = 10^seq(0, 7, length.out = 25),
               nfold = 5L, seed = 1L, knot_step = 5L)
panel_windows <- list(s1 = c(as.Date("2022-09-01"), as.Date("2023-02-01")),
                      s2 = c(as.Date("2023-09-01"), as.Date("2024-04-12")))

pi_EH <- ce_gamma_pmf(flu_cfg$mean_hosp,       flu_cfg$sd_hosp)
g_flu <- flu_gi()
hosps_all <- load_flu_hosps()

# Build TF design + CV + fit at every lambda on the grid.
fit_tf_full_grid_season <- function(season) {
  sw  <- flu_seasons[[season]]
  win <- hosps_all[hosps_all$date >= sw$beta - 30L & hosps_all$date <= sw$eos, ]
  win$day <- seq_len(nrow(win)) - 1L
  X_hat <- gi_deconvolve_exposures(win$hosps_1d, pi_EY = pi_EH,
              severity_rate = flu_cfg$severity_rate,
              link = flu_cfg$deconv_link, burn_in = 30L)$X_hat
  d_tf <- suppressWarnings(build_gi_design_tf(
    E_star = X_hat, g = g_flu, pi_EY = pi_EH, max_time = nrow(win) - 1L,
    inf_inc = win$hosps_1d, knot_step = tf_cfg$knot_step,
    first_rt_index = sw$beta, dates = win$date, dow_dates = win$date))
  d_tf <- gi_enforce_likelihood_start(d_tf, win$day[which(d_tf$valid_mask)],
                                      win$day[win$date == sw$beta + 21L])
  d_tf <- gi_apply_severity_to_design(d_tf, rep(flu_cfg$severity_rate, nrow(win)))

  cat(sprintf("\n[%s] DoW + CV + fit-all-lambdas (k=%d, %d lams)\n",
              season, tf_cfg$k_diff, length(tf_cfg$lam_grid)))
  t0 <- Sys.time()
  dow <- .gi_tf_estimate_dow(d_tf, stats::median(tf_cfg$lam_grid), tf_cfg$k_diff)
  d_t <- as.numeric(dow[d_tf$dow_valid])
  d_tf$Z_full <- d_tf$Z_full * d_t; d_tf$off_valid <- d_tf$off_valid * d_t
  d_tf$use_dow <- FALSE; d_tf$dow_valid <- NULL
  sel  <- gi_select_lambda_cv_tf(d_tf, tf_cfg$lam_grid, k_diff = tf_cfg$k_diff,
                                 nfold = tf_cfg$nfold, seed = tf_cfg$seed)
  prob <- gi_build_tf_problem(d_tf, k_diff = tf_cfg$k_diff, active_cols = sel$active_cols)
  curves <- lapply(tf_cfg$lam_grid, function(lam) gi_solve_tf(prob, lam)$theta)
  cat(sprintf("  done in %.1fs.  lam_1se=%.2e\n",
              as.numeric(Sys.time() - t0, units = "secs"), sel$best_lam_1se))
  list(curves = curves, sel = sel, dates = win$date, season = season)
}

oracle_lam_vs_spline <- function(res, spline_rt, win) {
  ix <- which(res$dates >= win[1] & res$dates <= win[2])
  spline_in_win <- spline_rt[match(res$dates[ix], spline_rt$date), "Rt_mean", drop = TRUE]
  mae <- vapply(res$curves, function(th) {
    th_w <- th[ix]; mask <- is.finite(th_w) & is.finite(spline_in_win)
    if (sum(mask) < 5L) return(NA_real_); mean(abs(th_w[mask] - spline_in_win[mask]))
  }, numeric(1))
  list(lam = tf_cfg$lam_grid[which.min(mae)], mae = min(mae, na.rm = TRUE))
}

load_spline <- function(season) {
  f <- list.files(file.path(flu_results_dir, "girt", season, "retrospective", "results"),
                  full.names = TRUE)[1]
  rt <- readRDS(f)$rt_df
  tibble(date = as.Date(rt$date), Rt_mean = rt$Rt_mean)
}

assemble <- function(season) {
  res <- fit_tf_full_grid_season(season); win <- panel_windows[[season]]
  spline_rt <- load_spline(season)
  lam_1se <- res$sel$best_lam_1se; idx_1se <- which(res$sel$lam_grid == lam_1se)[1]
  oracle  <- oracle_lam_vs_spline(res, spline_rt, win)
  idx_orc <- which(tf_cfg$lam_grid == oracle$lam)[1]
  cat(sprintf("  [%s] lam_1se=%.2e  oracle_lam=%.2e  oracle_MAE=%.4f\n",
              season, lam_1se, oracle$lam, oracle$mae))

  curve_df <- function(th, method) tibble(
      date = as.Date(res$dates), Rt_mean = as.numeric(th), method = method,
      season = season) |> filter(date >= win[1], date <= win[2])
  spl_df <- spline_rt |> filter(date >= win[1], date <= win[2]) |>
    mutate(method = "MechRt (spline)", season = season)

  list(rt = bind_rows(spl_df,
         curve_df(res$curves[[idx_1se]], "MechRt-TF (1se)"),
         curve_df(res$curves[[idx_orc]], "MechRt-TF (oracle vs spline)")),
       summary = tibble(season = season, lam_1se = lam_1se,
                        lam_oracle = oracle$lam, oracle_MAE_vs_spline = oracle$mae))
}

s1 <- assemble("s1")
cat("\n--- Summary ---\n"); print(as.data.frame(s1$summary))

method_colors <- c("MechRt (spline)" = "#0072B2", "MechRt-TF (1se)" = "#D55E00",
                   "MechRt-TF (oracle vs spline)" = "#009E73")
method_ltype  <- c("MechRt (spline)" = "solid", "MechRt-TF (1se)" = "dashed",
                   "MechRt-TF (oracle vs spline)" = "dotted")

make_panel <- function(df, title, base_size = 12, line_width = 0.8, hline = TRUE,
                       labels = ggplot2::waiver()) {
  p <- ggplot(df, aes(date, Rt_mean, color = method, linetype = method))
  if (hline) p <- p + geom_hline(yintercept = 1, linetype = "dotted", color = "grey55")
  p +
    geom_line(linewidth = line_width) +
    scale_color_manual(values = method_colors, name = NULL, labels = labels) +
    scale_linetype_manual(values = method_ltype, name = NULL, labels = labels) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b\n%Y") +
    labs(title = title, x = NULL, y = expression(R[t])) +
    theme_bw(base_size = base_size) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = rel(0.9)),
          legend.key.width = unit(2, "lines"),
          plot.title = element_text(hjust = 0.5, face = "bold"))
}

# Polished standalone Season 1 figure: cleaner title, no Rt=1 reference line.
out_dir <- file.path(figures_dir, "extra_methods", "tf"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ttl <- expression(paste("girt: trend filtering vs. spline ",
                        R[t], " estimates (NHSN flu data)"))
s1_labels <- c("MechRt (spline)" = "Spline",
               "MechRt-TF (1se)" = "Trend filtering (1se)",
               "MechRt-TF (oracle vs spline)" = "Trend filtering (oracle)")
ggsave(file.path(out_dir, "tf_vs_spline.pdf"),
       make_panel(s1$rt, ttl, base_size = 13, line_width = 1.0, hline = FALSE,
                  labels = s1_labels),
       width = 8, height = 3.6)
cat(sprintf("Wrote %s\n", file.path(out_dir, "tf_vs_spline.pdf")))
