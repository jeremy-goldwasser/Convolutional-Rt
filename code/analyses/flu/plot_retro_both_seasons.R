# code/analyses/flu/plot_retro_both_seasons.R
#
# Both-seasons RETROSPECTIVE figures, regenerated with girt (log deconv):
#   figures/real/flu/both_seasons/comparison_rt_both_seasons_with_hosps_s2eos.pdf
#   figures/real/flu/both_seasons/comparison_rt_both_seasons_mechrt_simband_vs_pointwise.pdf
# MechRt Rt refit live via girt (fit_flu_retro_girt); EpiNow2 read from canonical tree.

suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(patchwork) })
sd <- dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))
if (!length(sd) || is.na(sd)) sd <- "code/analyses/flu"
source(file.path(sd, "_common.R")); source(file.path(sd, "_realtime.R"))

col_ours <- "#2D6CA2"; col_ep <- "#5D3A1E"
season_labels <- c(s1 = "2022/23", s2 = "2023/24")
rt_scale_y <- function() scale_y_continuous(breaks = seq(0, 5, by = 0.1))
panel_windows <- list(s1 = c(as.Date("2022-09-01"), as.Date("2023-02-01")),
                      s2 = c(as.Date("2023-09-01"), as.Date("2024-04-12")))
shared_x_scale <- function(w) scale_x_date(date_breaks = "1 month", date_labels = "%b\n%Y",
                                           limits = w, expand = expansion(mult = c(0.005, 0.005)))
band_line_key <- function(a) function(data, params, size)
  grid::grobTree(grid::rectGrob(gp = grid::gpar(col = NA, fill = scales::alpha(data$colour, a))),
                 ggplot2::draw_key_path(data, params, size))

pEH    <- gi_discrete_gamma_delay(flu_cfg$mean_hosp, flu_cfg$sd_hosp)
g      <- flu_gi()
hosps  <- load_flu_hosps()

load_epinow2 <- function(season) {
  d <- file.path(flu_results_dir, season, "retrospective", "epinow2", "results")
  f <- list.files(d, pattern = "__gp__summary\\.rds$", full.names = TRUE)
  if (!length(f)) stop("no EpiNow2 retro summary for ", season); readRDS(f[[1]])$rt_df
}

run_season_girt <- function(season) {
  b <- flu_seasons[[season]]$beta; e <- flu_seasons[[season]]$eos
  fit <- fit_flu_retro_girt(hosps, e, b, b, pEH, g, level = 0.90)
  rt <- fit$rt_df; rt <- rt[rt$date >= b & rt$date <= e, ]
  list(season = season, rt_ours = rt, rt_ep = load_epinow2(season), hosps = fit$hosps)
}
all_res <- list(s1 = run_season_girt("s1"), s2 = run_season_girt("s2"))

# ---- Figure 1: pointwise Rt (MechRt vs EpiNow2) + daily hosps, 2x2 ----------
plot_cmp <- function(res, w) {
  z <- qnorm(0.975) / qnorm(0.95)
  o <- res$rt_ours[res$rt_ours$date >= w[1] & res$rt_ours$date <= w[2], c("date","Rt_mean","Rt_lo","Rt_hi")]
  o$Rt_lo <- o$Rt_mean - z * (o$Rt_mean - o$Rt_lo); o$Rt_hi <- o$Rt_mean + z * (o$Rt_hi - o$Rt_mean)
  ep <- res$rt_ep[res$rt_ep$date >= w[1] & res$rt_ep$date <= w[2], c("date","Rt_mean","Rt_sd")]
  ep$lo <- ep$Rt_mean - qnorm(0.975) * ep$Rt_sd; ep$hi <- ep$Rt_mean + qnorm(0.975) * ep$Rt_sd
  ggplot() +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_ribbon(data = ep, aes(date, ymin = lo, ymax = hi), fill = col_ep, alpha = 0.18) +
    geom_ribbon(data = o, aes(date, ymin = Rt_lo, ymax = Rt_hi), fill = col_ours, alpha = 0.30) +
    geom_line(data = ep, aes(date, Rt_mean, colour = "EpiNow2"), linewidth = 0.9, key_glyph = band_line_key(0.18)) +
    geom_line(data = o, aes(date, Rt_mean, colour = "MechRt"), linewidth = 0.9, key_glyph = band_line_key(0.30)) +
    scale_colour_manual(values = c("MechRt" = col_ours, "EpiNow2" = col_ep), name = NULL, breaks = c("MechRt","EpiNow2")) +
    rt_scale_y() + labs(x = NULL, y = expression(R[t]), title = sprintf("%s flu season", season_labels[[res$season]])) +
    theme_minimal(base_size = 16) +
    theme(axis.text = element_text(size = 16), axis.title.y = element_text(size = 18),
          legend.position = "bottom", plot.title = element_text(hjust = 0.5, size = 20)) +
    shared_x_scale(w) + guides(colour = guide_legend(override.aes = list(linewidth = 1.2)))
}
plot_hosps <- function(res, w) {
  d <- res$hosps[res$hosps$date >= w[1] & res$hosps$date <= w[2], ]
  ggplot(d, aes(date, hosps_1d)) + geom_col(fill = "grey40", width = 0.9) + shared_x_scale(w) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + labs(x = NULL, y = "Daily hosps") +
    theme_minimal(base_size = 16) +
    theme(axis.text = element_text(size = 16), axis.title.y = element_text(size = 18), panel.grid.minor = element_blank())
}
strip_x <- function(p) p + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
top1 <- strip_x(plot_cmp(all_res$s1, panel_windows$s1)); top2 <- strip_x(plot_cmp(all_res$s2, panel_windows$s2))
fig1 <- (top1 + top2) / (plot_hosps(all_res$s1, panel_windows$s1) + plot_hosps(all_res$s2, panel_windows$s2)) +
  plot_layout(heights = c(2.4, 1), guides = "collect") &
  theme(legend.position = "bottom", legend.text = element_text(size = 20),
        legend.key.width = unit(1.4, "cm"))
out1 <- file.path(flu_figures_dir, "both_seasons", "comparison_rt_both_seasons_with_hosps_s2eos.pdf")
dir.create(dirname(out1), recursive = TRUE, showWarnings = FALSE)
ggsave(out1, fig1, width = 14, height = 7.4, dpi = 300); cat("wrote", out1, "\n")

# ---- Figure 2: MechRt pointwise vs simultaneous band ------------------------
plot_bands <- function(res, w) {
  o <- res$rt_ours[res$rt_ours$date >= w[1] & res$rt_ours$date <= w[2],
                   c("date","Rt_mean","Rt_lo","Rt_hi","sim_lo","sim_hi")]
  ggplot() +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_ribbon(data = o, aes(date, ymin = sim_lo, ymax = sim_hi, fill = "Simultaneous"), alpha = 0.18) +
    geom_ribbon(data = o, aes(date, ymin = Rt_lo, ymax = Rt_hi, fill = "Pointwise"), alpha = 0.40) +
    geom_line(data = o, aes(date, Rt_mean), colour = col_ours, linewidth = 0.9) +
    scale_fill_manual(values = c("Pointwise" = col_ours, "Simultaneous" = col_ours),
                      name = "MechRt 90% band", breaks = c("Pointwise","Simultaneous")) +
    rt_scale_y() + labs(x = NULL, y = expression(R[t]), title = sprintf("%s flu season", season_labels[[res$season]])) +
    theme_minimal(base_size = 16) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b", limits = w, expand = expansion(mult = c(0.005, 0.005))) +
    theme(plot.title = element_text(hjust = 0.5, size = 20), axis.text = element_text(size = 17),
          axis.title.y = element_text(size = 20), legend.position = "bottom") +
    guides(fill = guide_legend(override.aes = list(alpha = c(0.40, 0.18))))
}
fig2 <- (plot_bands(all_res$s1, panel_windows$s1) + plot_bands(all_res$s2, panel_windows$s2)) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.title = element_text(size = 18), legend.text = element_text(size = 18))
out2 <- file.path(flu_figures_dir, "both_seasons", "comparison_rt_both_seasons_mechrt_simband_vs_pointwise.pdf")
ggsave(out2, fig2, width = 14, height = 5.4, dpi = 300); cat("wrote", out2, "\n")
