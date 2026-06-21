# code/analyses/flu/plot_spaghetti_constant.R
#   -> figures/real/flu/both_seasons/s1_nowcast_s2_spaghetti_constant.pdf
# MechRt real-time constant-tail vs tapered-linear, girt (log). All MechRt from girt tree.

suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(patchwork); library(ggtext); library(viridisLite) })
sd <- dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))
if (!length(sd) || is.na(sd)) sd <- "code/analyses/flu"
source(file.path(sd, "_common.R"))

windows <- list(s1 = c(as.Date("2022-09-24"), as.Date("2023-03-01")),
                s2 = c(as.Date("2023-09-24"), as.Date("2024-04-12")))
col_retro <- "#2D6CA2"; col_tapered <- "#2D6CA2"; col_constant <- "#D55E00"
shared_x_scale <- function(w) scale_x_date(date_breaks = "1 month", date_labels = "%b\n%Y",
                                           limits = w, expand = expansion(mult = c(0.005, 0.005)))
gtree <- function(s) file.path(flu_results_dir, "girt", s)

retro_ss <- function(s) { r <- readRDS(list.files(file.path(gtree(s), "retrospective", "results"), full.names = TRUE)[1])$rt_df
  data.frame(date = as.Date(r$date), Rt_mean = r$Rt_mean) }
read_weekly <- function(s, pat) do.call(rbind, lapply(
  list.files(file.path(gtree(s), "real_time", "results"), pattern = pat, full.names = TRUE), function(f) {
    r <- readRDS(f); rt <- r$rt_df; if (is.null(rt) || !nrow(rt)) return(NULL)
    we <- as.Date(r$meta$week_end_date); rt <- rt[as.Date(rt$date) <= we, ]
    data.frame(date = as.Date(rt$date), Rt_mean = rt$Rt_mean, week_end = we) }))

build_left <- function() {
  w <- windows$s1; lv <- c("Retrospective","Real-time, linear","Real-time, constant")
  rs <- retro_ss("s1"); rs <- rs[rs$date >= w[1] & rs$date <= w[2], ]
  tap <- read_weekly("s1", "natural_linear_taper_cv_1se")
  con <- read_weekly("s1", "regression_constant_notaper_cv_1se")
  seg <- function(df) df[df$date <= df$week_end & df$date >= df$week_end - 6L & df$date >= w[1] & df$date <= w[2], ]
  nc  <- function(df) { d <- df[df$date == df$week_end, c("date","Rt_mean")]; d[d$date >= w[1] & d$date <= w[2], ] }
  retro_df <- transform(rs, series = factor("Retrospective", levels = lv))
  seg_df <- rbind(transform(seg(tap)[,c("date","Rt_mean","week_end")], series = "Real-time, linear"),
                  transform(seg(con)[,c("date","Rt_mean","week_end")], series = "Real-time, constant"))
  seg_df$series <- factor(seg_df$series, levels = lv)
  nc_df <- rbind(transform(nc(tap), series = "Real-time, linear"), transform(nc(con), series = "Real-time, constant"))
  nc_df$series <- factor(nc_df$series, levels = lv)
  ggplot() + geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_line(data = retro_df, aes(date, Rt_mean, colour = series, linetype = series), linewidth = 1.0) +
    geom_line(data = seg_df, aes(date, Rt_mean, colour = series, linetype = series, group = interaction(series, week_end)), linewidth = 1.0) +
    geom_point(data = nc_df, aes(date, Rt_mean, colour = series), size = 1.7, show.legend = FALSE) +
    scale_colour_manual(values = c("Retrospective"=col_retro,"Real-time, linear"=col_tapered,"Real-time, constant"=col_constant), breaks = lv, name = NULL) +
    scale_linetype_manual(values = c("Retrospective"="solid","Real-time, linear"="dashed","Real-time, constant"="dashed"), breaks = lv, name = NULL) +
    shared_x_scale(w) + scale_y_continuous(breaks = seq(0,5,by=0.1)) +
    labs(x = NULL, y = expression(R[t]), title = "2022/23: constant vs tapered linear tail") +
    guides(colour = guide_legend(nrow = 1, byrow = TRUE), linetype = guide_legend(nrow = 1, byrow = TRUE)) +
    theme_minimal(base_size = 18) +
    theme(plot.title = element_markdown(hjust = 0.5, size = 23), legend.position = "bottom", legend.direction = "horizontal",
          legend.text = element_text(size = 19), legend.key.width = grid::unit(3.0, "lines"),
          axis.text = element_text(size = 17), axis.title.y = element_text(size = 20), panel.grid.minor.y = element_blank())
}
build_right <- function() {
  w <- windows$s2; df <- read_weekly("s2", "regression_constant_notaper_cv_1se")
  df <- df[df$date >= w[1] & df$date <= w[2], ]
  ep <- df |> group_by(week_end) |> slice_max(date, n = 1L) |> ungroup()
  ggplot() + geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_line(data = df, aes(date, Rt_mean, group = week_end, colour = week_end), linewidth = 0.55, alpha = 0.8) +
    geom_point(data = ep, aes(date, Rt_mean, colour = week_end), size = 1.4, alpha = 0.95) +
    scale_colour_viridis_c(option = "viridis", trans = "date", name = "Vintage") +
    shared_x_scale(w) + scale_y_continuous(breaks = seq(0,5,by=0.1)) +
    labs(x = NULL, y = expression(R[t]), title = "2023/24: real-time vintages") +
    guides(colour = guide_colourbar(direction = "horizontal", title.position = "left", title.vjust = 0.9, title.hjust = 0,
                                    barwidth = grid::unit(22, "lines"), barheight = grid::unit(0.9, "lines"))) +
    theme_minimal(base_size = 18) +
    theme(plot.title = element_markdown(hjust = 0.5, size = 23), legend.position = "bottom", legend.direction = "horizontal",
          legend.text = element_text(size = 18), legend.title = element_text(size = 20, margin = margin(r = 14)),
          axis.text = element_text(size = 17), axis.title.y = element_text(size = 20))
}
stitched <- (build_left() | build_right()) +
  plot_annotation(title = "MechRt real-time Rt: constant-tail extrapolation",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold")))
out <- file.path(flu_figures_dir, "both_seasons", "s1_nowcast_s2_spaghetti_constant.pdf")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
ggsave(out, stitched, width = 18, height = 7.4, dpi = 300); cat("wrote", out, "\n")
