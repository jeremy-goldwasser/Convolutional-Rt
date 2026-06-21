# code/analyses/sim_flu/plot_combo_4methods.R
#   -> figures/sim/flu/wiggly/retro_realtime_4methods_wiggly_ci.pdf
#   -> figures/sim/flu/noisy_smooth/retro_realtime_4methods_smooth_ci.pdf
# Unified 4-methods retro + realtime combo with MechRt 95% CI ribbon.  The cached
# weekly fits already carry post-conformal Rt_lo/Rt_hi; retro is Wald (overdisp).
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(ggplot2)
  library(patchwork); library(cowplot); library(ggtext)
})
sd_ <- dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))
if (!length(sd_) || is.na(sd_)) sd_ <- "code/analyses/sim_flu"
source(file.path(sd_, "_common.R"))

methods_keep <- c("MechRt", "EpiNow2", "estimateR", "EpiEstim")
legend_breaks <- c(methods_keep, "Ground Truth")
legend_lwidth <- c(MechRt = 1.3, EpiNow2 = 1.3, estimateR = 1.3, EpiEstim = 1.3, "Ground Truth" = 1.3)
legend_scale <- scale_colour_manual(values = method_palette, name = NULL, breaks = legend_breaks,
  guide = guide_legend(nrow = 1, override.aes = list(
    linetype = c(rep("solid", length(methods_keep)), "dashed"),
    linewidth = unname(legend_lwidth))))
x_lo <- as.Date("2022-10-01"); x_hi <- as.Date("2023-02-01")
x_breaks <- seq(x_lo, x_hi, by = "1 month")
base_theme <- theme_minimal(base_size = 22) +
  theme(plot.title = element_markdown(hjust = 0.5, size = 30, margin = margin(t = 2, b = 4)),
        legend.position = "bottom", legend.text = element_text(size = 28),
        legend.key.width = grid::unit(3.5, "lines"), axis.title.y = element_text(size = 26),
        axis.text.x = element_text(size = 24), axis.text.y = element_text(size = 24),
        panel.grid.minor = element_blank(), plot.margin = margin(t = 5, r = 18, b = 5, l = 5))

z95 <- qnorm(0.975)

render <- function(scenario, out_name, variant_label) {
  pth <- sim_flu_paths(scenario)
  sim_all <- readRDS(pth$sim_rds); sim_all$date <- as.Date(sim_all$date)
  truth_col <- pth$truth_col   # rt_raw or rt_loess
  variant_filter <- if (scenario == "wiggly") "raw" else "loess"   # match retro_sim_rt_estimates$variant
  ep_pattern <- sprintf("__%s_gp__summary\\.rds$", variant_filter)

  # ---- retro panel
  retro_truth <- sim_all |> transmute(date, Rt_truth = .data[[truth_col]]) |>
    filter(date >= x_lo, date <= x_hi)
  rt_all <- readRDS(file.path(pth$retro_dir, "retro_sim_rt_estimates.rds"))
  rt_all$date <- as.Date(rt_all$date)
  label_map <- c("MechRt" = "MechRt", "EpiNow2 GP" = "EpiNow2", "estimateR" = "estimateR", "EpiEstim" = "EpiEstim")
  retro_lines <- rt_all |>
    filter(as.character(method) %in% names(label_map), variant == variant_filter,
           date >= x_lo, date <= x_hi) |>
    mutate(Method = factor(unname(label_map[as.character(method)]), levels = methods_keep)) |>
    select(Method, date, Rt_mean, Rt_lo, Rt_hi)
  # If EpiNow2 retro is missing from the aggregate, splice in from epinow2_retro/.
  if (!"EpiNow2" %in% as.character(retro_lines$Method)) {
    fs <- list.files(pth$ep_retro, ep_pattern, full.names = TRUE)
    if (length(fs)) {
      ep <- readRDS(fs[1])$rt_df
      retro_lines <- bind_rows(retro_lines, tibble(
        Method = factor("EpiNow2", levels = methods_keep),
        date = as.Date(ep$date), Rt_mean = ep$Rt_mean,
        Rt_lo = if (!is.null(ep$Rt_sd)) ep$Rt_mean - z95 * ep$Rt_sd else NA_real_,
        Rt_hi = if (!is.null(ep$Rt_sd)) ep$Rt_mean + z95 * ep$Rt_sd else NA_real_) |>
        filter(date >= x_lo, date <= x_hi))
    }
  }
  retro_mechrt_ci <- retro_lines |> filter(Method == "MechRt")

  p_retro <- ggplot() +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_ribbon(data = retro_mechrt_ci, aes(date, ymin = Rt_lo, ymax = Rt_hi),
                fill = method_palette[["MechRt"]], alpha = 0.18, colour = NA) +
    geom_line(data = retro_lines |> filter(Method != "MechRt"),
              aes(date, Rt_mean, colour = Method, group = Method), linewidth = 1.3, alpha = 0.95) +
    geom_line(data = retro_lines |> filter(Method == "MechRt"),
              aes(date, Rt_mean, colour = Method, group = Method), linewidth = 1.3, alpha = 0.95) +
    geom_line(data = retro_truth, aes(date, Rt_truth, colour = "Ground Truth"),
              linewidth = 1.3, linetype = "dashed") +
    legend_scale +
    scale_x_date(breaks = x_breaks, date_labels = "%b", limits = c(x_lo, x_hi), expand = c(0, 0)) +
    labs(title = "**(a)** Retrospective (End-of-Season)", x = NULL, y = expression(R[t])) +
    coord_cartesian(ylim = c(0.6, 1.35)) + base_theme

  # ---- realtime panel (trailing-7-day segments)
  D_MAX <- 6L
  man <- read.csv(pth$manifest); man$week_end_date <- as.Date(man$week_end_date)
  load_mech <- do.call(rbind, lapply(seq_len(nrow(man)), function(i) {
    f <- file.path(pth$weekly_dir, man$result_filename[i])
    if (!file.exists(f)) return(NULL)
    r <- readRDS(f); rt <- r$rt_df
    data.frame(date = as.Date(rt$date), week_end = man$week_end_date[i],
               Rt_mean = rt$Rt_mean, Rt_lo = rt$Rt_lo, Rt_hi = rt$Rt_hi, stringsAsFactors = FALSE)
  })) |> mutate(d_to_edge = as.integer(week_end - date)) |>
    filter(d_to_edge >= 0L, d_to_edge <= D_MAX)
  fs <- list.files(pth$ep_rt, pattern = ep_pattern, full.names = TRUE)
  load_ep <- do.call(rbind, lapply(fs, function(f) {
    r <- readRDS(f); rt <- r$rt_df
    data.frame(date = as.Date(rt$date), week_end = as.Date(r$meta$week_end_date),
               Rt_mean = rt$Rt_mean, stringsAsFactors = FALSE)
  })) |> mutate(d_to_edge = as.integer(week_end - date)) |>
    filter(d_to_edge >= 0L, d_to_edge <= D_MAX)
  load_other <- function(method) {
    path <- file.path(pth$other_rt, sprintf("%s__%s.rds", method, variant_filter))
    if (!file.exists(path)) return(NULL)
    readRDS(path)$rt |>
      mutate(date = as.Date(date), week_end = as.Date(week_end),
             d_to_edge = as.integer(week_end - date)) |>
      filter(d_to_edge >= 0L, d_to_edge <= D_MAX) |>
      select(date, week_end, d_to_edge, Rt_mean)
  }
  rt_segs <- bind_rows(
    load_mech |> mutate(Method = "MechRt"),
    load_ep   |> mutate(Method = "EpiNow2"),
    load_other("estimateR") |> mutate(Method = "estimateR"),
    load_other("EpiEstim")  |> mutate(Method = "EpiEstim")
  ) |> filter(Method %in% methods_keep)
  seg_window <- as.Date(c("2022-10-01", "2023-01-28"))
  rt_plot <- rt_segs |> filter(date >= seg_window[1], date <= seg_window[2])
  truth_seg <- sim_all |> transmute(date, Rt_truth = .data[[truth_col]]) |>
    filter(date >= seg_window[1], date <= seg_window[2])
  rt_plot$Method <- factor(rt_plot$Method, levels = methods_keep)
  rt_mechrt_ci <- rt_plot |> filter(Method == "MechRt")

  p_rt <- ggplot() +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60") +
    geom_ribbon(data = rt_mechrt_ci, aes(date, ymin = Rt_lo, ymax = Rt_hi, group = week_end),
                fill = method_palette[["MechRt"]], alpha = 0.18, colour = NA) +
    geom_line(data = rt_plot |> filter(Method != "MechRt"),
              aes(date, Rt_mean, group = interaction(Method, week_end), colour = Method),
              linewidth = 1.0, alpha = 0.9) +
    geom_line(data = rt_plot |> filter(Method == "MechRt"),
              aes(date, Rt_mean, group = interaction(Method, week_end), colour = Method),
              linewidth = 1.0, alpha = 0.9) +
    geom_line(data = truth_seg, aes(date, Rt_truth, colour = "Ground Truth"),
              linewidth = 1.3, linetype = "dashed") +
    legend_scale +
    scale_x_date(breaks = x_breaks, date_labels = "%b", limits = c(x_lo, x_hi), expand = c(0, 0)) +
    labs(title = "**(b)** Real-Time (Trailing 7 Days)", x = NULL, y = expression(R[t])) +
    coord_cartesian(ylim = c(0.6, 1.35)) + base_theme

  legend_grob <- cowplot::get_plot_component(
    p_retro + theme(legend.position = "bottom", legend.justification = "center",
                    legend.box.margin = margin(0, 0, 0, 0)),
    "guide-box-bottom", return_all = TRUE)
  panels <- (p_retro + theme(legend.position = "none") + p_rt + theme(legend.position = "none")) +
    plot_layout(widths = c(1, 1)) +
    plot_annotation(title = expression(bold(R[t]) ~ bold("Estimates, Flu Simulation")),
                    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 36,
                                                            margin = margin(t = 2, b = 6))))
  combined <- cowplot::plot_grid(panels, legend_grob, ncol = 1, rel_heights = c(1, 0.08))
  dir.create(pth$fig_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(pth$fig_dir, out_name), combined, width = 20, height = 6.6, dpi = 300)
  cat("wrote", out_name, "\n")
}

render("wiggly",       "retro_realtime_4methods_wiggly_ci.pdf", "Wiggly")
render("noisy_smooth", "retro_realtime_4methods_smooth_ci.pdf", "Smooth")
