# code/analyses/weekly/run_retro_realtime_combined.R
#
# 1x2 PDF combining a retrospective and a real-time view of weekly girt -- the
# girt port of figures/mechrt/weekly_demo/weekly_retro_realtime_combined.pdf.
#
#   left  : daily vs weekly vs truth on the flu "wiggly" simulation
#           (retrospective; weekly at the 1se rule, daily at lam_min).
#   right : s2 real-time weekly fit at the 2024-01-28 cutoff with the
#           viridis-by-gamma family of Rt curves, retro overlay (black) and
#           FV-min gamma (blue) highlighted.
#
# Reads only cached fits (run run_sim_weekly.R + run_s2_weekly_curves.R first).
# Output: figures/extra_methods/weekly/weekly_retro_realtime_combined.pdf

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(ggplot2)
  library(patchwork); library(scales); library(ggnewscale); library(ggtext)
})
options(warn = 1)

sd <- dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))
if (!length(sd) || is.na(sd)) sd <- "code/analyses/weekly"
source(file.path(sd, "_common.R"))

sim_results <- file.path(wk_results_dir, "sim")
s2_results  <- file.path(wk_results_dir, "s2")
dir.create(wk_figures_dir, recursive = TRUE, showWarnings = FALSE)

# ----- LEFT PANEL: sim retrospective ----------------------------------------
w_sim <- readRDS(file.path(sim_results, "weekly_fit.rds"))
d_sim <- readRDS(file.path(sim_results, "daily_fit.rds"))

left_start <- as.Date("2022-10-01"); left_end <- as.Date("2023-02-01")
clip <- function(df, a, b) df |> filter(!is.na(date), date >= a, date <= b)

rt_daily_w  <- clip(d_sim$rt_daily,      left_start, left_end)
rt_weekly_w <- clip(w_sim$rt_at_lam_1se, left_start, left_end)
truth_w     <- clip(d_sim$rt_daily |> select(date, Rt_truth) |>
                      rename(Rt_mean = Rt_truth), left_start, left_end)

lvls <- c("Daily", "Weekly", "Truth")
lines_df <- bind_rows(
  rt_daily_w  |> transmute(date, Rt_mean, Rt_lo, Rt_hi, fit = "Daily"),
  rt_weekly_w |> transmute(date, Rt_mean, Rt_lo, Rt_hi, fit = "Weekly")
) |> mutate(fit = factor(fit, levels = lvls))
truth_df <- truth_w |> transmute(date, Rt_mean, fit = factor("Truth", levels = lvls))

pal       <- c(Daily = "#1F4E79", Weekly = "#B7472A", Truth = "#2B2B2B")
linetypes <- c(Daily = "solid",   Weekly = "solid",   Truth = "dashed")

draw_key_band_line <- function(data, params, size) {
  is_line_only <- !is.null(data$linetype) && as.character(data$linetype) == "dashed"
  if (is_line_only) return(ggplot2::draw_key_path(data, params, size))
  grid::grobTree(
    grid::rectGrob(gp = grid::gpar(col = NA, fill = scales::alpha(data$colour, 0.18))),
    ggplot2::draw_key_path(data, params, size))
}

all_lines_df <- bind_rows(lines_df, truth_df |> mutate(Rt_lo = NA_real_, Rt_hi = NA_real_))

p_left <- ggplot() +
  geom_ribbon(data = lines_df, aes(x = date, ymin = Rt_lo, ymax = Rt_hi, fill = fit),
              alpha = 0.18, color = NA) +
  geom_line(data = all_lines_df, aes(x = date, y = Rt_mean, color = fit, linetype = fit),
            linewidth = 1.0, key_glyph = draw_key_band_line) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "grey45", linewidth = 0.6) +
  scale_color_manual(values = pal, breaks = lvls, drop = FALSE) +
  scale_fill_manual (values = pal, breaks = c("Daily", "Weekly"), drop = FALSE, guide = "none") +
  scale_linetype_manual(values = linetypes, breaks = lvls, drop = FALSE) +
  labs(title = "**(a)** Retrospective", x = NULL, y = expression(R[t]),
       color = NULL, linetype = NULL) +
  coord_cartesian(xlim = c(left_start, left_end)) +
  theme_minimal(base_size = 18) +
  theme(plot.title = element_markdown(size = 22, hjust = 0.5, margin = margin(b = 4)),
        legend.position = "bottom", legend.text = element_text(size = 21),
        legend.key.width = unit(3.6, "lines"), axis.title.y = element_text(size = 22),
        axis.text = element_text(size = 18), panel.grid.minor = element_blank())

# ----- RIGHT PANEL: s2 real-time at 2024-01-28 ------------------------------
cell_data  <- readRDS(file.path(s2_results, "weekly_combined_rt_curves_data.rds"))
right_cell <- Filter(function(x) format(x$week_end_date) == "2024-01-28", cell_data)[[1]]
retro_s2   <- readRDS(file.path(s2_results, "weekly_fit.rds"))$rt_df

WEEK_END  <- right_cell$week_end_date
win_start <- WEEK_END - 42

rt_plot <- right_cell$rt_per_gamma |> filter(date >= win_start, date <= WEEK_END)
retro_w <- retro_s2 |> filter(date >= win_start, date <= WEEK_END)
idx_min <- which.min(abs(log10(right_cell$gamma_grid) - log10(right_cell$gamma_min)))
rt_min  <- rt_plot |> filter(gamma_index == idx_min)

ymin <- min(c(rt_plot$Rt_mean, retro_w$Rt_mean), na.rm = TRUE) - 0.03
ymax <- max(c(rt_plot$Rt_mean, retro_w$Rt_mean), na.rm = TRUE) + 0.03

highlight_cols <- c("Retrospective" = "#2B2B2B", "FV Min-Rule" = "#2D6CA2")

p_right <- ggplot() +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_line(data = rt_plot, aes(x = date, y = Rt_mean, group = gamma_index, colour = gamma),
            linewidth = 0.5, alpha = 0.78) +
  scale_colour_viridis_c(
    trans = "log10", option = "viridis", name = expression(gamma),
    breaks = 10^seq(-2, 5, by = 2),
    labels = trans_format("log10", math_format(10^.x)),
    guide = guide_colourbar(direction = "vertical", position = "right",
                            title.position = "top", title.hjust = 0.5,
                            barwidth = unit(1.0, "lines"), barheight = unit(13, "lines"),
                            order = 2)) +
  new_scale_colour() +
  geom_line(data = retro_w, aes(x = date, y = Rt_mean, colour = "Retrospective"), linewidth = 1.4) +
  geom_line(data = rt_min,  aes(x = date, y = Rt_mean, colour = "FV Min-Rule"),   linewidth = 1.4) +
  geom_vline(xintercept = WEEK_END, linetype = "dotted", colour = "grey40") +
  scale_colour_manual(values = highlight_cols, name = NULL, breaks = names(highlight_cols),
                      guide = guide_legend(order = 1, position = "bottom",
                                           direction = "horizontal",
                                           keywidth = unit(3.6, "lines"))) +
  scale_x_date(breaks = rev(seq(WEEK_END, win_start, by = -7)),
               date_labels = "%b %d", minor_breaks = NULL) +
  scale_y_continuous(breaks = seq(0, 5, by = 0.1)) +
  coord_cartesian(ylim = c(ymin, ymax)) +
  labs(title = "**(b)** Real-time", x = NULL, y = NULL) +
  theme_minimal(base_size = 18) +
  theme(plot.title = element_markdown(size = 22, hjust = 0.5, margin = margin(b = 4)),
        axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 21), legend.title = element_text(size = 24),
        panel.grid.minor = element_blank())

# ----- Combine --------------------------------------------------------------
combined <- p_left + p_right +
  plot_layout(widths = c(1.15, 1)) +
  plot_annotation(
    title = "girt with weekly hospitalization data",
    theme = theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5,
                                            margin = margin(b = 8))))

out_path <- file.path(wk_figures_dir, "weekly_retro_realtime_combined.pdf")
ggsave(out_path, combined, width = 17, height = 5.5, dpi = 300)
cat(sprintf("Saved: %s\n", out_path))
