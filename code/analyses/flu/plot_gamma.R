# code/analyses/flu/plot_gamma.R
#   -> figures/real/flu/s2/combined_fv_errors.pdf
#   -> figures/real/flu/s2/combined_rt_curves_2x2.pdf
# Gamma-sweep figures from the girt gamma data (results/real/flu/girt/s2/.../gamma/).

suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(patchwork)
  library(scales); library(viridisLite); library(ggnewscale) })
sd <- dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))
if (!length(sd) || is.na(sd)) sd <- "code/analyses/flu"
source(file.path(sd, "_common.R"))
gamma_dir <- file.path(flu_results_dir, "girt", "s2", "real_time", "gamma")
out_dir   <- file.path(flu_figures_dir, "s2"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

lambda_expr <- function(lam) { if (!is.finite(lam) || lam <= 0) return(bquote(lambda == .(lam)))
  e <- floor(log10(lam)); m <- round(lam / 10^e); if (m == 10) { m <- 1; e <- e + 1 }; bquote(lambda == .(m) %*% 10^.(e)) }
panel_title_expr <- function(we, lam) bquote(.(format(we, "%b %d, %Y")) ~ "  " ~ .(lambda_expr(lam)))

build_rt_panel <- function(date, rds, y_floor = NULL, base_size = 18, title_bold = TRUE, barw = 19) {
  d <- readRDS(file.path(gamma_dir, date, rds)); WE <- d$week_end_date; ws <- WE - 42
  rt <- d$rt_per_gamma %>% filter(date >= ws, date <= WE)
  retro <- d$retro_rt_curve %>% filter(date >= ws, date <= WE)
  idx <- which.min(abs(log10(d$gamma_grid) - log10(d$best_gamma_min)))
  rt_min <- rt %>% filter(gamma_index == idx)
  ymin <- min(c(rt$Rt_mean, retro$Rt_mean), na.rm = TRUE) - 0.03; ymax <- max(c(rt$Rt_mean, retro$Rt_mean), na.rm = TRUE) + 0.03
  if (!is.null(y_floor)) ymin <- y_floor
  ggplot() + geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_line(data = rt, aes(date, Rt_mean, group = gamma_index, colour = gamma), linewidth = 0.45, alpha = 0.75) +
    scale_colour_viridis_c(trans = "log10", option = "viridis", name = expression(gamma),
      breaks = 10^seq(0, 5, by = 1), labels = trans_format("log10", math_format(10^.x)),
      guide = guide_colourbar(direction = "horizontal", title.position = "left", title.vjust = 0.85,
                              barwidth = unit(barw, "lines"), barheight = unit(0.8, "lines"), order = 2)) +
    new_scale_colour() +
    geom_line(data = retro, aes(date, Rt_mean, colour = "Retrospective"), linewidth = 1.3) +
    geom_line(data = rt_min, aes(date, Rt_mean, colour = "FV Min-Rule"), linewidth = 1.3) +
    geom_vline(xintercept = WE, linetype = "dotted", colour = "grey40") +
    scale_colour_manual(values = c("Retrospective" = "black", "FV Min-Rule" = "#2D6CA2"), name = NULL,
      breaks = c("Retrospective", "FV Min-Rule"), guide = guide_legend(order = 1, keywidth = unit(3, "lines"))) +
    scale_x_date(breaks = rev(seq(WE, ws, by = -7)), date_labels = "%b %d", minor_breaks = NULL) +
    scale_y_continuous(breaks = seq(0, 5, by = 0.1)) + coord_cartesian(ylim = c(ymin, ymax)) +
    labs(x = NULL, y = expression(R[t]), title = panel_title_expr(WE, d$lam_used)) +
    theme_minimal(base_size = base_size) +
    theme(plot.title = element_text(face = if (title_bold) "bold" else "plain", hjust = 0.5, size = 20),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 16), axis.text.y = element_text(size = 16),
          axis.title.y = element_text(size = 20), panel.grid.minor = element_blank())
}
build_fv_panel <- function(date, rds, show_y = TRUE) {
  d <- readRDS(file.path(gamma_dir, date, rds))
  fv <- data.frame(gamma = d$gamma_grid, fv_mean = d$fv_mean, fv_se = d$fv_se)
  idx <- which.min(abs(log10(d$gamma_grid) - log10(d$best_gamma_min)))
  mk <- data.frame(rule = "FV Min-Rule", gamma = d$gamma_grid[idx], fv_mean = d$fv_mean[idx])
  ggplot(fv, aes(gamma, fv_mean)) +
    geom_errorbar(aes(ymin = fv_mean - fv_se, ymax = fv_mean + fv_se), width = 0.12, linewidth = 0.55, colour = "grey40") +
    geom_line(colour = "grey30", linewidth = 0.7) + geom_point(colour = "grey30", size = 2.1) +
    geom_point(data = mk, aes(gamma, fv_mean, colour = rule), size = 5.6, stroke = 1.8, shape = 21, fill = "white") +
    scale_x_log10(breaks = 10^seq(0, 5, by = 1), labels = trans_format("log10", math_format(10^.x))) +
    scale_colour_manual(values = c("FV Min-Rule" = "#2D6CA2"), name = NULL) +
    labs(x = expression(gamma), y = if (show_y) "FV error" else NULL, title = format(d$week_end_date, "%b %d, %Y")) +
    theme_minimal(base_size = 15) +
    theme(plot.title = element_text(hjust = 0.5, size = 17), axis.title.x = element_text(size = 16),
          axis.text = element_text(size = 13),
          axis.text.y = if (show_y) element_text(size = 13) else element_blank(),
          axis.title.y = if (show_y) element_text(size = 16) else element_blank(), panel.grid.minor = element_blank())
}

# ---- combined_fv_errors.pdf (1x2): 2023-12-31/s1retrolam, 2024-01-28/1se -----
fv <- (build_fv_panel("2023-12-31", "data_s1retrolam.rds") | build_fv_panel("2024-01-28", "data_1se.rds")) +
  plot_annotation(title = "Real-time influenza FV error across taper strength",
    subtitle = expression("Blue circle: FV min-rule " * gamma),
    theme = theme(plot.title = element_text(hjust = 0.5, size = 20, margin = margin(b = 2)),
                  plot.subtitle = element_text(hjust = 0.5, size = 13, colour = "grey30", margin = margin(b = 4)))) &
  theme(legend.position = "none")
ggsave(file.path(out_dir, "combined_fv_errors.pdf"), fv, width = 12, height = 4.4, dpi = 300)
cat("wrote combined_fv_errors.pdf\n")

# ---- combined_rt_curves_2x2.pdf: top=s1retrolam, bottom=retrolam x 2 dates ---
rt2x2 <- (build_rt_panel("2023-12-31","data_s1retrolam.rds", y_floor=0.8) | build_rt_panel("2024-01-28","data_s1retrolam.rds")) /
         (build_rt_panel("2023-12-31","data_retrolam.rds",   y_floor=0.8) | build_rt_panel("2024-01-28","data_retrolam.rds")) +
  plot_annotation(title = "Real-time influenza Rt across taper strength",
    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 22))) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.title = element_text(size = 21), legend.text = element_text(size = 20))
ggsave(file.path(out_dir, "combined_rt_curves_2x2.pdf"), rt2x2, width = 13, height = 8.5, dpi = 300)
cat("wrote combined_rt_curves_2x2.pdf\n")
