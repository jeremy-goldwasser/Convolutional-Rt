# code/analyses/sim_flu/plot_mechrt_lambda_cv.R
#   -> figures/sim/flu/wiggly/mechrt_lambda_cv_combo.pdf
# 1x2 combo: MechRt fits across the smoothing-penalty grid + 5-fold Poisson-deviance CV curve.
# Reads cached diagnostic objects mechrt_lambda_grid.rds + mechrt_cv_curve.rds.
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(patchwork); library(cowplot); library(ggtext) })
sd_ <- dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))
if (!length(sd_) || is.na(sd_)) sd_ <- "code/analyses/sim_flu"
source(file.path(sd_, "_common.R"))
p <- sim_flu_paths("wiggly")

grid <- readRDS(file.path(p$retro_dir, "mechrt_lambda_grid.rds"))
cvc  <- readRDS(file.path(p$retro_dir, "mechrt_cv_curve.rds"))
lam_min <- cvc$lam_min; log10_lam_min <- log10(lam_min)
cat(sprintf("CV-min lambda = %.4e (log10 = %.3f)\n", lam_min, log10_lam_min))

LAM_LIMS <- c(-2, 8)
lam_scale <- scale_colour_viridis_c(option = "viridis", direction = 1, limits = LAM_LIMS,
  breaks = seq(-2, 8, by = 2), oob = scales::squish, name = expression(log[10] ~ lambda),
  guide = guide_colourbar(direction = "horizontal", title.position = "left", title.vjust = 0.9,
    barwidth = grid::unit(13, "cm"), barheight = grid::unit(0.55, "cm"), ticks.colour = "white", frame.colour = NA))
base_theme <- theme_minimal(base_size = 22) +
  theme(plot.title = element_markdown(hjust = 0.5, size = 30, margin = margin(t = 2, b = 6)),
        axis.title.x = element_text(size = 26, margin = margin(t = 6)), axis.title.y = element_text(size = 26),
        axis.text.x = element_text(size = 24), axis.text.y = element_text(size = 24),
        panel.grid.minor = element_blank(), plot.margin = margin(t = 5, r = 18, b = 5, l = 5))
hi_dark <- "#D55E00"

x_lo <- as.Date("2022-10-01"); x_hi <- as.Date("2023-02-01")
fan <- grid$plugin %>% filter(date >= x_lo, date <= x_hi)
hl  <- grid$cvmin  %>% filter(date >= x_lo, date <= x_hi)
p_fan <- ggplot(fan, aes(date, Rt_mean, group = lam, colour = log10_lam)) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey55", linewidth = 0.5) +
  geom_line(linewidth = 0.7, alpha = 0.75) +
  geom_line(data = hl, aes(date, Rt_mean), colour = "white", linewidth = 2.6, alpha = 0.85, inherit.aes = FALSE) +
  geom_line(data = hl, aes(date, Rt_mean), colour = hi_dark, linewidth = 1.5, inherit.aes = FALSE) +
  lam_scale +
  scale_x_date(breaks = seq(x_lo, x_hi, by = "1 month"), date_labels = "%b", limits = c(x_lo, x_hi), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0.6, 1.4, by = 0.2)) + coord_cartesian(ylim = c(0.6, 1.4)) +
  labs(title = "**(a)** MechRt fits across the penalty grid", x = NULL, y = expression(R[t])) + base_theme

cv <- cvc$cv; mn <- cv %>% filter(lam == lam_min)
y_top <- max(cv$cv_mean + ifelse(is.na(cv$cv_se), 0, cv$cv_se), na.rm = TRUE)
p_cv <- ggplot(cv, aes(log10_lam, cv_mean)) +
  geom_vline(xintercept = log10_lam_min, linetype = "dashed", colour = hi_dark, linewidth = 0.7) +
  geom_errorbar(aes(ymin = cv_mean - cv_se, ymax = cv_mean + cv_se), width = 0.22, colour = "grey75", linewidth = 0.5, na.rm = TRUE) +
  geom_line(colour = "grey55", linewidth = 0.8) + geom_point(aes(colour = log10_lam), size = 3.0) +
  geom_point(data = mn, colour = "white", size = 7.0, inherit.aes = FALSE, aes(log10_lam, cv_mean)) +
  geom_point(data = mn, aes(log10_lam, cv_mean, colour = log10_lam), size = 4.4, inherit.aes = FALSE) +
  geom_point(data = mn, aes(log10_lam, cv_mean), shape = 1, size = 7.4, colour = hi_dark, stroke = 1.3, inherit.aes = FALSE) +
  annotate("text", x = log10_lam_min - 0.4, y = y_top, hjust = 1, vjust = 1, label = "CV min", size = 7, fontface = "bold", colour = hi_dark) +
  lam_scale + scale_x_continuous(breaks = seq(-2, 8, by = 2), expand = c(0.02, 0)) +
  labs(title = "**(b)** Penalty chosen by 5-fold CV", x = expression(log[10] ~ lambda), y = "CV Poisson deviance") + base_theme

cbar <- cowplot::get_plot_component(
  p_fan + theme(legend.position = "bottom", legend.justification = "center",
                legend.title = element_text(size = 26), legend.text = element_text(size = 22),
                legend.box.margin = margin(0, 0, 0, 0)), "guide-box-bottom", return_all = TRUE)
panels <- (p_fan + theme(legend.position = "none") + p_cv + theme(legend.position = "none")) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(title = expression(bold("MechRt smoothing penalty ") * bold(lambda) * bold(": selection on the Flu Simulation")),
                  theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 34, margin = margin(t = 2, b = 6))))
combined <- cowplot::plot_grid(panels, cbar, ncol = 1, rel_heights = c(1, 0.12))
dir.create(p$fig_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(p$fig_dir, "mechrt_lambda_cv_combo.pdf"), combined, width = 20, height = 7.0, dpi = 300)
cat("wrote mechrt_lambda_cv_combo.pdf\n")
