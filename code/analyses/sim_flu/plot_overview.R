# code/analyses/sim_flu/plot_overview.R
#   -> figures/sim/flu/wiggly/tuning_sim_overview_wiggly.pdf
#   -> figures/sim/flu/noisy_smooth/tuning_sim_overview_smooth.pdf
# Pure sim viz (true Rt + hosps).  Ported from code/sim/flu{3,4}/overview_plot.R.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(patchwork) })
sd_ <- dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))
if (!length(sd_) || is.na(sd_)) sd_ <- "code/analyses/sim_flu"
source(file.path(sd_, "_common.R"))

col_wiggly <- "#B7472A"; col_smooth <- "#1F4E79"
variant_levels <- c("Wiggly", "Smooth"); variant_cols <- c(Wiggly = col_wiggly, Smooth = col_smooth)
base_theme <- theme_minimal(base_size = 16) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 18), legend.key.width = unit(1.6, "cm"),
        plot.title = element_text(hjust = 0.5, size = 20), plot.subtitle = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14),
        axis.text.y = element_text(size = 14))

make_rt_panel <- function(rt_long) {
  ggplot(rt_long, aes(date, Rt, colour = variant)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_line(linewidth = 0.9) +
    scale_colour_manual(values = variant_cols) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b") +
    scale_y_continuous(breaks = scales::breaks_width(0.1)) +
    labs(title = expression("True"~italic(R)[t]), x = NULL, y = NULL) +
    base_theme + theme(legend.position = "none",
                       plot.title = element_text(hjust = 0.5, size = 17, face = "plain"),
                       axis.text.x = element_text(angle = 0, hjust = 0.5, size = 15))
}
make_hosp_panel <- function(hosp_long) {
  ggplot() +
    geom_point(data = filter(hosp_long, kind == "observed"),
               aes(date, value, colour = variant), size = 0.5, alpha = 0.5, show.legend = FALSE) +
    geom_line(data = filter(hosp_long, kind == "expected"),
              aes(date, value, colour = variant), linewidth = 0.9, show.legend = FALSE) +
    scale_colour_manual(values = variant_cols, guide = "none") +
    scale_x_date(date_breaks = "1 month", date_labels = "%b") +
    labs(title = "Hospitalizations", x = NULL, y = NULL) +
    base_theme + theme(plot.title = element_text(hjust = 0.5, size = 17, face = "plain"),
                       axis.text.x = element_text(angle = 0, hjust = 0.5, size = 15))
}

render <- function(scenario, out_name, variant_label) {
  p <- sim_flu_paths(scenario)
  sim_df <- readRDS(p$sim_rds)
  curve_df <- sim_df |> filter(rt_raw != 1 | rt_loess != 1) |> filter(date >= as.Date("2022-09-01"))
  raw_col  <- if (variant_label == "Wiggly") "rt_raw" else "rt_loess"
  obs_col  <- if (variant_label == "Wiggly") "obs_cases_raw" else "obs_cases_loess"
  exp_col  <- if (variant_label == "Wiggly") "expected_cases_raw" else "expected_cases_loess"
  rt_long  <- tibble(date = as.Date(curve_df$date), variant = factor(variant_label, levels = variant_levels),
                     Rt = curve_df[[raw_col]])
  hosp_long <- bind_rows(
    tibble(date = as.Date(curve_df$date), variant = factor(variant_label, levels = variant_levels),
           kind = "observed", value = curve_df[[obs_col]]),
    tibble(date = as.Date(curve_df$date), variant = factor(variant_label, levels = variant_levels),
           kind = "expected", value = curve_df[[exp_col]]))
  combined <- (make_rt_panel(rt_long) | make_hosp_panel(hosp_long))
  dir.create(p$fig_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(p$fig_dir, out_name), combined, width = 11, height = 3.2, dpi = 300)
  cat("wrote", out_name, "\n")
}

render("wiggly",       "tuning_sim_overview_wiggly.pdf", "Wiggly")
render("noisy_smooth", "tuning_sim_overview_smooth.pdf", "Smooth")
