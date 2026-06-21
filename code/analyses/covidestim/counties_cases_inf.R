# code/analyses/covidestim/counties_cases_inf.R
#
# 2x2 panel: reported cases vs CovidEstim latent infections for four counties.
# Pure data plot (no Rt fitting), FROM SCRATCH from covidestim_us.rds.
# Output: figures/real/covidestim/cases_infections_panel.pdf

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(patchwork); library(scales) })
script_dir <- dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))
if (!length(script_dir) || is.na(script_dir)) script_dir <- "code/analyses/covidestim"
source(file.path(script_dir, "_common.R"))

counties <- list(
  list(fips = "06037", title = "Los Angeles County, CA"),
  list(fips = "17031", title = "Cook County, IL (Chicago)"),
  list(fips = "48201", title = "Harris County, TX (Houston)"),
  list(fips = "12086", title = "Miami-Dade County, FL"))
start_date <- as.Date("2020-03-15"); retro_end <- as.Date("2021-12-01")
col_cases <- "#0072B2"; col_inf <- "#E69F00"

us <- readRDS(file.path(ce_data_dir, "covidestim_us.rds"))

make_panel <- function(spec, with_legend = FALSE) {
  cty <- us |> filter(fips == spec$fips, date >= start_date, date <= retro_end) |> arrange(date)
  df <- bind_rows(
    tibble(date = cty$date, value = pmax(0, round(cty$cases.fitted)), series = "Reported cases"),
    tibble(date = cty$date, value = cty$infections, series = "Latent infections (CovidEstim)")) |>
    mutate(series = factor(series, levels = c("Reported cases", "Latent infections (CovidEstim)")))
  x_min <- min(df$date); x_max <- max(df$date)
  ggplot(df, aes(x = date, y = value, color = series)) +
    geom_line(linewidth = 0.7) +
    scale_color_manual(values = c("Reported cases" = col_cases, "Latent infections (CovidEstim)" = col_inf)) +
    scale_x_date(breaks = seq(x_min, x_max, by = "4 months"), date_labels = "%b %Y", limits = c(x_min, x_max)) +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(title = spec$title, x = NULL, y = "Daily count", color = NULL) +
    theme_minimal(base_size = 16) +
    theme(plot.title = element_text(hjust = 0.5, size = 18), axis.text = element_text(size = 14),
          axis.title = element_text(size = 16), legend.position = if (with_legend) "bottom" else "none",
          legend.text = element_text(size = 18), legend.key.width = unit(1.6, "cm"))
}
panels <- lapply(seq_along(counties), function(i) make_panel(counties[[i]], with_legend = (i == 1L)))
combined <- (panels[[1]] | panels[[2]]) / (panels[[3]] | panels[[4]]) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.text = element_text(size = 18))
combined <- combined + plot_annotation(title = "Reported cases vs. CovidEstim latent infections",
  theme = theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")))
ggsave(file.path(ce_fig_dir, "cases_infections_panel.pdf"), combined, width = 14, height = 7.5, dpi = 300)
cat(sprintf("wrote %s\n", file.path(ce_fig_dir, "cases_infections_panel.pdf")))
