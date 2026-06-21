# code/analyses/covidestim/counties_retro.R
#
# Four-county retrospective Rt: girt vs CovidEstim, generated FROM SCRATCH from
# covidestim_us.rds (county-level cases.fitted / infections / Rt).
# Output: figures/real/covidestim/Rt_vs_covidestim_panel.pdf

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(patchwork) })
script_dir <- dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))
if (!length(script_dir) || is.na(script_dir)) script_dir <- "code/analyses/covidestim"
source(file.path(script_dir, "_common.R"))

counties <- list(
  list(fips = "17031", label = "Cook_IL",        title = "Cook County, IL (Chicago)",   pop = 5150233L),
  list(fips = "48201", label = "Harris_TX",      title = "Harris County, TX (Houston)",  pop = 4731145L),
  list(fips = "06037", label = "Los_Angeles_CA", title = "Los Angeles County, CA",       pop = 10014009L),
  list(fips = "12086", label = "Miami_Dade_FL",  title = "Miami-Dade County, FL",        pop = 2716940L))

mean_latent <- 2.88; sd_latent <- 1.10
mean_infectious <- 3.65; sd_infectious <- 1.40
mean_case <- 10.2; sd_case <- 3.6
start_date <- as.Date("2020-03-15"); beta_start_date <- as.Date("2020-04-01")
likelihood_start <- beta_start_date + 21L; retro_end <- as.Date("2021-12-01")
pi_lat <- ce_gamma_pmf(mean_latent, sd_latent)
pi_IR  <- ce_gamma_pmf(mean_infectious, sd_infectious)
pi_EC  <- ce_gamma_pmf(mean_case, sd_case); d_max_EC <- length(pi_EC)
# girt GI: covidestim MechRt convention (lag-1 + nominal mip).
g_cov  <- gi_from_compartmental(pi_lat, pi_IR,
                                mean_infectious_nominal = mean_infectious,
                                lag_one = TRUE)$g

us <- readRDS(file.path(ce_data_dir, "covidestim_us.rds"))

fit_one <- function(spec) {
  cty <- us |> filter(fips == spec$fips, date >= start_date, date <= retro_end) |> arrange(date) |>
    mutate(cases = pmax(0, round(cases.fitted)), day = row_number() - 1L)
  cty$ascertain <- ce_ascertainment(cty$cases.fitted, cty$infections, pi_EC)
  gf <- fit_covidestim_girt(cases = cty$cases, ascertain = cty$ascertain, dates = cty$date,
                            pi_EC = pi_EC, g = g_cov,
                            beta_start_date = beta_start_date, likelihood_start = likelihood_start,
                            use_dow = FALSE)
  rt_g <- gi_extract_rt(gf$fit, gf$design, lam = gf$lambda,
                        level = 0.95, overdispersion = FALSE) |>
    mutate(date = cty$date[day + 1L]) |> select(date, Rt_mean) |> filter(date >= beta_start_date)
  rj <- rt_g |> left_join(cty |> select(date, Rt_cov = Rt), by = "date")
  cat(sprintf("%-16s lambda=%.4g | r(girt,covidestim)=%.3f\n", spec$label, gf$lambda,
              cor(rj$Rt_mean, rj$Rt_cov, use = "complete.obs")))
  list(spec = spec, rt = rt_g, cty = cty)
}
results <- lapply(counties, fit_one)

global_y <- do.call(c, lapply(results, function(it) {
  xs <- min(it$cty$date) + d_max_EC; xe <- max(it$cty$date) - d_max_EC
  o <- it$rt |> rename(Rt_ours = Rt_mean) |> left_join(it$cty |> select(date, Rt_cov = Rt), by = "date") |>
    filter(date >= xs, date <= xe); c(o$Rt_ours, o$Rt_cov) }))
y_lo <- floor(min(global_y, na.rm = TRUE) * 10) / 10; y_hi <- ceiling(max(global_y, na.rm = TRUE) * 10) / 10

panel_one <- function(it, with_legend = FALSE) {
  xs <- min(it$cty$date) + d_max_EC; xe <- max(it$cty$date) - d_max_EC
  ov <- it$rt |> rename(Rt_ours = Rt_mean) |> left_join(it$cty |> select(date, Rt_cov = Rt), by = "date") |>
    filter(date >= xs, date <= xe)
  ggplot(ov, aes(x = date)) +
    geom_line(aes(y = Rt_ours, color = "MechRt"), linewidth = 1.0) +
    geom_line(aes(y = Rt_cov,  color = "CovidEstim"), linewidth = 0.9, na.rm = TRUE) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
    scale_color_manual(values = c("MechRt" = "#2D6CA2", "CovidEstim" = "#1b9e77")) +
    scale_x_date(breaks = seq(xs, xe, by = "4 months"), date_labels = "%b %Y", limits = c(xs, xe)) +
    scale_y_continuous(breaks = seq(y_lo, y_hi, by = 0.2), minor_breaks = NULL) +
    coord_cartesian(ylim = c(y_lo, y_hi)) +
    labs(title = it$spec$title, x = NULL, y = expression(R[t]), color = NULL) +
    theme_minimal(base_size = 17) +
    theme(plot.title = element_text(hjust = 0.5, size = 19), axis.title.y = element_text(size = 18),
          axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),
          legend.position = if (with_legend) "bottom" else "none",
          legend.text = element_text(size = 20), legend.key.width = unit(2.0, "cm"),
          panel.grid.minor = element_blank())
}
panels <- lapply(seq_along(results), function(i) panel_one(results[[i]], with_legend = (i == 1L)))
combined <- (panels[[1]] | panels[[2]]) / (panels[[3]] | panels[[4]]) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
combined <- combined + plot_annotation(
  title = expression(bold(paste("Estimated ", bolditalic(R)[bold(t)], ": MechRt vs CovidEstim"))),
  theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 22)))
ggsave(file.path(ce_fig_dir, "Rt_vs_covidestim_panel.pdf"), combined, width = 14, height = 6.5, dpi = 300)
cat(sprintf("wrote %s\n", file.path(ce_fig_dir, "Rt_vs_covidestim_panel.pdf")))
