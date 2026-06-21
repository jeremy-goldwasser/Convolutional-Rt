# code/analyses/covidestim/delay_sensitivity_LA_GI.R
#
# Los Angeles Rt sensitivity to the generation interval, FROM SCRATCH.
# Builds LA from covidestim_us.rds (fips 06037), recomputes the GI-target variant
# grid (E->I / I->R delays spanning GI means 3.5..6.0 d), fits each with girt.
# Output: figures/real/covidestim/LA_Rt_GI_sensitivity.pdf

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(viridisLite) })
script_dir <- dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))
if (!length(script_dir) || is.na(script_dir)) script_dir <- "code/analyses/covidestim"
source(file.path(script_dir, "_common.R"))

N_LA <- 10040000L
mean_case <- 10.2; sd_case <- 3.6
start_date <- as.Date("2020-03-15"); beta_start_date <- as.Date("2020-04-01")
likelihood_start <- beta_start_date + 21L; retro_end <- as.Date("2021-12-01")
pi_EC <- ce_gamma_pmf(mean_case, sd_case); d_max_EC <- length(pi_EC)

# --- GI-target variant grid (recomputed from scratch) -------------------------
# GI mean = mean of the sum-normalized infectious kernel zeta^EI.
compute_gi_mean <- function(mL, sL, mI, sI) {
  Kn <- gi_from_compartmental(ce_gamma_pmf(mL, sL), ce_gamma_pmf(mI, sI))$g
  sum(seq_along(Kn) * Kn)
}
GI_TARGETS <- c(3.5, 4.0, 4.5, 5.0, 5.5, 6.0)
SD_LAT <- 1.10; SD_INF <- 1.40; RATIO <- 1.27; MEAN_LAT_CAP <- 3.5; MEAN_INF_CAP <- 5.2
tune_for_gi <- function(target) {
  f1 <- function(mL) compute_gi_mean(mL, SD_LAT, mL * RATIO, SD_INF) - target
  mL <- tryCatch(uniroot(f1, c(0.3, 12), extendInt = "yes")$root, error = function(e) NA_real_)
  if (!is.na(mL) && mL <= MEAN_LAT_CAP && mL * RATIO <= MEAN_INF_CAP)
    return(list(mean_lat = mL, sd_lat = SD_LAT, mean_inf = mL * RATIO, sd_inf = SD_INF))
  mL <- MEAN_LAT_CAP
  f2 <- function(mR) compute_gi_mean(mL, SD_LAT, mR, SD_INF) - target
  mR <- tryCatch(uniroot(f2, c(0.3, 30), extendInt = "yes")$root, error = function(e) NA_real_)
  if (!is.na(mR) && mR <= MEAN_INF_CAP)
    return(list(mean_lat = mL, sd_lat = SD_LAT, mean_inf = mR, sd_inf = SD_INF))
  mR <- MEAN_INF_CAP
  f3 <- function(sR) compute_gi_mean(mL, SD_LAT, mR, sR) - target
  sR <- uniroot(f3, c(0.5, 5), extendInt = "yes")$root
  list(mean_lat = mL, sd_lat = SD_LAT, mean_inf = mR, sd_inf = sR)
}
variants <- lapply(GI_TARGETS, function(g) c(list(target_gi = g, is_current = abs(g - 4.5) < 1e-6), tune_for_gi(g)))

# --- LA data + (shared) ascertainment -----------------------------------------
us <- readRDS(file.path(ce_data_dir, "covidestim_us.rds"))
la <- us |> filter(fips == "06037", date >= start_date, date <= retro_end) |> arrange(date) |>
  mutate(cases = pmax(0, round(cases.fitted)), day = row_number() - 1L)
ascertain <- ce_ascertainment(la$cases.fitted, la$infections, pi_EC)

rt_long <- bind_rows(lapply(variants, function(v) {
  pi_lat <- ce_gamma_pmf(v$mean_lat, v$sd_lat); pi_IR <- ce_gamma_pmf(v$mean_inf, v$sd_inf)
  g_v <- gi_from_compartmental(pi_lat, pi_IR,
                               mean_infectious_nominal = v$mean_inf,
                               lag_one = TRUE)$g
  gf <- fit_covidestim_girt(cases = la$cases, ascertain = ascertain, dates = la$date,
                            pi_EC = pi_EC, g = g_v,
                            beta_start_date = beta_start_date, likelihood_start = likelihood_start,
                            use_dow = FALSE)
  rt_g <- gi_extract_rt(gf$fit, gf$design, lam = gf$lambda,
                        level = 0.95, overdispersion = TRUE) |>
    mutate(date = la$date[day + 1L]) |> filter(date >= beta_start_date)
  cat(sprintf("GI=%.1f (mL=%.2f mI=%.2f): lambda=%.4g\n", v$target_gi, v$mean_lat, v$mean_inf, gf$lambda))
  tibble(target_gi = v$target_gi, is_current = v$is_current, date = rt_g$date, Rt_mean = rt_g$Rt_mean)
}))

x_start <- min(la$date) + d_max_EC; x_end <- max(la$date) - d_max_EC
rt_long <- rt_long |> filter(date >= x_start, date <= x_end) |> mutate(grp = factor(target_gi))
targets <- sort(unique(rt_long$target_gi))

p <- ggplot(rt_long, aes(x = date, y = Rt_mean, group = grp, colour = target_gi, linewidth = is_current)) +
  geom_line(alpha = 0.95) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  scale_colour_viridis_c(option = "D", end = 0.92, breaks = targets,
                         labels = sprintf("%.1f", targets), name = "Generation interval (mean days)") +
  scale_linewidth_manual(values = c(`TRUE` = 1.8, `FALSE` = 0.7), guide = "none") +
  scale_x_date(breaks = seq(x_start, x_end, by = "3 months"), date_labels = "%b %Y", limits = c(x_start, x_end)) +
  scale_y_continuous(breaks = seq(0, 5, by = 0.2), minor_breaks = NULL) +
  labs(title = expression(paste("Los Angeles ", R[t], ": sensitivity to generation interval")),
       x = NULL, y = expression(R[t])) +
  guides(colour = guide_colourbar(direction = "horizontal", title.position = "left", title.vjust = 0.9,
    title.hjust = 0, barwidth = grid::unit(22, "lines"), barheight = grid::unit(0.9, "lines"))) +
  theme_minimal(base_size = 16) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 18), axis.text = element_text(size = 14),
        axis.title = element_text(size = 16), legend.position = "bottom", legend.direction = "horizontal",
        legend.text = element_text(size = 14), legend.title = element_text(size = 16, margin = margin(r = 12)))
ggsave(file.path(ce_fig_dir, "LA_Rt_GI_sensitivity.pdf"), p, width = 11, height = 4.8, dpi = 300)
cat(sprintf("wrote %s\n", file.path(ce_fig_dir, "LA_Rt_GI_sensitivity.pdf")))
