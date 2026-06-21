# code/analyses/covidestim/ca_retro.R
#
# California retrospective Rt: girt vs CovidEstim, generated FROM SCRATCH.
# Inputs (raw covidestim outputs): covidestim_ca_daily.rds (state aggregate) and
# covidestim_ca_state_rt_band.rds (covidestim's own CA Rt + 80% band).
# Output: figures/real/covidestim/CA_Rt_vs_covidestim.pdf

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })
script_dir <- dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))
if (!length(script_dir) || is.na(script_dir)) script_dir <- "code/analyses/covidestim"
source(file.path(script_dir, "_common.R"))

N_CA <- 39540000L
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

# --- raw CA state data ---------------------------------------------------------
ca <- readRDS(file.path(ce_data_dir, "covidestim_ca_daily.rds")) |>
  filter(date >= start_date, date <= retro_end) |> arrange(date) |>
  mutate(cases = pmax(0, round(cases.fitted)), day = row_number() - 1L)
ca$ascertain <- ce_ascertainment(ca$cases.fitted, ca$infections, pi_EC)

cat("Fitting CA with girt (from scratch)...\n")
gf <- fit_covidestim_girt(cases = ca$cases, ascertain = ca$ascertain, dates = ca$date,
                          pi_EC = pi_EC, g = g_cov,
                          beta_start_date = beta_start_date, likelihood_start = likelihood_start,
                          use_dow = TRUE)

# 80% pointwise band (overdispersion ON), as in the comparison figure
rt_pw80 <- gi_extract_rt(gf$fit, gf$design, lam = gf$lambda,
                         level = 0.80, overdispersion = TRUE) |>
  mutate(date = ca$date[day + 1L]) |> select(date, Rt_mean, Rt_lo, Rt_hi) |>
  filter(date >= beta_start_date)

ce_band  <- readRDS(file.path(ce_data_dir, "covidestim_ca_state_rt_band.rds"))
x_start  <- min(ca$date) + d_max_EC; x_end <- max(ca$date) - d_max_EC
rt_overlay <- rt_pw80 |> rename(Rt_ours = Rt_mean) |> left_join(ce_band, by = "date") |>
  filter(date >= x_start, date <= x_end)
yv <- c(rt_overlay$Rt_ours, rt_overlay$Rt_cov, rt_overlay$Rt_lo, rt_overlay$Rt_hi,
        rt_overlay$Rt_cov_lo, rt_overlay$Rt_cov_hi)
y_lo <- floor(min(yv, na.rm = TRUE) * 10) / 10; y_hi <- ceiling(max(yv, na.rm = TRUE) * 10) / 10

p <- ggplot(rt_overlay, aes(x = date)) +
  geom_ribbon(aes(ymin = Rt_cov_lo, ymax = Rt_cov_hi, fill = "CovidEstim"), alpha = 0.18, na.rm = TRUE) +
  geom_ribbon(aes(ymin = Rt_lo, ymax = Rt_hi, fill = "MechRt"), alpha = 0.38) +
  geom_line(aes(y = Rt_ours, color = "MechRt"), linewidth = 0.8) +
  geom_line(aes(y = Rt_cov,  color = "CovidEstim"), linewidth = 1.0, na.rm = TRUE) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("MechRt" = "#2D6CA2", "CovidEstim" = "#1b9e77")) +
  scale_fill_manual (values = c("MechRt" = "#2D6CA2", "CovidEstim" = "#1b9e77")) +
  scale_x_date(breaks = seq(as.Date("2020-04-01"), as.Date("2021-10-01"), by = "3 months"),
               date_labels = "%b %Y", limits = c(x_start, x_end)) +
  scale_y_continuous(breaks = seq(y_lo, y_hi, by = 0.1), minor_breaks = NULL) +
  coord_cartesian(ylim = c(y_lo, y_hi)) +
  labs(title = expression(paste("California ", R[t], ": MechRt vs CovidEstim")),
       x = NULL, y = expression(R[t]), color = NULL, fill = NULL) +
  theme_minimal(base_size = 16) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 18), legend.position = "bottom",
        legend.text = element_text(size = 18), legend.key.width = unit(1.6, "cm"))
ggsave(file.path(ce_fig_dir, "CA_Rt_vs_covidestim.pdf"), p, width = 11, height = 4.8, dpi = 300)
cat(sprintf("CA: lambda=%.4g | Pearson r (girt vs covidestim) = %.3f\n",
            gf$lambda, cor(rt_overlay$Rt_ours, rt_overlay$Rt_cov, use = "complete.obs")))
cat(sprintf("wrote %s\n", file.path(ce_fig_dir, "CA_Rt_vs_covidestim.pdf")))
