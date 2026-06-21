# code/examples/demo.R
#
# Self-contained demo for girt.  Simulate an outbreak with a knowingly
# non-trivial R_t trajectory (growth, peak, decline through 1, dip below 1,
# small rebound), recover R_t with girt retrospectively and in real time, and
# plot the result.  No external data; depends only on girt.R.

suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tibble); library(patchwork) })

# Locate + load girt.R relative to this script (so the demo Just Works no
# matter the working dir).
.this <- (function() {
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) normalizePath(f) else NULL
})()
girt_root <- if (!is.null(.this)) dirname(dirname(.this)) else file.path(getwd(), "code")
source(file.path(girt_root, "girt.R"))

set.seed(2)

# ---- Truth: 170 displayed days, two-wave R_t -------------------------------
#   Rt(t) = 1 + 0.50 * gaussian(t, 55, 25)       <- peak ~1.50 around day 55
#         - 0.28 * gaussian(t, 120, 22)          <- dip ~0.72
#         + 0.30 * gaussian(t, 155, 18)          <- rebound to ~1.30
# 120-day pre-window burn-in (renewal kernel + reporting delay both settle).
# Window ends at day 170 (within the rebound) so right-edge counts are still
# healthy -- below ~5 cases/day the convolution likelihood has too little
# signal to resolve R_t.
n_show   <- 170L
burn_in  <- 120L
n_total  <- n_show + burn_in
tt_all   <- seq.int(-burn_in + 1L, n_show)
rt_fn <- function(t) {
  1.0 +
  0.50 * exp(-((t -  55) / 25)^2) -
  0.28 * exp(-((t - 120) / 22)^2) +
  0.30 * exp(-((t - 155) / 18)^2)
}
Rt_true <- rt_fn(tt_all)
# Pre-window Rt is whatever rt_fn gives at t <= 0 (essentially 1, by design).

# Generation interval + reporting delay (generic; not tied to any pathogen).
g     <- gi_discrete_gamma_delay(mean = 3.5, sd = 1.8)$pmf
mean_EY <- 5.7; sd_EY <- 2.3; pi_EY <- gi_discrete_gamma_delay(mean_EY, sd_EY)$pmf
rho   <- 0.035   # severity rate; tuned so peak observed Y > 250 (healthy SNR)

# Simulate exposures via the renewal equation X_t ~ Pois(Rt * sum_k g_k X_{t-k}).
# Start with a small steady seed; the burn-in lets transmission settle into a
# realistic shape before we expose the displayed window.
X <- numeric(n_total); X[seq_len(min(5L, length(g)))] <- 30
for (t in (min(5L, length(g)) + 1L):n_total) {
  k  <- 1:min(length(g), t - 1L)
  lam <- Rt_true[t] * sum(g[k] * X[t - k])
  X[t] <- rpois(1L, lam)
}

# Observations Y_t = Pois(rho * sum_k pi_EY[k] X_{t-k}).
muY <- vapply(seq_len(n_total), function(t) {
  k <- 1:min(length(pi_EY), t - 1L); if (t == 1L) 0 else rho * sum(pi_EY[k] * X[t - k])
}, numeric(1))
Y <- rpois(n_total, pmax(muY, 0))

# Crop to the displayed window (days 1..n_show).  Pre-history disappears here
# because we will pass first_rt_date 3 weeks in -- the design's burn-in handles
# the right-edge of the kernel automatically.
keep    <- (burn_in + 1L):n_total
dates   <- as.Date("2023-01-01") + 0:(n_show - 1L)
Y_obs   <- Y[keep]
Rt_show <- Rt_true[keep]
truth   <- tibble(date = dates, Rt_true = Rt_show)
cat(sprintf("simulated: %d days; peak exposures %.0f, peak observed Y %.0f\n",
            n_show, max(X[keep]), max(Y_obs)))

# ---- Retrospective fit ------------------------------------------------------
retro <- fit_girt_retrospective(
  obs_inc = Y_obs, dates = dates, g = g,
  mean_EY = mean_EY, sd_EY = sd_EY, severity = rho,
  first_rt_date         = dates[21],   # start Rt 3 weeks in
  likelihood_start_date = dates[35],   # let reporting delay fill
  knot_step             = 5L
)

# ---- Real-time fit at a mid-epidemic vintage --------------------------------
W  <- 90L                               # vintage = day 90 (clearly past first peak, declining)
vintage_date <- dates[W]
rt <- fit_girt_realtime(
  obs_inc = Y_obs[1:W], dates = dates[1:W], g = g,
  mean_EY = mean_EY, sd_EY = sd_EY, severity = rho,
  first_rt_date         = dates[21],
  likelihood_start_date = dates[35],
  knot_step             = 5L
)

# ---- Metrics ----------------------------------------------------------------
eval_retro <- retro$rt_df |> inner_join(truth, by = "date") |>
  filter(date >= dates[35], date <= dates[n_show - 3L])      # drop edge burn-in
cov95 <- mean(eval_retro$Rt_true >= eval_retro$Rt_lo &
              eval_retro$Rt_true <= eval_retro$Rt_hi)
mae   <- mean(abs(eval_retro$Rt_mean - eval_retro$Rt_true))
edge  <- tail(rt$rt_df$Rt_mean, 1L)
cat(sprintf("retro: MAE = %.3f  cov95 = %.2f  lambda = %.2e\n", mae, cov95, retro$lam))
cat(sprintf("realtime: vintage %s  edge Rt_hat = %.3f (truth %.3f)  gamma = %.2e\n",
            vintage_date, edge, Rt_true[burn_in + W], rt$gamma))

# ---- Plot -------------------------------------------------------------------
# Two stacked panels: Rt curves on top, observed cases on bottom.
rt_seg <- rt$rt_df |> filter(date >= vintage_date - 27, date <= vintage_date)
band   <- retro$rt_df |> filter(date >= dates[21])

p_rt <- ggplot() +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey55") +
  geom_ribbon(data = band, aes(date, ymin = Rt_lo, ymax = Rt_hi),
              fill = "#0072B2", alpha = 0.22) +
  geom_line(data = truth, aes(date, Rt_true), colour = "grey10", linewidth = 1.4) +
  geom_line(data = band,  aes(date, Rt_mean), colour = "#0072B2", linewidth = 1.1) +
  geom_line(data = rt_seg, aes(date, Rt_mean), colour = "#D55E00", linewidth = 1.4) +
  geom_point(data = tail(rt_seg, 1L), aes(date, Rt_mean), colour = "#D55E00", size = 2.6) +
  geom_vline(xintercept = vintage_date, linetype = "dashed", colour = "#D55E00", alpha = 0.5) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  labs(subtitle = sprintf(
        "black = truth, blue = retrospective (95%% band), orange = real-time nowcast (vintage %s)",
        format(vintage_date)),
       x = NULL, y = expression(R[t])) +
  theme_bw(base_size = 14) + theme(plot.subtitle = element_text(colour = "grey30"))

p_obs <- ggplot(tibble(date = dates, Y = Y_obs), aes(date, Y)) +
  geom_col(fill = "grey45", width = 0.9) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  labs(x = NULL, y = "Reported cases") +
  theme_bw(base_size = 14)

combined <- (p_rt / p_obs) + plot_layout(heights = c(2, 1)) +
  plot_annotation(title = "girt demo: two-wave R_t recovery from simulated cases",
                  theme = theme(plot.title = element_text(size = 17, face = "bold", hjust = 0.5)))

out_pdf <- file.path(girt_root, "examples", "demo_rt.pdf")
ggsave(out_pdf, combined, width = 11, height = 7)
cat(sprintf("wrote %s\n", out_pdf))
