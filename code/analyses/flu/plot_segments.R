# code/analyses/flu/plot_segments.R
#   -> figures/real/flu/s2/realtime_segments_with_ci_and_epinow2.pdf
# Real-time s2 nowcast segments with split-conformal CIs (calibrated from the
# girt daily vintages) + EpiNow2 overlay.  MechRt = girt (log); EpiNow2 canonical.

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(scales) })
sd <- dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))
if (!length(sd) || is.na(sd)) sd <- "code/analyses/flu"
source(file.path(sd, "_common.R"))

CI_LEVEL <- 0.95; BUFFER_DAYS <- 14L; D_MAX_C <- 6L; Z95 <- qnorm(0.975)
END_VINTAGE <- as.Date("2024-01-28")
COL_RETRO <- "#111111"; COL_RT_LINE <- "#2D6CA2"; COL_INFL <- "#2D6CA2"
COL_EPI_LINE <- "#5D3A1E"; COL_EPI_FILL <- "#5D3A1E"; ALPHA_DISP <- 1 - CI_LEVEL

s2_g  <- file.path(flu_results_dir, "girt", "s2")    # MechRt (girt)
s2_c  <- file.path(flu_results_dir, "s2")            # EpiNow2 (canonical)

# daily girt vintages (meta$week_end_date is the vintage date)
all_vints <- do.call(rbind, lapply(
  list.files(file.path(s2_g, "real_time", "daily_fits"), pattern = "^daily_.*\\.rds$", full.names = TRUE),
  function(f) { r <- readRDS(f); rt <- r$rt_df
    data.frame(vintage = as.Date(r$meta$week_end_date), date = as.Date(rt$date),
               Rt_mean = rt$Rt_mean, Rt_lo = rt$Rt_lo, Rt_hi = rt$Rt_hi) })) |>
  mutate(sigma = (Rt_hi - Rt_lo) / (2 * Z95))
cat(sprintf("Loaded %d daily girt vintages (%s -> %s)\n",
            length(unique(all_vints$vintage)), min(all_vints$vintage), max(all_vints$vintage)))

s2_retro <- { r <- readRDS(list.files(file.path(s2_g, "retrospective", "results"), full.names = TRUE)[1])$rt_df
  data.frame(date = as.Date(r$date), Rt_mean = r$Rt_mean) }

slice_seg <- function(f) { r <- readRDS(f); rt <- r$rt_df
  data.frame(date = as.Date(rt$date), Rt_mean = rt$Rt_mean, Rt_lo = rt$Rt_lo, Rt_hi = rt$Rt_hi,
             week_end = as.Date(r$meta$week_end_date)) }
s2_raw <- do.call(rbind, lapply(
  list.files(file.path(s2_g, "real_time", "results"), pattern = "__natural_linear_taper_cv_1se\\.rds$", full.names = TRUE),
  slice_seg)) |>
  mutate(d_to_edge = as.integer(week_end - date)) |>
  filter(d_to_edge >= 0, d_to_edge <= D_MAX_C, week_end <= END_VINTAGE)
weekly_W <- sort(unique(s2_raw$week_end))

.cq <- function(s, alpha) { s <- s[is.finite(s)]; n <- length(s); if (!n) return(NA_real_)
  k <- ceiling((1 - alpha) * (n + 1)); if (k > n) return(NA_real_); sort(s)[k] }
conformal_q_at_W <- function(W) {
  rtT <- all_vints |> filter(vintage == W) |> transmute(date, Rt_anchor = Rt_mean)
  if (!nrow(rtT)) return(NULL)
  cal <- all_vints |> filter(vintage < W) |> mutate(d_to_edge = as.integer(vintage - date)) |>
    filter(d_to_edge >= 0L, d_to_edge <= D_MAX_C, date <= W - BUFFER_DAYS) |>
    inner_join(rtT, by = "date") |> transmute(d_to_edge, s_abs = abs(Rt_mean - Rt_anchor))
  if (!nrow(cal)) return(NULL)
  cal |> group_by(d_to_edge) |> summarise(q_abs = .cq(s_abs, ALPHA_DISP), .groups = "drop") |> mutate(week_end = W)
}
q_by_W <- bind_rows(lapply(weekly_W, conformal_q_at_W))

s2_segs <- s2_raw |> left_join(q_by_W, by = c("week_end", "d_to_edge")) |>
  mutate(Rt_lo_inf = Rt_mean - q_abs, Rt_hi_inf = Rt_mean + q_abs)
x_min <- min(s2_segs$date); x_max <- max(s2_segs$date)
s2_retro_win <- s2_retro[s2_retro$date >= x_min & s2_retro$date <= x_max, ]

cov_df <- s2_segs |> left_join(s2_retro |> rename(Rt_retro = Rt_mean), by = "date") |>
  filter(!is.na(Rt_retro), is.finite(Rt_lo_inf)) |>
  mutate(in_inf = Rt_retro >= Rt_lo_inf & Rt_retro <= Rt_hi_inf)
cov_infl <- mean(cov_df$in_inf); n_conf <- nrow(cov_df)

z_disp <- qnorm(0.5 + CI_LEVEL / 2)
ep_segs <- do.call(rbind, lapply(
  list.files(file.path(s2_c, "real_time", "epinow2", "results"), pattern = "__gp__summary\\.rds$", full.names = TRUE),
  function(f) { r <- readRDS(f); rt <- r$rt_df
    data.frame(week_end = as.Date(r$meta$week_end_date), date = as.Date(rt$date), Rt_mean = rt$Rt_mean, Rt_sd = rt$Rt_sd) })) |>
  mutate(d_to_edge = as.integer(week_end - date)) |>
  filter(d_to_edge >= 0, d_to_edge <= D_MAX_C, week_end <= END_VINTAGE) |>
  mutate(Rt_lo = Rt_mean - z_disp * Rt_sd, Rt_hi = Rt_mean + z_disp * Rt_sd)
conf_vints <- s2_segs |> filter(is.finite(Rt_lo_inf)) |> pull(week_end) |> unique()
s2_segs_c <- s2_segs |> filter(week_end %in% conf_vints); ep_segs_c <- ep_segs |> filter(week_end %in% conf_vints)
ep_retro <- readRDS(list.files(file.path(s2_c, "retrospective", "epinow2", "results"),
  pattern = "__gp__summary\\.rds$", full.names = TRUE)[[1]])$rt_df |> transmute(date = as.Date(date), Rt_retro = Rt_mean)
cov_ep <- ep_segs_c |> inner_join(ep_retro, by = "date") |> summarise(cov = mean(Rt_retro >= Rt_lo & Rt_retro <= Rt_hi), n = dplyr::n())

band_line_key <- function(a) function(data, params, size)
  grid::grobTree(grid::rectGrob(gp = grid::gpar(col = NA, fill = scales::alpha(data$colour, a))),
                 ggplot2::draw_key_path(data, params, size))
p <- ggplot() +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_ribbon(data = ep_segs_c, aes(date, ymin = Rt_lo, ymax = Rt_hi, group = week_end), fill = COL_EPI_FILL, alpha = 0.18) +
  geom_ribbon(data = filter(s2_segs_c, is.finite(Rt_lo_inf)), aes(date, ymin = Rt_lo_inf, ymax = Rt_hi_inf, group = week_end), fill = COL_INFL, alpha = 0.30) +
  geom_line(data = s2_retro_win, aes(date, Rt_mean, colour = "Retrospective"), linewidth = 1.0) +
  geom_line(data = ep_segs_c, aes(date, Rt_mean, group = week_end, colour = "EpiNow2"), linewidth = 0.9, key_glyph = band_line_key(0.18)) +
  geom_point(data = filter(ep_segs_c, d_to_edge == 0L), aes(date, Rt_mean), colour = COL_EPI_LINE, size = 1.5, show.legend = FALSE) +
  geom_line(data = s2_segs_c, aes(date, Rt_mean, group = week_end, colour = "MechRt"), linewidth = 0.9, key_glyph = band_line_key(0.30)) +
  geom_point(data = filter(s2_segs_c, d_to_edge == 0L), aes(date, Rt_mean), colour = COL_RT_LINE, size = 1.5, show.legend = FALSE) +
  scale_colour_manual(values = c("MechRt" = COL_RT_LINE, "EpiNow2" = COL_EPI_LINE, "Retrospective" = COL_RETRO),
                      name = NULL, breaks = c("MechRt", "EpiNow2", "Retrospective")) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.2))) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d", minor_breaks = NULL) +
  labs(title = "Real-time nowcasts", x = NULL, y = expression(R[t])) +
  theme_minimal(base_size = 20) +
  theme(plot.title = element_text(hjust = 0, size = 20, margin = margin(b = 8)), axis.title.y = element_text(size = 22),
        axis.text = element_text(size = 18), legend.text = element_text(size = 20), legend.key.width = grid::unit(1.6, "cm"),
        panel.grid.minor = element_blank(), legend.position = "bottom")
out <- file.path(flu_figures_dir, "s2", "realtime_segments_with_ci_and_epinow2.pdf")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
ggsave(out, p, width = 13, height = 6.0, dpi = 300)
cat(sprintf("wrote %s  (MechRt conformal cov=%.0f%% n=%d ; EpiNow2 cov=%.0f%% n=%d)\n",
            out, 100 * cov_infl, n_conf, 100 * cov_ep$cov, cov_ep$n))
