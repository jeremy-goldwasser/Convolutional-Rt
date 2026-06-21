# code/analyses/sim_flu/plot_rt_evolution.R
#   -> figures/sim/flu/wiggly/rt_evolution_raw.pdf
# Appendix figure: EpiNow2 GP posterior-mean Rt over MCMC iterations (wiggly).
# Reads snapshot poller output from results/sim/flu/wiggly/epinow2_snapshots/<run>/snapshots/.
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(viridis) })
sd_ <- dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))
if (!length(sd_) || is.na(sd_)) sd_ <- "code/analyses/sim_flu"
source(file.path(sd_, "_common.R"))
p <- sim_flu_paths("wiggly")
dir.create(p$fig_dir, recursive = TRUE, showWarnings = FALSE)

run_dirs <- list.dirs(p$snap_root, recursive = FALSE)
raw_runs <- run_dirs[grepl("__raw$", basename(run_dirs))]
if (!length(raw_runs)) stop("no __raw snapshot run under ", p$snap_root)
run_dir <- raw_runs[order(basename(raw_runs), decreasing = TRUE)][1L]
cat("using snapshot run:", basename(run_dir), "\n")

idx <- readRDS(file.path(run_dir, "rt_index_map.rds"))
snap_files <- sort(list.files(file.path(run_dir, "snapshots"), pattern = "^snap_.*\\.rds$", full.names = TRUE))
if (!length(snap_files)) stop("no snapshots in ", run_dir)

rows <- lapply(snap_files, function(f) {
  s <- readRDS(f)
  data.frame(snapshot_id = as.integer(sub(".*snap_(\\d+)_.*", "\\1", basename(f))),
             n_draws_total = s$n_draws_total, elapsed_sec = s$elapsed_sec,
             stan_index = s$stan_index, Rt_mean = unname(s$mean), stringsAsFactors = FALSE)
})
df <- do.call(rbind, rows) |>
  dplyr::left_join(idx[, c("stan_index","date","type")], by = "stan_index") |>
  dplyr::filter(type %in% c("estimate", "estimate based on partial data")) |>
  dplyr::mutate(date = as.Date(date))

snap_lookup <- df |> dplyr::distinct(snapshot_id, n_draws_total, elapsed_sec) |> dplyr::arrange(n_draws_total)
iter_breaks <- c(min(snap_lookup$n_draws_total), pretty(snap_lookup$n_draws_total, n = 4))
iter_breaks <- unique(iter_breaks[iter_breaks >= min(snap_lookup$n_draws_total) & iter_breaks <= max(snap_lookup$n_draws_total)])

p_fig <- ggplot(df, aes(x = date, y = Rt_mean, group = snapshot_id, color = n_draws_total)) +
  geom_line(alpha = 0.55, linewidth = 0.4) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  scale_color_viridis_c(name = "MCMC iterations", option = "viridis",
                        breaks = iter_breaks, labels = as.character(iter_breaks),
                        guide = guide_colorbar(barwidth = 12, barheight = 0.6,
                                               title.position = "left", title.vjust = 0.9)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  labs(title = expression(paste("EpiNow2 ", R[t], ": evolution during MCMC")), x = NULL, y = expression(R[t])) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size = 17), axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16), legend.position = "bottom", legend.box = "horizontal",
        legend.title = element_text(size = 15), legend.text = element_text(size = 13))
ggsave(file.path(p$fig_dir, "rt_evolution_raw.pdf"), p_fig, width = 7, height = 3.6, device = "pdf")
cat("wrote rt_evolution_raw.pdf\n")
