# code/analyses/rtestim/plot_sim.R
#   -> figures/sim/rtestim/sim_combined_4x2.pdf
# Pure simulation viz (true Rt + infections/reports) from the cached sim_combined.rds.

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(patchwork); library(scales) })
sd <- dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))
if (!length(sd) || is.na(sd)) sd <- "code/analyses/rtestim"
source(file.path(sd, "_common.R"))
out_dir <- rt_figures_dir; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sc_levels <- c("Scenario 1: piecewise constant","Scenario 2: piecewise exponential",
               "Scenario 3: piecewise linear","Scenario 4: periodic")
sc_labels <- c("1. Piecewise constant","2. Piecewise exponential","3. Piecewise linear","4. Periodic")
sim <- readRDS(file.path(rt_data_dir, "sim_combined.rds")) |> mutate(scenario = factor(scenario, sc_levels, sc_labels))
long <- sim |> select(time, scenario, seir_incidence, seir_reports) |>
  pivot_longer(c(seir_incidence, seir_reports), names_to = "series", values_to = "count") |>
  mutate(series = factor(series, c("seir_incidence","seir_reports"),
                         c("True daily infections","Observed reported cases")))
pal <- c("True daily infections"="#0072B2","Observed reported cases"="#D55E00")
rt_df <- sim |> distinct(time, scenario, Rt)
base_theme <- theme_bw(base_size = 20) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(color="grey92", linewidth=0.3),
        panel.border = element_rect(color="grey60", linewidth=0.5), strip.background = element_blank(),
        plot.title.position = "plot", plot.margin = margin(4,8,4,8), panel.spacing.x = unit(1.2,"lines"),
        axis.text = element_text(size=18), axis.title = element_text(size=21))
p_rt <- ggplot(rt_df, aes(time, Rt)) +
  geom_hline(yintercept=1, color="grey55", linetype="22", linewidth=0.55) +
  geom_line(linewidth=0.9, color="grey15") + facet_wrap(~scenario, nrow=1) +
  scale_x_continuous(expand=expansion(mult=c(0.02,0.02)), breaks=c(0,150,300)) +
  scale_y_continuous(limits=c(0.3,2.85), breaks=c(0.5,1,1.5,2,2.5), expand=expansion(mult=c(0.02,0.05))) +
  labs(x=NULL, y=expression("True"~R[t])) + base_theme +
  theme(strip.text=element_text(face="bold", size=20, margin=margin(b=6)),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_text(margin=margin(r=8)))
p_counts <- ggplot(long, aes(time, pmax(count,1), color=series)) +
  geom_line(linewidth=0.7, alpha=0.9) + facet_wrap(~scenario, nrow=1) +
  scale_x_continuous(expand=expansion(mult=c(0.02,0.02)), breaks=c(0,150,300)) +
  scale_y_log10(labels=label_log(), breaks=10^(0:6), expand=expansion(mult=c(0.02,0.05))) +
  scale_color_manual(values=pal) + labs(x="Time (days)", y="Daily count", color=NULL) + base_theme +
  theme(strip.text=element_blank(), legend.position="bottom", legend.text=element_text(size=23),
        legend.key.width=unit(3.0,"lines"), legend.spacing.x=unit(0.6,"lines"), legend.margin=margin(t=8,b=2),
        axis.title.y=element_text(margin=margin(r=8))) +
  guides(color=guide_legend(override.aes=list(linewidth=1.7)))
ggsave(file.path(out_dir, "sim_combined_4x2.pdf"), (p_rt / p_counts) + plot_layout(heights=c(1,1.6)), width=16, height=8.5)
cat("wrote sim_combined_4x2.pdf\n")
