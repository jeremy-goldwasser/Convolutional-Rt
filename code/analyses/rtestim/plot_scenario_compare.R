# code/analyses/rtestim/plot_scenario_compare.R
#   -> figures/sim/rtestim/scenario2_compare.pdf
# Scenario-2 method comparison (zoom), girt MechRt + cached baselines + EpiNow2 (GP). No jumps for S2.

suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(tibble) })
sd <- dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))
if (!length(sd) || is.na(sd)) sd <- "code/analyses/rtestim"
source(file.path(sd, "_common.R"))
gtree <- file.path(rt_results_dir, "girt"); fit_dir <- rt_results_dir
plot_dir <- rt_figures_dir; dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
date0 <- rt_cfg$date0; z95 <- qnorm(0.975)

sc_name <- "Scenario 2: piecewise exponential"; tag <- "fake__Scenario_2_piecewise_exponential"
x_lo <- 70L; x_hi <- 130L; epinow2_dir <- "epinow2_s2"; epinow2_label <- "EpiNow2 GP"
method_colors <- c(Truth="black", MechRt="#2D6CA2", `MechRt (jumps)`="#8FBEDB", EpiNow2="#5D3A1E",
                   estimateR="#D55E00", EpiEstim="#009E73", rtestim="#CC79A7", EpiLPS="#882255")
method_ltys   <- c(Truth="solid", MechRt="22", `MechRt (jumps)`="22", EpiNow2="11",
                   estimateR="solid", EpiEstim="solid", rtestim="solid", EpiLPS="solid")
fg_methods <- c("MechRt","MechRt (jumps)","EpiNow2")

rt_cache <- readRDS(file.path(gtree, "retro_rt_estimates.rds")); sim <- readRDS(file.path(rt_data_dir, "sim_combined.rds"))
rt_other <- rt_cache |> filter(method != "EpiNow2 GP", scenario == sc_name, day >= x_lo, day <= x_hi) |> mutate(method = as.character(method))
ef <- readRDS(file.path(fit_dir, epinow2_dir, sprintf("fit_epinow2_%s__summary.rds", tag)))$rt_df
if ("type" %in% names(ef)) ef <- ef |> filter(type == "estimate")
rt_epi <- ef |> mutate(date = as.Date(date), day = as.integer(date - date0) + 1L,
    Rt_lo = if ("Rt_lo_90" %in% names(ef)) Rt_lo_90 else Rt_mean - z95*Rt_sd,
    Rt_hi = if ("Rt_hi_90" %in% names(ef)) Rt_hi_90 else Rt_mean + z95*Rt_sd) |>
  filter(day >= x_lo, day <= x_hi) |> transmute(method="EpiNow2", scenario=sc_name, si_type="fake", date, Rt_mean, Rt_lo, Rt_hi, day)
rt_all <- bind_rows(rt_other, rt_epi)
truth <- sim |> filter(scenario == sc_name, si_type == "fake") |> transmute(day = as.integer(time), Rt_truth = Rt) |> filter(day >= x_lo, day <= x_hi)
method_order <- c("Truth","MechRt","EpiNow2","EpiEstim","estimateR","EpiLPS","rtestim")
rt_all$method <- factor(rt_all$method, levels = method_order[-1L])
truth_lab <- truth |> mutate(method = factor("Truth", levels = method_order))
rt_bg <- rt_all |> filter(!(as.character(method) %in% fg_methods)); rt_fg <- rt_all |> filter(as.character(method) %in% fg_methods)
legend_w <- c(Truth=1.3, MechRt=1.0, `MechRt (jumps)`=1.0, EpiNow2=0.85, EpiEstim=0.6, estimateR=0.6, EpiLPS=0.6, rtestim=0.6)
legend_labels <- setNames(method_order, method_order); legend_labels["EpiNow2"] <- epinow2_label

p <- ggplot() +
  geom_hline(yintercept=1, linetype="dashed", colour="grey50") +
  geom_line(data=truth_lab, aes(day, Rt_truth, color=method, linetype=method), linewidth=1.3) +
  geom_line(data=rt_bg, aes(day, Rt_mean, color=method, group=method, linetype=method), linewidth=0.7, alpha=0.8) +
  geom_line(data=rt_fg |> filter(method=="EpiNow2"), aes(day, Rt_mean, color=method, group=method, linetype=method), linewidth=1.0, alpha=0.95) +
  geom_line(data=rt_fg |> filter(method!="EpiNow2"), aes(day, Rt_mean, color=method, group=method, linetype=method), linewidth=1.1, alpha=0.95) +
  scale_color_manual(values=method_colors, name=NULL, breaks=method_order, labels=legend_labels,
    guide=guide_legend(override.aes=list(linewidth=unname(legend_w[method_order]), linetype=unname(method_ltys[method_order]), alpha=1))) +
  scale_linetype_manual(values=method_ltys, name=NULL, breaks=method_order, labels=legend_labels, guide="none") +
  coord_cartesian(xlim=c(x_lo, x_hi)) + labs(title=sc_name, x="Day", y=expression(R[t])) +
  theme_minimal(base_size=13) +
  theme(legend.position="right", legend.direction="vertical", plot.title=element_text(size=18, hjust=0.5),
        legend.text=element_text(size=14), legend.key.width=unit(1.2,"cm"))
ggsave(file.path(plot_dir, "scenario2_compare.pdf"), p, width=9.5, height=2.8, dpi=300)
cat("wrote scenario2_compare.pdf\n")
