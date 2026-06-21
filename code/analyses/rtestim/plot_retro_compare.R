# code/analyses/rtestim/plot_retro_compare.R
#   -> figures/sim/rtestim/retro_compare_paper.pdf
# Paper retrospective comparison: girt MechRt (+jumps S1-S3) vs EpiNow2/EpiEstim/estimateR.
# MechRt + baselines from the girt bundle; EpiNow2 from the canonical per-scenario kernel.

suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(tibble); library(scales) })
sd <- dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))
if (!length(sd) || is.na(sd)) sd <- "code/analyses/rtestim"
source(file.path(sd, "_common.R"))

gtree    <- file.path(rt_results_dir, "girt")             # girt MechRt bundle + jumps
fit_dir  <- rt_results_dir                                # canonical EpiNow2 + baseline fits
plot_dir <- rt_figures_dir; dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
date0 <- rt_cfg$date0; burn_in <- 30L; tail_drop <- 10L; n_days <- 300L; eval_hi <- n_days - tail_drop + 1L

rt_all <- readRDS(file.path(gtree, "retro_rt_estimates.rds")) |>
  filter(method != "EpiNow2 GP") |> mutate(method = droplevels(factor(method))) |> filter(day >= burn_in)
sim <- readRDS(file.path(rt_data_dir, "sim_combined.rds"))
truth_plot <- sim |> distinct(scenario, si_type, time, Rt) |> transmute(scenario, si_type, day = time, Rt_truth = Rt) |> filter(day >= burn_in)

en2_dir_for <- c("Scenario 1: piecewise constant"="epinow2_grw","Scenario 2: piecewise exponential"="epinow2_s2",
                 "Scenario 3: piecewise linear"="epinow2_grw","Scenario 4: periodic"="epinow2_s4")
en2_canon <- bind_rows(lapply(unique(rt_all$scenario), function(sc) {
  tag <- paste0("fake__", gsub("[:,]", "", gsub(" ", "_", sc)))
  d <- readRDS(file.path(fit_dir, en2_dir_for[[sc]], sprintf("fit_epinow2_%s__summary.rds", tag)))$rt_df
  tibble(method="EpiNow2", scenario=sc, si_type="fake", date=d$date, Rt_mean=d$Rt_mean,
         Rt_lo=d$Rt_lo_90, Rt_hi=d$Rt_hi_90, day=as.integer(d$date - date0) + 1L) }))
jtags <- c("Scenario 1: piecewise constant"="fake__Scenario_1_piecewise_constant",
           "Scenario 2: piecewise exponential"="fake__Scenario_2_piecewise_exponential",
           "Scenario 3: piecewise linear"="fake__Scenario_3_piecewise_linear")
mechrt_jumps <- bind_rows(lapply(names(jtags), function(sc) {
  f <- file.path(gtree, sprintf("fit_ours_jump_%s.rds", jtags[[sc]])); if (!file.exists(f)) return(NULL)
  rt <- readRDS(f)$rt
  tibble(method="MechRt (jumps)", scenario=sc, si_type="fake", date=rt$date, Rt_mean=rt$Rt_mean,
         Rt_lo=rt$Rt_lo, Rt_hi=rt$Rt_hi, day=as.integer(rt$date - date0) + 1L) }))

paper_methods <- c("Truth","MechRt","MechRt (jumps)","EpiNow2","EpiEstim","estimateR")
sc_levels <- c("Scenario 1: piecewise constant","Scenario 2: piecewise exponential","Scenario 3: piecewise linear","Scenario 4: periodic")
sc_labels <- c("1. Piecewise constant","2. Piecewise exponential","3. Piecewise linear","4. Periodic")
ylim_paper <- tibble(scenario=sc_levels, ymin=c(0.65,0.50,0.45,0.35), ymax=c(2.20,2.85,2.55,2.45))

rt_paper <- bind_rows(
  rt_all |> filter(as.character(method) %in% c("MechRt","EpiEstim","estimateR")) |> select(method,scenario,si_type,date,Rt_mean,Rt_lo,Rt_hi,day),
  mechrt_jumps, en2_canon) |>
  filter(day >= burn_in, day <= eval_hi) |> left_join(ylim_paper, by="scenario") |>
  mutate(Rt_mean = pmin(pmax(Rt_mean, ymin), ymax), method = factor(method, paper_methods),
         scenario = factor(scenario, sc_levels, sc_labels))
truth_paper <- truth_plot |> filter(day >= burn_in, day <= eval_hi) |>
  mutate(scenario = factor(scenario, sc_levels, sc_labels), method = factor("Truth", paper_methods))
ylim_blank <- bind_rows(ylim_paper |> transmute(scenario, day=burn_in, Rt_mean=ymin),
                        ylim_paper |> transmute(scenario, day=burn_in, Rt_mean=ymax)) |>
  mutate(scenario = factor(scenario, sc_levels, sc_labels))
paper_colors <- c(Truth="black", MechRt="#2D6CA2", `MechRt (jumps)`="#8FBEDB", EpiNow2="#5D3A1E", EpiEstim="#009E73", estimateR="#D55E00")
paper_ltys   <- c(Truth="solid", MechRt="22", `MechRt (jumps)`="62", EpiNow2="11", EpiEstim="solid", estimateR="solid")
legend_w     <- c(Truth=1.5, MechRt=1.3, `MechRt (jumps)`=1.3, EpiNow2=1.2, EpiEstim=1.05, estimateR=1.05)
fg <- c("MechRt","MechRt (jumps)","EpiNow2")
rt_bg <- rt_paper |> filter(!(as.character(method) %in% fg)); rt_fg <- rt_paper |> filter(as.character(method) %in% fg)

p <- ggplot() +
  geom_hline(yintercept=1, linetype="dashed", colour="grey55", linewidth=0.5) +
  geom_blank(data=ylim_blank, aes(day, Rt_mean)) +
  geom_line(data=truth_paper, aes(day, Rt_truth, colour=method, linetype=method), linewidth=1.5) +
  geom_line(data=rt_bg, aes(day, Rt_mean, colour=method, linetype=method, group=method), linewidth=1.05, alpha=0.95) +
  geom_line(data=rt_fg |> filter(method=="EpiNow2"), aes(day, Rt_mean, colour=method, linetype=method, group=method), linewidth=1.2, alpha=0.95) +
  geom_line(data=rt_fg |> filter(method!="EpiNow2"), aes(day, Rt_mean, colour=method, linetype=method, group=method), linewidth=1.3, alpha=0.95) +
  facet_wrap(~scenario, scales="free_y", ncol=2) +
  scale_colour_manual(values=paper_colors, name=NULL, breaks=paper_methods,
    guide=guide_legend(override.aes=list(linewidth=unname(legend_w[paper_methods]), linetype=unname(paper_ltys[paper_methods]), alpha=1))) +
  scale_linetype_manual(values=paper_ltys, name=NULL, breaks=paper_methods, guide="none") +
  scale_x_continuous(expand=expansion(mult=c(0.01,0.02)), breaks=c(50,150,250)) +
  scale_y_continuous(breaks=scales::pretty_breaks(n=4)) + coord_cartesian(xlim=c(burn_in, eval_hi)) +
  labs(x="Time (days)", y=expression(R[t])) + theme_bw(base_size=20) +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_line(colour="grey92", linewidth=0.3),
        panel.border=element_rect(colour="grey60", linewidth=0.5), panel.spacing=unit(1.0,"lines"),
        strip.background=element_blank(), strip.text=element_text(size=20, margin=margin(b=6)),
        axis.text=element_text(size=18), axis.title=element_text(size=21), axis.title.y=element_text(margin=margin(r=8)),
        legend.position="bottom", legend.text=element_text(size=22), legend.key.width=unit(2.8,"lines"),
        legend.spacing.x=unit(0.5,"lines"), legend.margin=margin(t=8,b=2), plot.margin=margin(4,8,4,8)) +
  guides(colour=guide_legend(nrow=1, byrow=TRUE, override.aes=list(linewidth=unname(legend_w[paper_methods]), linetype=unname(paper_ltys[paper_methods]), alpha=1)))
ggsave(file.path(plot_dir, "retro_compare_paper.pdf"), p, width=14, height=9)
cat("wrote retro_compare_paper.pdf\n")
