# code/analyses/rtestim/table.R
#   -> results/sim/rtestim/{seir_full_compare_overdisp,seir_final_table}.csv
# Final rtestim performance table (MAE | CE, x1e-2), MechRt via girt (log, overdisp=TRUE).
# MechRt smooth/jump: girt re-extract at CV lam_min/lam_1se, overdispersion=TRUE.
# EpiNow2 (per-scenario canonical kernel) + estimateR/EpiEstim/rtestim/EpiLPS: cached.
# BurnIn=50 window; per-scenario rule = 1se (S1-S3), min (S4).

sd <- dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))
if (!length(sd) || is.na(sd)) sd <- "code/analyses/rtestim"
source(file.path(sd, "_common.R"))
suppressPackageStartupMessages({ library(tidyr); library(readr) })

z95 <- qnorm(0.975); z90 <- qnorm(0.95); alpha_grid <- c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
windows <- list("30" = c(30L, 291L), "50" = c(50L, 291L))
rt_results <- rt_results_dir
out_dir <- rt_results_dir; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
tag_for <- function(sc) sprintf("fake__%s", gsub(":? ", "_", sc))
scenario_jumps <- list("Scenario 1: piecewise constant"=120L, "Scenario 2: piecewise exponential"=100L,
                       "Scenario 3: piecewise linear"=c(75L,150L,225L))
rule_by_sc <- c("Scenario 1: piecewise constant"="1se","Scenario 2: piecewise exponential"="1se",
                "Scenario 3: piecewise linear"="1se","Scenario 4: periodic"="min")
epinow2_dir <- c("Scenario 1: piecewise constant"="epinow2_grw","Scenario 2: piecewise exponential"="epinow2_s2",
                 "Scenario 3: piecewise linear"="epinow2_grw","Scenario 4: periodic"="epinow2_s4")
epinow2_kernel <- c("Scenario 1: piecewise constant"="GRW","Scenario 2: piecewise exponential"="GP",
                    "Scenario 3: piecewise linear"="GRW","Scenario 4: periodic"="GP")

emp_cov <- function(rt) { rt <- rt[is.finite(rt$Rt_truth) & is.finite(rt$Rt_mean) & is.finite(rt$sigma), ]
  vapply(alpha_grid, function(a) { z <- qnorm(1 - (1 - a)/2)
    mean(rt$Rt_truth >= rt$Rt_mean - z*rt$sigma & rt$Rt_truth <= rt$Rt_mean + z*rt$sigma) }, numeric(1)) }
eval_band <- function(band, truth_sc, win_lo, win_hi) {
  band$sigma <- (band$Rt_hi - band$Rt_lo) / (2 * z95)
  ev <- merge(band, truth_sc, by = "date")
  ev <- ev[ev$date >= rt_cfg$date0 + (win_lo - 1L) & ev$date <= rt_cfg$date0 + (win_hi - 1L), ]
  emp <- emp_cov(ev)
  data.frame(MAE = mean(abs(ev$Rt_mean - ev$Rt_truth)), CE = mean(abs(emp - alpha_grid)))
}

sim <- readRDS(file.path(rt_data_dir, "sim_combined.rds")); sim$date <- rt_cfg$date0 + (sim$time - 1L)
dl  <- rt_delays("fake"); scenarios <- unique(sim$scenario)
epinow2_band <- function(sc) { r <- readRDS(file.path(rt_results, epinow2_dir[[sc]], sprintf("fit_epinow2_%s__summary.rds", tag_for(sc))))$rt_df
  data.frame(date = as.Date(r$date), Rt_mean = r$Rt_mean, Rt_lo = r$Rt_mean - z95*r$Rt_sd, Rt_hi = r$Rt_mean + z95*r$Rt_sd) }

# ---- overdisp comparison (girt MechRt smooth/jump at min/1se, od=TRUE) -------
rows <- list()
for (sc in scenarios) {
  d <- sim[sim$scenario == sc & sim$si_type == "fake", ]; d <- d[order(d$time), ]
  truth_sc <- data.frame(date = d$date, Rt_truth = d$Rt); y <- as.numeric(d$seir_reports); dates <- d$date
  dc_s <- rt_design_cv(y, dates, dl)
  bs_min <- rt_band(dc_s, dc_s$lam_min, TRUE); bs_1se <- rt_band(dc_s, dc_s$lam_1se, TRUE)
  jt <- scenario_jumps[[sc]]; bj_min <- bj_1se <- NULL
  if (!is.null(jt)) { dc_j <- rt_design_cv(y, dates, dl, jump_times = jt)
    bj_min <- rt_band(dc_j, dc_j$lam_min, TRUE); bj_1se <- rt_band(dc_j, dc_j$lam_1se, TRUE) }
  eb <- epinow2_band(sc)
  for (b in names(windows)) { w <- windows[[b]]
    add <- function(method, variant, rule, band) rows[[length(rows)+1L]] <<- cbind(
      data.frame(Scenario=sc, Method=method, Variant=variant, Rule=rule, BurnIn=b), eval_band(band, truth_sc, w[1], w[2]))
    add("MechRt","smooth","min", bs_min); add("MechRt","smooth","1se", bs_1se)
    if (!is.null(jt)) { add("MechRt","jump","min", bj_min); add("MechRt","jump","1se", bj_1se) }
    add(sprintf("EpiNow2 (%s)", epinow2_kernel[[sc]]), NA, NA, eb)
  }
}
od <- do.call(rbind, rows); od$MAE <- round(100*od$MAE, 2); od$CE <- round(100*od$CE, 2)
write_csv(od, file.path(out_dir, "seir_full_compare_overdisp.csv"))

# ---- final table: MechRt/jump/EpiNow2 (BurnIn50, per-rule) + cached baselines ----
mech <- do.call(rbind, lapply(scenarios, function(sc) {
  s <- od[od$BurnIn=="50" & od$Method=="MechRt" & od$Scenario==sc & od$Rule==rule_by_sc[[sc]], ]
  data.frame(Row = ifelse(s$Variant=="jump","MechRt + jump","MechRt"), Scenario=sc, MAE=s$MAE, CE=s$CE) }))
en2 <- do.call(rbind, lapply(scenarios, function(sc) {
  s <- od[od$BurnIn=="50" & grepl("^EpiNow2", od$Method) & od$Scenario==sc, ]
  data.frame(Row="EpiNow2", Scenario=sc, MAE=s$MAE, CE=s$CE) }))
truth_all <- do.call(rbind, lapply(scenarios, function(sc){ d<-sim[sim$scenario==sc & sim$si_type=="fake",]; data.frame(date=d$date, Rt_truth=d$Rt, scenario=sc)}))
load_other <- function(lbl, tmpl) do.call(rbind, lapply(scenarios, function(sc) {
  x <- readRDS(file.path(rt_results, sprintf(tmpl, tag_for(sc)))); rt <- if (is.data.frame(x)) x else x$rt
  band <- data.frame(date=as.Date(rt$date), Rt_mean=rt$Rt_mean, Rt_lo=rt$Rt_lo, Rt_hi=rt$Rt_hi)
  m <- eval_band(band, truth_all[truth_all$scenario==sc, c("date","Rt_truth")], 50L, 291L)
  data.frame(Row=lbl, Scenario=sc, MAE=100*m$MAE, CE=100*m$CE) }))
others <- rbind(load_other("estimateR","fit_estimater_%s.rds"), load_other("EpiEstim","fit_epiestim_%s.rds"),
                load_other("rtestim","fit_rtestim_%s.rds"), load_other("EpiLPS","fit_epilps_%s.rds"))
all_df <- rbind(mech, others, en2)
sc_short <- c("Scenario 1: piecewise constant"="S1","Scenario 2: piecewise exponential"="S2",
              "Scenario 3: piecewise linear"="S3","Scenario 4: periodic"="S4")
wide <- all_df |> mutate(Scen = sc_short[Scenario]) |>
  pivot_longer(c(MAE, CE), names_to="Metric", values_to="Value") |>
  unite("col", Scen, Metric, sep="_") |> pivot_wider(id_cols=Row, names_from=col, values_from=Value) |>
  mutate(Row = factor(Row, levels=c("MechRt","MechRt + jump","EpiNow2","estimateR","EpiEstim","rtestim","EpiLPS"))) |>
  arrange(Row) |> select(Row, S1_MAE,S1_CE,S2_MAE,S2_CE,S3_MAE,S3_CE,S4_MAE,S4_CE)
write_csv(wide, file.path(out_dir, "seir_final_table.csv"))
cat("wrote", file.path(out_dir, "seir_final_table.csv"), "\n")
print(as.data.frame(wide), row.names = FALSE, digits = 3)
