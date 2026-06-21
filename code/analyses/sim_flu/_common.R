# code/analyses/sim_flu/_common.R
#
# Paths + shared style for the sim-flu deliverables.  Tree layout (post-restructure):
#   data/sim/flu/{wiggly,noisy_smooth}/tuning_sim_results.rds + supporting files
#   results/sim/flu/{wiggly,noisy_smooth}/{retro/,realtime/,epinow2_{retro,realtime}/,other_methods_realtime/}
#     wiggly only also has epinow2_snapshots/<run>/snapshots/, table/
#   figures/sim/flu/{wiggly,noisy_smooth}/
#
# Both scenarios use IDENTITY deconvolution (mechrt default, verified from source).
# girt engine matches mechrt fits byte-for-byte under the same convention; see
# (verify.R was here pre-flatten; deleted).  Plot scripts here READ the canonical
# cached weekly/retro fits (whose Rt_lo/Rt_hi are already post-conformal for
# real-time fits, post-Wald for retro), so they do NOT re-fit.
#
# GI convention (for any re-fit):  g = zeta^EI/sum(zeta^EI), no lag (sim never used a lag-1 pre_I);
# effective mip = sum(zeta^EI) = sum(P_still_I).  Same compartmental params as real flu
# (mean_lat 2.0/sd 1.2, mean_inf 2.75/sd 1.0, mean_EH 5.7/sd 2.3, severity 0.015,
# knot_step 5).  No DoW.  Retro cv_select_rule = "min"; realtime cv_min + FV-min gamma + taper.

.find_repo_root <- function() {
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) tryCatch(dirname(normalizePath(f[1])), error = function(e) getwd()) else getwd()
  for (i in 1:9) { if (all(c("data","results","figures","code") %in% list.files(d))) return(normalizePath(d)); p <- dirname(d); if (p == d) break; d <- p }
  stop("repo root not found")
}
repo_root   <- .find_repo_root()
data_dir    <- file.path(repo_root, "data")
results_dir <- file.path(repo_root, "results")
figures_dir <- file.path(repo_root, "figures")

# scenario keys: "wiggly" (= flu3, variant rt_raw) and "noisy_smooth" (= flu4, variant rt_loess)
sim_flu_paths <- function(scenario = c("wiggly", "noisy_smooth")) {
  scenario <- match.arg(scenario)
  base_r <- file.path(results_dir, "sim", "flu", scenario)
  list(
    scenario  = scenario,
    truth_col = if (scenario == "wiggly") "rt_raw" else "rt_loess",
    data_dir  = file.path(data_dir,    "sim", "flu", scenario),
    sim_rds   = file.path(data_dir,    "sim", "flu", scenario, "tuning_sim_results.rds"),
    retro_dir = file.path(base_r, "retro"),
    weekly_dir = file.path(base_r, "realtime", "weekly"),
    manifest  = file.path(base_r, "realtime", "manifest.csv"),
    ep_retro  = file.path(base_r, "epinow2_retro"),
    ep_rt     = file.path(base_r, "epinow2_realtime"),
    other_rt  = file.path(base_r, "other_methods_realtime"),
    snap_root = file.path(base_r, "epinow2_snapshots"),  # wiggly only
    table_dir = file.path(base_r, "table"),
    fig_dir   = file.path(figures_dir, "sim", "flu", scenario)
  )
}

# Shared method colors / labels (matches the originals).
method_palette <- c("MechRt"    = "#2D6CA2",
                    "EpiNow2"   = "#5D3A1E",
                    "estimateR" = "#D55E00",
                    "EpiEstim"  = "#009E73",
                    "rtestim"   = "#CC79A7",
                    "EpiLPS"    = "#882255",
                    "Ground Truth" = "grey20")
